/*
  ZynAddSubFX - a software synthesizer

  ADnoteParameters.cpp - Parameters for ADnote (ADsynth)
  Copyright (C) 2002-2005 Nasca Octavian Paul
  Author: Nasca Octavian Paul

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "ADnoteParameters.h"
#include "EnvelopeParams.h"
#include "LFOParams.h"
#include "SEQParams.h"
#include "../Misc/Time.h"
#include "../Misc/XMLwrapper.h"
#include "../DSP/FFTwrapper.h"
#include "../Synth/OscilGen.h"
#include "../Synth/Resonance.h"
#include "FilterParams.h"

#include <rtosc/ports.h>
#include <rtosc/port-sugar.h>

namespace zyn {

using rtosc::Ports;
using rtosc::RtData;

#define rObject ADnoteVoiceParam

#undef rChangeCb
#define rChangeCb if (obj->time) { obj->last_update_timestamp = obj->time->time(); }
static const Ports voicePorts = {
    //Send Messages To Oscillator Realtime Table
    {"OscilSmp/", rDoc("Primary Oscillator"),
        &OscilGen::ports,
        rBOIL_BEGIN
            if(obj->OscilGn == NULL) return;
        data.obj = obj->OscilGn;
        SNIP
            OscilGen::realtime_ports.dispatch(msg, data);
        if(data.matches == 0)
            data.forward();
        rBOIL_END},
    {"FMSmp/", rDoc("Modulating Oscillator"),
        &OscilGen::ports,
        rBOIL_BEGIN
            if(obj->FmGn == NULL) return;
        data.obj = obj->FmGn;
        SNIP
            OscilGen::realtime_ports.dispatch(msg, data);
        if(data.matches == 0)
            data.forward();
        rBOIL_END},
    rRecurp(FreqLfo, "Frequency LFO"),
    rRecurp(FreqSeq, "Frequency Sequencer"),
    rRecurp(AmpLfo, "Amplitude LFO"),
    rRecurp(AmpSeq, "Amplitude Sequencer"),
    rRecurp(FilterLfo, "Filter LFO"),
    rRecurp(FilterSeq, "Filter Sequencer"),
    rRecurp(FreqEnvelope,   "Frequency Envelope"),
    rRecurp(AmpEnvelope,    "Amplitude Envelope"),
    rRecurp(FilterEnvelope, "Filter Envelope"),
    rRecurp(FMFreqEnvelope, "Modulator Frequency Envelope"),
    rRecurp(FMAmpEnvelope,  "Modulator Amplitude Envelope"),
    rRecurp(VoiceFilter,    "Optional Voice Filter"),

//    rToggle(Enabled,       rShort("enable"), "Voice Enable"),
    rParamI(Unison_size, rShort("size"), rMap(min, 1), rMap(max, 50),
        rDefault(1), "Number of subvoices"),
    rParamZyn(Unison_phase_randomness, rShort("ph.rnd."), rDefault(127),
        "Phase Randomness"),
    rParamZyn(Unison_frequency_spread, rShort("detune"), rDefault(60),
        "Subvoice detune"),
    rParamZyn(Unison_stereo_spread,    rShort("spread"), rDefault(64),
        "Subvoice L/R Separation"),
    rParamZyn(Unison_vibratto,         rShort("vib."), rDefault(64),
        "Subvoice vibratto"),
    rParamZyn(Unison_vibratto_speed,   rShort("speed"), rDefault(64),
        "Subvoice vibratto speed"),
    rOption(Unison_invert_phase,       rShort("inv."),
            rOptions(none, random, 50%, 33%, 25%), rDefault(none),
            "Subvoice Phases"),
    rOption(Type,            rShort("type"), rOptions(Sound,White,Pink,DC),
        rLinear(0,127), rDefault(Sound), "Type of Sound"),
    rParamZyn(PDelay,        rShort("delay"), rDefault(0),
        "Voice Startup Delay"),
    rToggle(Presonance,      rShort("enable"), rDefault(true),
        "Resonance Enable"),
    rParamI(Pextoscil,       rDefault(-1),     rShort("ext."),
            rMap(min, -1), rMap(max, 16), "External Oscillator Selection"),
    rParamI(PextFMoscil,     rDefault(-1),     rShort("ext."),
            rMap(min, -1), rMap(max, 16), "External FM Oscillator Selection"),
    rParamZyn(Poscilphase,   rShort("phase"),  rDefault(64),
        "Oscillator Phase"),
    rParamZyn(PFMoscilphase, rShort("phase"),  rDefault(64),
        "FM Oscillator Phase"),
    rToggle(Pfilterbypass,   rShort("bypass"), rDefault(false),
        "Filter Bypass"),

    //Freq Stuff
    rToggle(Pfixedfreq,      rShort("fixed"),  rDefault(false),
        "If frequency is fixed"),
    rParamZyn(PfixedfreqET,  rShort("e.t."),   rDefault(0),
        "Equal Temperament Parameter"),
    rParamZyn(PBendAdjust,   rShort("bend"),   rDefault(88) /* 64 + 24 */,
        "Pitch bend adjustment"),
    rParamZyn(POffsetHz,     rShort("offset"), rDefault(64),
        "Voice constant offset"),
    //nominally -8192..8191
    rParamI(PDetune,              rShort("fine"),       rLinear(0, 16383),
        rDefault(8192),           "Fine Detune"),
    rParamI(PCoarseDetune,        rShort("coarse"),     rDefault(0),
        "Coarse Detune"),
    rParamZyn(PDetuneType,        rShort("type"),
        rOptions(L35cents, L10cents, E100cents, E1200cents), rDefault(L35cents),
        "Magnitude of Detune"),
    rToggle(PFreqEnvelopeEnabled, rShort("enable"),     rDefault(false),
        "Frequency Envelope Enable"),
    rToggle(PFreqLfoEnabled,      rShort("enable"),     rDefault(false),
        "Frequency LFO Enable"),
    rToggle(PFreqSeqEnabled,      rShort("enable"),     rDefault(false),
        "Frequency Sequencer Enable"),

    //Amplitude Stuff
    rParamZyn(PPanning,                  rShort("pan."), rDefault(64),
        "Panning"),
    {"PVolume::i", rShort("vol.") rLinear(0,127)
        rDoc("Volume"), NULL,
        [](const char *msg, RtData &d)
        {
            rObject *obj = (rObject *)d.obj;
            if (!rtosc_narguments(msg))
                d.reply(d.loc, "i", (int)roundf(127.0f * (1.0f + obj->volume/60.0f)));
            else
                obj->volume = -60.0f * (1.0f - rtosc_argument(msg, 0).i / 127.0f);
        }},
    {"volume::f", rShort("volume") rProp(parameter) rUnit(dB) rDefault(-12.75) rLinear(-60.0f, 0.0f)
        rDoc("Part Volume"), NULL,
        [](const char *msg, RtData &d)
        {
            rObject *obj = (rObject *)d.obj;
            if(!rtosc_narguments(msg)) {
               d.reply(d.loc, "f", obj->volume);
            } else if (rtosc_narguments(msg) == 1 && rtosc_type(msg, 0) == 'f') {
               obj->volume = rtosc_argument(msg, 0).f;
               d.broadcast(d.loc, "f", ((obj->volume)));
            }
        }},
    rToggle(PVolumeminus,                rShort("inv."), rDefault(false),
        "Signal Inverter"), //do we really need this??
    rToggle(PAAEnabled,        rShort("enable"), rDefault(false),
        "AntiAliasing Enable"),
    rParamZyn(PAmpVelocityScaleFunction, rShort("sense"), rDefault(127),
        "Velocity Sensing"),
    rToggle(PAmpEnvelopeEnabled, rShort("enable"), rDefault(false),
        "Amplitude Envelope Enable"),
    rToggle(PAmpLfoEnabled,      rShort("enable"), rDefault(false),
        "Amplitude LFO Enable"),
    rToggle(PAmpSeqEnabled,      rShort("enable"), rDefault(false),
        "Amplitude Sequencer Enable"),

    //Filter Stuff
    rToggle(PFilterEnabled,         rShort("enable"), rDefault(false),
        "Filter Enable"),
    rToggle(PFilterEnvelopeEnabled, rShort("enable"), rDefault(false),
        "Filter Envelope Enable"),
    rToggle(PFilterLfoEnabled,      rShort("enable"), rDefault(false),
        "Filter LFO Enable"),
    rToggle(PFilterSeqEnabled,      rShort("enable"), rDefault(false),
        "Filter Sequencer Enable"),
    rParamZyn(PFilterVelocityScale,         rShort("v.scale"), rDefault(0),
        "Filter Velocity Magnitude"),
    rParamZyn(PFilterVelocityScaleFunction, rShort("v.sense"), rDefault(64),
        "Filter Velocity Function Shape"),


    //Modulator Stuff
    rOption(PFMEnabled, rShort("mode"), rOptions(none, mix, ring, phase,
                frequency, pulse), rLinear(0,127), rDefault(none), "Modulator mode"),
    rParamI(PFMVoice,                   rShort("voice"), rDefault(-1),
        "Modulator Oscillator Selection"),
    rParamF(FMvolume,                   rShort("vol."),  rLinear(0.0, 100.0),
        rDefault(70.0),                 "Modulator Magnitude"),
    rParamZyn(PFMVolumeDamp,            rShort("damp."), rDefault(64),
        "Modulator HF dampening"),
    rParamZyn(PFMVelocityScaleFunction, rShort("sense"), rDefault(64),
        "Modulator Velocity Function"),
    //nominally -8192..8191
    rParamI(PFMDetune,                  rShort("fine"),
            rLinear(0, 16383), rDefault(8192), "Modulator Fine Detune"),
    rParamI(PFMCoarseDetune,            rShort("coarse"), rDefault(0),
            "Modulator Coarse Detune"),
    rParamZyn(PFMDetuneType,            rShort("type"),
              rOptions(L35cents, L10cents, E100cents, E1200cents),
              rDefault(L35cents),
              "Modulator Detune Magnitude"),
    rToggle(PFMFixedFreq,               rShort("fixed"),  rDefault(false),
            "Modulator Frequency Fixed"),
    rToggle(PFMFreqEnvelopeEnabled,  rShort("enable"), rDefault(false),
            "Modulator Frequency Envelope"),
    rToggle(PFMAmpEnvelopeEnabled,   rShort("enable"), rDefault(false),
            "Modulator Amplitude Envelope"),



    //weird stuff for PCoarseDetune
    {"detunevalue:",  rMap(unit,cents) rDoc("Get detune in cents"), NULL,
        [](const char *, RtData &d)
        {
            rObject *obj = (rObject *)d.obj;
            unsigned detuneType =
            obj->PDetuneType == 0 ? *(obj->GlobalPDetuneType)
            : obj->PDetuneType;
            //TODO check if this is accurate or if PCoarseDetune is utilized
            //TODO do the same for the other engines
            d.reply(d.loc, "f", getdetune(detuneType, 0, obj->PDetune));
        }},
    {"octave::c:i", rProp(parameter) rShort("octave") rLinear(-8, 7) rDoc("Octave note offset"),
        NULL,
        [](const char *msg, RtData &d)
        {
            rObject *obj = (rObject *)d.obj;
            auto get_octave = [&obj](){
                int k=obj->PCoarseDetune/1024;
                if (k>=8) k-=16;
                return k;
            };
            if(!rtosc_narguments(msg)) {
                d.reply(d.loc, "i", get_octave());
            } else {
                int k=(int) rtosc_argument(msg, 0).i;
                if (k<0) k+=16;
                obj->PCoarseDetune = k*1024 + obj->PCoarseDetune%1024;
                d.broadcast(d.loc, "i", get_octave());
            }
        }},
    {"coarsedetune::c:i", rProp(parameter) rShort("coarse") rLinear(-64,63) rDefault(0)
        rDoc("Coarse note detune"), NULL,
        [](const char *msg, RtData &d)
        {
            rObject *obj = (rObject *)d.obj;
            auto get_coarse = [&obj](){
                int k=obj->PCoarseDetune%1024;
                if (k>=512) k-=1024;
                return k;
            };
            if(!rtosc_narguments(msg)) {
                d.reply(d.loc, "i", get_coarse());
            } else {
                int k=(int) rtosc_argument(msg, 0).i;
                if (k<0) k+=1024;
                obj->PCoarseDetune = k + (obj->PCoarseDetune/1024)*1024;
                d.broadcast(d.loc, "i", get_coarse());
            }
        }},
    {"PFMVolume::i", rShort("vol.") rLinear(0,127)
        rDoc("Modulator Magnitude"), NULL,
        [](const char *msg, RtData &d)
        {
            rObject *obj = (rObject *)d.obj;
            if (!rtosc_narguments(msg))
                d.reply(d.loc, "i", (int)roundf(127.0f * obj->FMvolume
                    / 100.0f));
            else
                obj->FMvolume = 100.0f * rtosc_argument(msg, 0).i / 127.0f;
        }},
    //weird stuff for PCoarseDetune
    {"FMdetunevalue:", rMap(unit,cents) rDoc("Get modulator detune"), NULL, [](const char *, RtData &d)
        {
            rObject *obj = (rObject *)d.obj;
            unsigned detuneType =
            obj->PFMDetuneType == 0 ? *(obj->GlobalPDetuneType)
            : obj->PFMDetuneType;
            //TODO check if this is accurate or if PCoarseDetune is utilized
            //TODO do the same for the other engines
            d.reply(d.loc, "f", getdetune(detuneType, 0, obj->PFMDetune));
        }},
    {"FMoctave::c:i", rProp(parameter) rShort("octave") rLinear(-8,7) rDoc("Octave note offset for modulator"), NULL,
        [](const char *msg, RtData &d)
        {
            rObject *obj = (rObject *)d.obj;
            auto get_octave = [&obj](){
                int k=obj->PFMCoarseDetune/1024;
                if (k>=8) k-=16;
                return k;
            };
            if(!rtosc_narguments(msg)) {
                d.reply(d.loc, "i", get_octave());
            } else {
                int k=(int) rtosc_argument(msg, 0).i;
                if (k<0) k+=16;
                obj->PFMCoarseDetune = k*1024 + obj->PFMCoarseDetune%1024;
                d.broadcast(d.loc, "i", get_octave());
            }
        }},
    {"FMcoarsedetune::c:i", rProp(parameter) rShort("coarse") rLinear(-64,63)
        rDoc("Coarse note detune for modulator"),
        NULL, [](const char *msg, RtData &d)
        {
            rObject *obj = (rObject *)d.obj;
            auto get_coarse = [&obj](){
                int k=obj->PFMCoarseDetune%1024;
                if (k>=512) k-=1024;
                return k;
            };
            if(!rtosc_narguments(msg)) {
                d.reply(d.loc, "i", get_coarse());
            } else {
                int k=(int) rtosc_argument(msg, 0).i;
                if (k<0) k+=1024;
                obj->PFMCoarseDetune = k + (obj->PFMCoarseDetune/1024)*1024;
                d.broadcast(d.loc, "i", get_coarse());
            }
        }},

    //Reader
    {"unisonFrequencySpreadCents:", rMap(unit,cents) rDoc("Unison Frequency Spread"),
        NULL, [](const char *, RtData &d)
        {
            rObject *obj = (rObject *)d.obj;
            d.reply(d.loc, "f", obj->getUnisonFrequencySpreadCents());
        }},
};
#undef rChangeCb

#undef  rObject
#define rObject ADnoteGlobalParam

#define rChangeCb if (obj->time) { obj->last_update_timestamp = obj->time->time(); }
static const Ports globalPorts = {
    rRecurp(Reson, "Resonance"),
    rRecurp(FreqLfo, "Frequency LFO"),
    rRecurp(FreqSeq, "Frequency Sequencer"),
    rRecurp(AmpLfo, "Amplitude LFO"),
    rRecurp(AmpSeq, "Amplitude Sequencer"),
    rRecurp(FilterLfo, "Filter LFO"),
    rRecurp(FilterSeq, "Filter Sequencer"),
    rRecurp(FreqEnvelope, "Frequency Envelope"),
    rRecurp(AmpEnvelope, "Frequency Envelope"),
    rRecurp(FilterEnvelope, "Frequency Envelope"),
    rRecurp(GlobalFilter, "Filter"),

    rToggle(PStereo, rShort("stereo"), rDefault(true), "Mono/Stereo Enable"),

    //Frequency
    //nominally -8192..8191
    rParamI(PDetune,              rShort("fine"),
            rLinear(0, 16383), rDefault(8192), "Fine Detune"),
    rParamI(PCoarseDetune,   rShort("coarse"), rDefault(0), "Coarse Detune"),
    rParamZyn(PDetuneType,   rShort("type"),
              rOptions(L35cents, L10cents, E100cents, E1200cents),
              rDefault(L10cents),
              "Detune Scaling Type"),
    rParamZyn(PBandwidth,    rShort("bw."), rDefault(64),
              "Relative Fine Detune Gain"),

    //Amplitude
    rParamZyn(PPanning, rShort("pan"), rDefault(64),
        "Panning of ADsynth (0 random, 1 left, 127 right)"),
    rParamF(Volume,                   rShort("vol"), rLinear(-47.9588f,32.0412f),
        rUnit(dB), rDefault(8.29f), "volume control"),
    rParamZyn(PAmpVelocityScaleFunction, rShort("sense"), rDefault(64),
        "Volume velocity sense"),
    {"PVolume::i", rShort("vol.") rLinear(0,127)
        rDoc("Volume"), NULL,
        [](const char *msg, RtData &d)
        {
            rObject *obj = (rObject *)d.obj;
            if (!rtosc_narguments(msg))
                d.reply(d.loc, "i", (int)roundf(96.0f * (1.0f + (obj->Volume - 12.0412f)/60.0f)));
            else
                obj->Volume = 12.0412f - 60.0f * (1.0f - rtosc_argument(msg, 0).i / 96.0f);
        }},
    rParamZyn(Fadein_adjustment, rDefault(FADEIN_ADJUSTMENT_SCALE),
        "Adjustment for anti-pop strategy."),
    rParamZyn(PPunchStrength, rShort("strength"),     rDefault(0),
        "Punch Strength"),
    rParamZyn(PPunchTime,     rShort("time"),         rDefault(60),
        "Length of Punch"),
    rParamZyn(PPunchStretch,  rShort("stretch"),      rDefault(64),
        "How Punch changes with note frequency"),
    rParamZyn(PPunchVelocitySensing, rShort("v.sns"), rDefault(72),
        "Punch Velocity control"),

    //Filter
    rParamZyn(PFilterVelocityScale,         rShort("scale"), rDefault(0),
        "Filter Velocity Magnitude"),
    rParamZyn(PFilterVelocityScaleFunction, rShort("sense"), rDefault(64),
        "Filter Velocity Function Shape"),


    //Resonance
    rToggle(Hrandgrouping, rDefault(false),
        "How randomness is applied to multiple voices using the same oscil"),

    //weird stuff for PCoarseDetune
    {"detunevalue:", rMap(unit,cents) rDoc("Get detune in cents"), NULL,
        [](const char *, RtData &d)
        {
            rObject *obj = (rObject *)d.obj;
            d.reply(d.loc, "f", getdetune(obj->PDetuneType, 0, obj->PDetune));
        }},
    {"octave::c:i", rProp(parameter) rShort("octave") rLinear(-8,7)
        rDoc("Octave note offset"), NULL,
        [](const char *msg, RtData &d)
        {
            rObject *obj = (rObject *)d.obj;
            auto get_octave = [&obj](){
                int k=obj->PCoarseDetune/1024;
                if (k>=8) k-=16;
                return k;
            };
            if(!rtosc_narguments(msg)) {
                d.reply(d.loc, "i", get_octave());
            } else {
                int k=(int) rtosc_argument(msg, 0).i;
                if (k<0) k+=16;
                obj->PCoarseDetune = k*1024 + obj->PCoarseDetune%1024;
                d.broadcast(d.loc, "i", get_octave());
            }
        }},
    {"coarsedetune::c:i", rProp(parameter) rShort("coarse") rLinear(-64, 63)
        rDoc("Coarse note detune"), NULL,
        [](const char *msg, RtData &d)
        {
            rObject *obj = (rObject *)d.obj;
            auto get_coarse = [&obj](){
                int k=obj->PCoarseDetune%1024;
                if (k>=512) k-=1024;
                return k;
            };
            if(!rtosc_narguments(msg)) {
                d.reply(d.loc, "i", get_coarse());
            } else {
                int k=(int) rtosc_argument(msg, 0).i;
                if (k<0) k+=1024;
                obj->PCoarseDetune = k + (obj->PCoarseDetune/1024)*1024;
                d.broadcast(d.loc, "i", get_coarse());
            }
        }},

};
#undef rChangeCb

#undef  rObject
#define rObject ADnoteParameters

#define rChangeCb obj->last_update_timestamp = obj->time->time();
static const Ports adPorts = {//XXX 16 should not be hard coded
    rSelf(ADnoteParameters),
    rPaste,
    rArrayPaste,
    rRecurs(VoicePar, NUM_VOICES),
    {"VoicePar#" STRINGIFY(NUM_VOICES) "/Enabled::T:F",
     rProp(parameter) rShort("enable") rDoc("Voice Enable")
     rDefault([true false false ...]),
     NULL, rArrayTCbMember(VoicePar, Enabled)},
    rRecur(GlobalPar, "Adnote Parameters"),
};
#undef rChangeCb
const Ports &ADnoteParameters::ports  = adPorts;
const Ports &ADnoteVoiceParam::ports  = voicePorts;
const Ports &ADnoteGlobalParam::ports = globalPorts;

ADnoteParameters::ADnoteParameters(const SYNTH_T &synth, FFTwrapper *fft_,
                                   const AbsTime *time_)
    :PresetsArray(), GlobalPar(time_), time(time_), last_update_timestamp(0)
{
    setpresettype("Padsynth");
    fft = fft_;


    for(int nvoice = 0; nvoice < NUM_VOICES; ++nvoice) {
        VoicePar[nvoice].GlobalPDetuneType = &GlobalPar.PDetuneType;
        VoicePar[nvoice].time = time_;
        EnableVoice(synth, nvoice, time_);
    }

    defaults();
}

ADnoteGlobalParam::ADnoteGlobalParam(const AbsTime *time_) :
        time(time_), last_update_timestamp(0)
{
    FreqEnvelope = new EnvelopeParams(0, 0, time_);
    FreqEnvelope->init(ad_global_freq);
    FreqLfo = new LFOParams(ad_global_freq, time_);
    FreqSeq = new SEQParams(ad_global_freq, time_);

    AmpEnvelope = new EnvelopeParams(64, 1, time_);
    AmpEnvelope->init(ad_global_amp);
    AmpLfo = new LFOParams(ad_global_amp, time_);
    AmpSeq = new SEQParams(ad_global_amp, time_);

    GlobalFilter   = new FilterParams(ad_global_filter, time_);
    FilterEnvelope = new EnvelopeParams(0, 1, time_);
    FilterEnvelope->init(ad_global_filter);
    FilterLfo = new LFOParams(ad_global_filter, time_);
    FilterSeq = new SEQParams(ad_global_filter, time_);
    Reson     = new Resonance();
}

void ADnoteParameters::defaults()
{
    //Default Parameters
    GlobalPar.defaults();

    for(int nvoice = 0; nvoice < NUM_VOICES; ++nvoice)
        defaults(nvoice);

    VoicePar[0].Enabled = 1;
}

void ADnoteGlobalParam::defaults()
{
    /* Frequency Global Parameters */
    PStereo = 1; //stereo
    PDetune = 8192; //zero
    PCoarseDetune = 0;
    PDetuneType   = 1;
    FreqEnvelope->defaults();
    FreqLfo->defaults();
    FreqSeq->defaults();
    PBandwidth = 64;

    /* Amplitude Global Parameters */
    Volume  = 8.29f;
    PPanning = 64; //center
    PAmpVelocityScaleFunction = 64;
    AmpEnvelope->defaults();
    AmpLfo->defaults();
    AmpSeq->defaults();
    Fadein_adjustment = FADEIN_ADJUSTMENT_SCALE;
    PPunchStrength = 0;
    PPunchTime     = 60;
    PPunchStretch  = 64;
    PPunchVelocitySensing = 72;
    Hrandgrouping = 0;

    /* Filter Global Parameters*/
    PFilterVelocityScale = 0;
    PFilterVelocityScaleFunction = 64;
    GlobalFilter->defaults();
    FilterEnvelope->defaults();
    FilterLfo->defaults();
    FilterSeq->defaults();
    Reson->defaults();
}

/*
 * Defaults a voice
 */
void ADnoteParameters::defaults(int n)
{
    VoicePar[n].defaults();
}

void ADnoteVoiceParam::defaults()
{
    Enabled = 0;

    Unison_size = 1;
    Unison_frequency_spread = 60;
    Unison_stereo_spread    = 64;
    Unison_vibratto = 64;
    Unison_vibratto_speed = 64;
    Unison_invert_phase   = 0;
    Unison_phase_randomness = 127;

    Type = 0;
    Pfixedfreq    = 0;
    PfixedfreqET  = 0;
    PBendAdjust = 88; // 64 + 24
    POffsetHz     = 64;
    Presonance    = 1;
    Pfilterbypass = 0;
    Pextoscil     = -1;
    PextFMoscil   = -1;
    Poscilphase   = 64;
    PFMoscilphase = 64;
    PDelay                    = 0;
    volume                    = -60.0f* (1.0f - 100.0f / 127.0f);
    PVolumeminus              = 0;
    PAAEnabled                = 0;
    PPanning                  = 64; //center
    PDetune                   = 8192; //8192=0
    PCoarseDetune             = 0;
    PDetuneType               = 0;
    PFreqLfoEnabled           = 0;
    PFreqSeqEnabled           = 0;
    PFreqEnvelopeEnabled      = 0;
    PAmpEnvelopeEnabled       = 0;
    PAmpLfoEnabled            = 0;
    PAmpSeqEnabled            = 0;
    PAmpVelocityScaleFunction = 127;
    PFilterEnabled            = 0;
    PFilterEnvelopeEnabled    = 0;
    PFilterLfoEnabled         = 0;
    PFilterSeqEnabled         = 0;
    PFilterVelocityScale = 0;
    PFilterVelocityScaleFunction = 64;
    PFMEnabled                = FMTYPE::NONE;
    PFMFixedFreq              = false;

    //I use the internal oscillator (-1)
    PFMVoice = -1;

    FMvolume       = 70.0;
    PFMVolumeDamp   = 64;
    PFMDetune       = 8192;
    PFMCoarseDetune = 0;
    PFMDetuneType   = 0;
    PFMFreqEnvelopeEnabled   = 0;
    PFMAmpEnvelopeEnabled    = 0;
    PFMVelocityScaleFunction = 64;

    OscilGn->defaults();
    FmGn->defaults();

    AmpEnvelope->defaults();
    AmpLfo->defaults();
    AmpSeq->defaults();

    FreqEnvelope->defaults();
    FreqLfo->defaults();
    FreqSeq->defaults();

    VoiceFilter->defaults();
    FilterEnvelope->defaults();
    FilterLfo->defaults();
    FilterSeq->defaults();

    FMFreqEnvelope->defaults();
    FMAmpEnvelope->defaults();
}



/*
 * Init the voice parameters
 */
void ADnoteParameters::EnableVoice(const SYNTH_T &synth, int nvoice,
                                   const AbsTime *time)
{
    VoicePar[nvoice].enable(synth, fft, GlobalPar.Reson, time);
}

void ADnoteVoiceParam::enable(const SYNTH_T &synth, FFTwrapper *fft,
                              Resonance *Reson, const AbsTime *time)
{
    OscilGn  = new OscilGen(synth, fft, Reson);
    FmGn    = new OscilGen(synth, fft, NULL);

    AmpEnvelope = new EnvelopeParams(64, 1, time);
    AmpEnvelope->init(ad_voice_amp);
    AmpLfo = new LFOParams(ad_voice_amp, time);
    AmpSeq = new SEQParams(ad_voice_amp, time);

    FreqEnvelope = new EnvelopeParams(0, 0, time);
    FreqEnvelope->init(ad_voice_freq);
    FreqLfo = new LFOParams(ad_voice_freq, time);
    FreqSeq = new SEQParams(ad_voice_freq, time);

    VoiceFilter    = new FilterParams(ad_voice_filter, time);
    FilterEnvelope = new EnvelopeParams(0, 0, time);
    FilterEnvelope->init(ad_voice_filter);
    FilterLfo = new LFOParams(ad_voice_filter, time);
    FilterSeq = new SEQParams(ad_voice_filter, time);

    FMFreqEnvelope = new EnvelopeParams(0, 0, time);
    FMFreqEnvelope->init(ad_voice_fm_freq);
    FMAmpEnvelope = new EnvelopeParams(64, 1, time);
    FMAmpEnvelope->init(ad_voice_fm_amp);
}

/*
 * Get the Multiplier of the fine detunes of the voices
 */
float ADnoteParameters::getBandwidthDetuneMultiplier() const
{
    float bw = (GlobalPar.PBandwidth - 64.0f) / 64.0f;
    bw = powf(2.0f, bw * powf(fabsf(bw), 0.2f) * 5.0f);

    return bw;
}

/*
 * Get the unison spread in cents for a voice
 */

float ADnoteParameters::getUnisonFrequencySpreadCents(int nvoice) const
{
    return VoicePar[nvoice].getUnisonFrequencySpreadCents();
}

float ADnoteVoiceParam::getUnisonFrequencySpreadCents(void) const {
    return powf(Unison_frequency_spread / 127.0f * 2.0f, 2.0f) * 50.0f; //cents
}

/*
 * Kill the voice
 */
void ADnoteParameters::KillVoice(int nvoice)
{
    VoicePar[nvoice].kill();
}

void ADnoteVoiceParam::kill()
{
    delete OscilGn;
    delete FmGn;

    delete AmpEnvelope;
    delete AmpLfo;
    delete AmpSeq;

    delete FreqEnvelope;
    delete FreqLfo;
    delete FreqSeq;

    delete VoiceFilter;
    delete FilterEnvelope;
    delete FilterLfo;
    delete FilterSeq;

    delete FMFreqEnvelope;
    delete FMAmpEnvelope;
}


ADnoteGlobalParam::~ADnoteGlobalParam()
{
    delete FreqEnvelope;
    delete FreqLfo;
    delete FreqSeq;
    delete AmpEnvelope;
    delete AmpLfo;
    delete AmpSeq;
    delete GlobalFilter;
    delete FilterEnvelope;
    delete FilterLfo;
    delete FilterSeq;
    delete Reson;
}

ADnoteParameters::~ADnoteParameters()
{
    for(int nvoice = 0; nvoice < NUM_VOICES; ++nvoice)
        KillVoice(nvoice);
}

void ADnoteParameters::add2XMLsection(XMLwrapper& xml, int n)
{
    int nvoice = n;
    if(nvoice >= NUM_VOICES)
        return;

    int oscilused = 0, fmoscilused = 0; //if the oscil or fmoscil are used by another voice

    for(int i = 0; i < NUM_VOICES; ++i) {
        if(VoicePar[i].Pextoscil == nvoice)
            oscilused = 1;
        if(VoicePar[i].PextFMoscil == nvoice)
            fmoscilused = 1;
    }

    xml.addparbool("enabled", VoicePar[nvoice].Enabled);
    if(((VoicePar[nvoice].Enabled == 0) && (oscilused == 0)
        && (fmoscilused == 0)) && (xml.minimal))
        return;

    VoicePar[nvoice].add2XML(xml, fmoscilused);
}

void ADnoteVoiceParam::add2XML(XMLwrapper& xml, bool fmoscilused)
{
    xml.addpar("type", Type);

    xml.addpar("unison_size", Unison_size);
    xml.addpar("unison_frequency_spread",
                Unison_frequency_spread);
    xml.addpar("unison_stereo_spread", Unison_stereo_spread);
    xml.addpar("unison_vibratto", Unison_vibratto);
    xml.addpar("unison_vibratto_speed", Unison_vibratto_speed);
    xml.addpar("unison_invert_phase", Unison_invert_phase);
    xml.addpar("unison_phase_randomness", Unison_phase_randomness);

    xml.addpar("delay", PDelay);
    xml.addparbool("resonance", Presonance);

    xml.addpar("ext_oscil", Pextoscil);
    xml.addpar("ext_fm_oscil", PextFMoscil);

    xml.addpar("oscil_phase", Poscilphase);
    xml.addpar("oscil_fm_phase", PFMoscilphase);

    xml.addparbool("filter_enabled", PFilterEnabled);
    xml.addparbool("filter_bypass", Pfilterbypass);

    xml.addpar("fm_enabled", (int)PFMEnabled);

    xml.beginbranch("OSCIL");
    OscilGn->add2XML(xml);
    xml.endbranch();


    xml.beginbranch("AMPLITUDE_PARAMETERS");
    xml.addpar("panning", PPanning);
    xml.addparreal("volume", volume);
    xml.addparbool("volume_minus", PVolumeminus);
    xml.addpar("velocity_sensing", PAmpVelocityScaleFunction);

    xml.addparbool("amp_envelope_enabled",
                    PAmpEnvelopeEnabled);
    if((PAmpEnvelopeEnabled != 0) || (!xml.minimal)) {
        xml.beginbranch("AMPLITUDE_ENVELOPE");
        AmpEnvelope->add2XML(xml);
        xml.endbranch();
    }
    xml.addparbool("amp_lfo_enabled", PAmpLfoEnabled);
    if((PAmpLfoEnabled != 0) || (!xml.minimal)) {
        xml.beginbranch("AMPLITUDE_LFO");
        AmpLfo->add2XML(xml);
        xml.endbranch();
    }
    xml.addparbool("amp_seq_enabled", PAmpSeqEnabled);
    if((PAmpSeqEnabled != 0) || (!xml.minimal)) {
        xml.beginbranch("AMPLITUDE_SEQ");
        AmpSeq->add2XML(xml);
        xml.endbranch();
    }
    xml.endbranch(); // AMPLITUDE_PARAMETERS

    xml.beginbranch("FREQUENCY_PARAMETERS");
    xml.addparbool("fixed_freq", Pfixedfreq);
    xml.addpar("fixed_freq_et", PfixedfreqET);
    xml.addpar("bend_adjust", PBendAdjust);
    xml.addpar("offset_hz", POffsetHz);
    xml.addpar("detune", PDetune);
    xml.addpar("coarse_detune", PCoarseDetune);
    xml.addpar("detune_type", PDetuneType);

    xml.addparbool("freq_envelope_enabled",
                    PFreqEnvelopeEnabled);
    if((PFreqEnvelopeEnabled != 0) || (!xml.minimal)) {
        xml.beginbranch("FREQUENCY_ENVELOPE");
        FreqEnvelope->add2XML(xml);
        xml.endbranch();
    }
    xml.addparbool("freq_lfo_enabled", PFreqLfoEnabled);
    if((PFreqLfoEnabled != 0) || (!xml.minimal)) {
        xml.beginbranch("FREQUENCY_LFO");
        FreqLfo->add2XML(xml);
        xml.endbranch();
    }
    xml.addparbool("freq_seq_enabled", PFreqSeqEnabled);
    if((PFreqSeqEnabled != 0) || (!xml.minimal)) {
        xml.beginbranch("FREQUENCY_SEQ");
        FreqSeq->add2XML(xml);
        xml.endbranch();
    }
    xml.endbranch(); // FREQUENCY_PARAMETERS


    if((PFilterEnabled != 0) || (!xml.minimal)) {
        xml.beginbranch("FILTER_PARAMETERS");
        xml.addpar("velocity_sensing_amplitude", PFilterVelocityScale);
        xml.addpar("velocity_sensing", PFilterVelocityScaleFunction);
        xml.beginbranch("FILTER");
        VoiceFilter->add2XML(xml);
        xml.endbranch();

        xml.addparbool("filter_envelope_enabled",
                        PFilterEnvelopeEnabled);
        if((PFilterEnvelopeEnabled != 0) || (!xml.minimal)) {
            xml.beginbranch("FILTER_ENVELOPE");
            FilterEnvelope->add2XML(xml);
            xml.endbranch();
        }

        xml.addparbool("filter_lfo_enabled",
                        PFilterLfoEnabled);
        if((PFilterLfoEnabled != 0) || (!xml.minimal)) {
            xml.beginbranch("FILTER_LFO");
            FilterLfo->add2XML(xml);
            xml.endbranch();
        }

        xml.addparbool("filter_seq_enabled",
                        PFilterSeqEnabled);
        if((PFilterSeqEnabled != 0) || (!xml.minimal)) {
            xml.beginbranch("FILTER_SEQ");
            FilterSeq->add2XML(xml);
            xml.endbranch();
        }
        xml.endbranch(); // FILTER_PARAMETERS
    }

    if((PFMEnabled != FMTYPE::NONE) || (fmoscilused != 0)
       || (!xml.minimal)) {
        xml.beginbranch("FM_PARAMETERS");
        xml.addpar("input_voice", PFMVoice);

        xml.addparreal("volume", FMvolume);
        xml.addpar("volume_damp", PFMVolumeDamp);
        xml.addpar("velocity_sensing",
                    PFMVelocityScaleFunction);

        xml.addparbool("amp_envelope_enabled",
                        PFMAmpEnvelopeEnabled);
        if((PFMAmpEnvelopeEnabled != 0) || (!xml.minimal)) {
            xml.beginbranch("AMPLITUDE_ENVELOPE");
            FMAmpEnvelope->add2XML(xml);
            xml.endbranch();
        }
        xml.beginbranch("MODULATOR");
        xml.addpar("detune", PFMDetune);
        xml.addpar("coarse_detune", PFMCoarseDetune);
        xml.addpar("detune_type", PFMDetuneType);

        xml.addparbool("freq_envelope_enabled",
                        PFMFreqEnvelopeEnabled);
        xml.addparbool("fixed_freq", PFMFixedFreq);
        if((PFMFreqEnvelopeEnabled != 0) || (!xml.minimal)) {
            xml.beginbranch("FREQUENCY_ENVELOPE");
            FMFreqEnvelope->add2XML(xml);
            xml.endbranch();
        }

        xml.beginbranch("OSCIL");
        FmGn->add2XML(xml);
        xml.endbranch();

        xml.endbranch(); // MODULATOR
        xml.endbranch(); // FM_PARAMETERS
    }
}

void ADnoteGlobalParam::add2XML(XMLwrapper& xml)
{
    xml.addparbool("stereo", PStereo);

    xml.beginbranch("AMPLITUDE_PARAMETERS");
    xml.addparreal("volume", Volume);
    xml.addpar("panning", PPanning);
    xml.addpar("velocity_sensing", PAmpVelocityScaleFunction);
    xml.addpar("fadein_adjustment", Fadein_adjustment);
    xml.addpar("punch_strength", PPunchStrength);
    xml.addpar("punch_time", PPunchTime);
    xml.addpar("punch_stretch", PPunchStretch);
    xml.addpar("punch_velocity_sensing", PPunchVelocitySensing);
    xml.addpar("harmonic_randomness_grouping", Hrandgrouping);

    xml.beginbranch("AMPLITUDE_ENVELOPE");
    AmpEnvelope->add2XML(xml);
    xml.endbranch();

    xml.beginbranch("AMPLITUDE_LFO");
    AmpLfo->add2XML(xml);
    xml.endbranch();

    xml.beginbranch("AMPLITUDE_SEQ");
    AmpSeq->add2XML(xml);
    xml.endbranch();
    
    xml.endbranch(); // AMPLITUDE_PARAMETERS

    xml.beginbranch("FREQUENCY_PARAMETERS");
    xml.addpar("detune", PDetune);

    xml.addpar("coarse_detune", PCoarseDetune);
    xml.addpar("detune_type", PDetuneType);

    xml.addpar("bandwidth", PBandwidth);

    xml.beginbranch("FREQUENCY_ENVELOPE");
    FreqEnvelope->add2XML(xml);
    xml.endbranch();

    xml.beginbranch("FREQUENCY_LFO");
    FreqLfo->add2XML(xml);
    xml.endbranch();

    xml.beginbranch("FREQUENCY_SEQ");
    FreqSeq->add2XML(xml);
    xml.endbranch();
    
    xml.endbranch(); // FREQUENCY_PARAMETERS


    xml.beginbranch("FILTER_PARAMETERS");
    xml.addpar("velocity_sensing_amplitude", PFilterVelocityScale);
    xml.addpar("velocity_sensing", PFilterVelocityScaleFunction);

    xml.beginbranch("FILTER");
    GlobalFilter->add2XML(xml);
    xml.endbranch();

    xml.beginbranch("FILTER_ENVELOPE");
    FilterEnvelope->add2XML(xml);
    xml.endbranch();

    xml.beginbranch("FILTER_LFO");
    FilterLfo->add2XML(xml);
    xml.endbranch();

    xml.beginbranch("FILTER_SEQ");
    FilterSeq->add2XML(xml);
    xml.endbranch();
    
    xml.endbranch(); // FILTER_PARAMETERS

    xml.beginbranch("RESONANCE");
    Reson->add2XML(xml);
    xml.endbranch();
}

void ADnoteParameters::add2XML(XMLwrapper& xml)
{
    GlobalPar.add2XML(xml);
    for(int nvoice = 0; nvoice < NUM_VOICES; ++nvoice) {
        xml.beginbranch("VOICE", nvoice);
        add2XMLsection(xml, nvoice);
        xml.endbranch();
    }
}


void ADnoteGlobalParam::getfromXML(XMLwrapper& xml)
{
    PStereo = xml.getparbool("stereo", PStereo);

    if(xml.enterbranch("AMPLITUDE_PARAMETERS")) {
        const bool upgrade_3_0_5 = (xml.fileversion() < version_type(3,0,5));
        const bool upgrade_3_0_3 = (xml.fileversion() < version_type(3,0,3)) ||
            (!xml.hasparreal("volume"));

        if (upgrade_3_0_3) {
            int vol = xml.getpar127("volume", 0);
            Volume = 12.0412f - 60.0f * ( 1.0f - vol / 96.0f);
        } else if (upgrade_3_0_5) {
            printf("file version less than 3.0.5\n");
            Volume = 12.0412f + xml.getparreal("volume", Volume);
        } else {
            Volume = xml.getparreal("volume", Volume);
        }
        PPanning = xml.getpar127("panning", PPanning);
        PAmpVelocityScaleFunction = xml.getpar127("velocity_sensing",
                                                   PAmpVelocityScaleFunction);

        Fadein_adjustment = xml.getpar127("fadein_adjustment", Fadein_adjustment);
        PPunchStrength = xml.getpar127("punch_strength", PPunchStrength);
        PPunchTime     = xml.getpar127("punch_time", PPunchTime);
        PPunchStretch  = xml.getpar127("punch_stretch", PPunchStretch);
        PPunchVelocitySensing = xml.getpar127("punch_velocity_sensing",
                                               PPunchVelocitySensing);
        Hrandgrouping = xml.getpar127("harmonic_randomness_grouping",
                                       Hrandgrouping);

        if(xml.enterbranch("AMPLITUDE_ENVELOPE")) {
            AmpEnvelope->getfromXML(xml);
            xml.exitbranch();
        }

        if(xml.enterbranch("AMPLITUDE_LFO")) {
            AmpLfo->getfromXML(xml);
            xml.exitbranch();
        }

        if(xml.enterbranch("AMPLITUDE_SEQ")) {
            AmpSeq->getfromXML(xml);
            xml.exitbranch();
        }

        xml.exitbranch(); // AMPLITUDE_PARAMETERS
    }

    if(xml.enterbranch("FREQUENCY_PARAMETERS")) {
        PDetune = xml.getpar("detune", PDetune, 0, 16383);
        PCoarseDetune = xml.getpar("coarse_detune", PCoarseDetune, 0, 16383);
        PDetuneType   = xml.getpar127("detune_type", PDetuneType);
        PBandwidth    = xml.getpar127("bandwidth", PBandwidth);

        xml.enterbranch("FREQUENCY_ENVELOPE");
        FreqEnvelope->getfromXML(xml);
        xml.exitbranch();

        xml.enterbranch("FREQUENCY_LFO");
        FreqLfo->getfromXML(xml);
        xml.exitbranch();

        if(xml.enterbranch("FREQUENCY_SEQ")) {
            FreqSeq->getfromXML(xml);
            xml.exitbranch();
        }

        xml.exitbranch(); // FREQUENCY_PARAMETERS
    }


    if(xml.enterbranch("FILTER_PARAMETERS")) {
        PFilterVelocityScale = xml.getpar127("velocity_sensing_amplitude",
                                              PFilterVelocityScale);
        PFilterVelocityScaleFunction = xml.getpar127(
            "velocity_sensing",
            PFilterVelocityScaleFunction);

        xml.enterbranch("FILTER");
        GlobalFilter->getfromXML(xml);
        xml.exitbranch();

        xml.enterbranch("FILTER_ENVELOPE");
        FilterEnvelope->getfromXML(xml);
        xml.exitbranch();

        xml.enterbranch("FILTER_LFO");
        FilterLfo->getfromXML(xml);
        xml.exitbranch();

        if(xml.enterbranch("FILTER_SEQ")) {
            FilterSeq->getfromXML(xml);
            xml.exitbranch();
        }
        
        xml.exitbranch(); // FILTER_PARAMETERS
    }

    if(xml.enterbranch("RESONANCE")) {
        Reson->getfromXML(xml);
        xml.exitbranch();
    }
}

void ADnoteParameters::getfromXML(XMLwrapper& xml)
{
    GlobalPar.getfromXML(xml);

    for(int nvoice = 0; nvoice < NUM_VOICES; ++nvoice) {
        VoicePar[nvoice].Enabled = 0;
        if(xml.enterbranch("VOICE", nvoice) == 0)
            continue;
        getfromXMLsection(xml, nvoice);
        xml.exitbranch();
    }
}

void ADnoteParameters::getfromXMLsection(XMLwrapper& xml, int n)
{
    int nvoice = n;
    if(nvoice >= NUM_VOICES)
        return;

    VoicePar[nvoice].getfromXML(xml, nvoice);
}

void ADnoteParameters::paste(ADnoteParameters &a)
{
    this->GlobalPar.paste(a.GlobalPar);
    for(int i=0; i<NUM_VOICES; ++i)
        this->VoicePar[i].paste(a.VoicePar[i]);

    if ( time ) {
        last_update_timestamp = time->time();
    }
}

void ADnoteParameters::pasteArray(ADnoteParameters &a, int nvoice)
{
    if(nvoice >= NUM_VOICES)
        return;

    VoicePar[nvoice].paste(a.VoicePar[nvoice]);

    if ( time ) {
        last_update_timestamp = time->time();
    }
}

#define copy(x) this->x = a.x
#define RCopy(x) this->x->paste(*a.x)
void ADnoteVoiceParam::paste(ADnoteVoiceParam &a)
{
    //Come on C++ get some darn reflection, this is horrible

    copy(Enabled);
    copy(Unison_size);
    copy(Unison_frequency_spread);
    copy(Unison_stereo_spread);
    copy(Unison_vibratto);
    copy(Unison_vibratto_speed);
    copy(Unison_invert_phase);
    copy(Unison_phase_randomness);
    copy(Type);
    copy(PDelay);
    copy(Presonance);
    copy(Pextoscil);
    copy(PextFMoscil);
    copy(Poscilphase);
    copy(PFMoscilphase);
    copy(PFilterEnabled);
    copy(Pfilterbypass);
    copy(PFMEnabled);
    copy(PFMFixedFreq);

    RCopy(OscilGn);


    copy(PPanning);
    copy(volume);
    copy(PVolumeminus);
    copy(PAmpVelocityScaleFunction);
    copy(PAmpEnvelopeEnabled);

    RCopy(AmpEnvelope);

    copy(PAmpLfoEnabled);

    RCopy(AmpLfo);

    copy(PAmpSeqEnabled);

    RCopy(AmpSeq);

    copy(Pfixedfreq);
    copy(PfixedfreqET);
    copy(PDetune);
    copy(PCoarseDetune);
    copy(PDetuneType);
    copy(PBendAdjust);
    copy(POffsetHz);
    copy(PFreqEnvelopeEnabled);

    RCopy(FreqEnvelope);

    copy(PFreqLfoEnabled);

    RCopy(FreqLfo);
    
    copy(PFreqSeqEnabled);

    RCopy(FreqSeq);

    RCopy(VoiceFilter);

    copy(PFilterEnvelopeEnabled);

    RCopy(FilterEnvelope);

    copy(PFilterLfoEnabled);
    copy(PFilterVelocityScale);
    copy(PFilterVelocityScaleFunction);

    RCopy(FilterLfo);
    
    copy(PFilterSeqEnabled);

    RCopy(FilterSeq);

    copy(PFMVoice);
    copy(FMvolume);
    copy(PFMVolumeDamp);
    copy(PFMVelocityScaleFunction);

    copy(PFMAmpEnvelopeEnabled);

    RCopy(FMAmpEnvelope);

    copy(PFMDetune);
    copy(PFMCoarseDetune);
    copy(PFMDetuneType);
    copy(PFMFreqEnvelopeEnabled);


    RCopy(FMFreqEnvelope);

    RCopy(FmGn);

    if ( time ) {
        last_update_timestamp = time->time();
    }
}

void ADnoteGlobalParam::paste(ADnoteGlobalParam &a)
{
    copy(PStereo);

    copy(Volume);
    copy(PPanning);
    copy(PAmpVelocityScaleFunction);

    copy(Fadein_adjustment);
    copy(PPunchStrength);
    copy(PPunchTime);
    copy(PPunchStretch);
    copy(PPunchVelocitySensing);
    copy(Hrandgrouping);

    RCopy(AmpEnvelope);
    RCopy(AmpLfo);
    RCopy(AmpSeq);

    copy(PDetune);
    copy(PCoarseDetune);
    copy(PDetuneType);
    copy(PBandwidth);

    RCopy(FreqEnvelope);
    RCopy(FreqLfo);
    RCopy(FreqSeq);

    copy(PFilterVelocityScale);
    copy(PFilterVelocityScaleFunction);

    RCopy(GlobalFilter);
    RCopy(FilterEnvelope);
    RCopy(FilterLfo);
    RCopy(FilterSeq);
    RCopy(Reson);

    if ( time ) {
        last_update_timestamp = time->time();
    }
}
#undef copy
#undef RCopy

void ADnoteVoiceParam::getfromXML(XMLwrapper& xml, unsigned nvoice)
{
    Enabled     = xml.getparbool("enabled", 0);
    Unison_size = xml.getpar127("unison_size", Unison_size);
    Unison_frequency_spread = xml.getpar127("unison_frequency_spread",
                                             Unison_frequency_spread);
    Unison_stereo_spread = xml.getpar127("unison_stereo_spread",
                                          Unison_stereo_spread);
    Unison_vibratto = xml.getpar127("unison_vibratto", Unison_vibratto);
    Unison_vibratto_speed = xml.getpar127("unison_vibratto_speed",
                                           Unison_vibratto_speed);
    Unison_invert_phase = xml.getpar127("unison_invert_phase",
                                         Unison_invert_phase);
    Unison_phase_randomness = xml.getpar127("unison_phase_randomness",
                        Unison_phase_randomness);

    Type       = xml.getpar127("type", Type);
    PDelay     = xml.getpar127("delay", PDelay);
    Presonance = xml.getparbool("resonance", Presonance);

    Pextoscil   = xml.getpar("ext_oscil", -1, -1, nvoice - 1);
    PextFMoscil = xml.getpar("ext_fm_oscil", -1, -1, nvoice - 1);

    Poscilphase    = xml.getpar127("oscil_phase", Poscilphase);
    PFMoscilphase  = xml.getpar127("oscil_fm_phase", PFMoscilphase);
    PFilterEnabled = xml.getparbool("filter_enabled", PFilterEnabled);
    Pfilterbypass  = xml.getparbool("filter_bypass", Pfilterbypass);
    PFMEnabled     = (FMTYPE)xml.getpar127("fm_enabled", (int)PFMEnabled);

    if(xml.enterbranch("OSCIL")) {
        OscilGn->getfromXML(xml);
        xml.exitbranch();
    }


    if(xml.enterbranch("AMPLITUDE_PARAMETERS")) {
        PPanning     = xml.getpar127("panning", PPanning);
        const bool upgrade_3_0_3 = (xml.fileversion() < version_type(3,0,3)) ||
            (!xml.hasparreal("volume"));

        if (upgrade_3_0_3) {
            int vol = xml.getpar127("volume", 0);
            volume    = -60.0f * ( 1.0f - vol / 127.0f);
        } else {
            volume    = xml.getparreal("volume", volume);
        }
        PVolumeminus = xml.getparbool("volume_minus", PVolumeminus);
        PAmpVelocityScaleFunction = xml.getpar127("velocity_sensing",
                                                   PAmpVelocityScaleFunction);

        PAmpEnvelopeEnabled = xml.getparbool("amp_envelope_enabled",
                                              PAmpEnvelopeEnabled);
        if(xml.enterbranch("AMPLITUDE_ENVELOPE")) {
            AmpEnvelope->getfromXML(xml);
            xml.exitbranch();
        }

        PAmpLfoEnabled = xml.getparbool("amp_lfo_enabled", PAmpLfoEnabled);
        if(xml.enterbranch("AMPLITUDE_LFO")) {
            AmpLfo->getfromXML(xml);
            xml.exitbranch();
        }

        PAmpSeqEnabled = xml.getparbool("amp_seq_enabled", PAmpSeqEnabled);
        if(xml.enterbranch("AMPLITUDE_SEQ")) {
            AmpSeq->getfromXML(xml);
            xml.exitbranch();
        }
        xml.exitbranch(); // AMPLITUDE_PARAMETERS
    }

    if(xml.enterbranch("FREQUENCY_PARAMETERS")) {
        Pfixedfreq    = xml.getparbool("fixed_freq", Pfixedfreq);
        PfixedfreqET  = xml.getpar127("fixed_freq_et", PfixedfreqET);
        PBendAdjust  = xml.getpar127("bend_adjust", PBendAdjust);
        POffsetHz  = xml.getpar127("offset_hz", POffsetHz);
        PDetune       = xml.getpar("detune", PDetune, 0, 16383);
        PCoarseDetune = xml.getpar("coarse_detune", PCoarseDetune, 0, 16383);
        PDetuneType   = xml.getpar127("detune_type", PDetuneType);
        PFreqEnvelopeEnabled = xml.getparbool("freq_envelope_enabled",
                                               PFreqEnvelopeEnabled);

        if(xml.enterbranch("FREQUENCY_ENVELOPE")) {
            FreqEnvelope->getfromXML(xml);
            xml.exitbranch();
        }

        PFreqLfoEnabled = xml.getparbool("freq_lfo_enabled", PFreqLfoEnabled);

        if(xml.enterbranch("FREQUENCY_LFO")) {
            FreqLfo->getfromXML(xml);
            xml.exitbranch();
        }
        
        PFreqSeqEnabled = xml.getparbool("freq_seq_enabled", PFreqSeqEnabled);
        if(xml.enterbranch("FREQUENCY_SEQ")) {
            FreqSeq->getfromXML(xml);
            xml.exitbranch();
        }
        
        xml.exitbranch(); // FREQUENCY_PARAMETERS
    }

    if(xml.enterbranch("FILTER_PARAMETERS")) {
        PFilterVelocityScale = xml.getpar127("velocity_sensing_amplitude",
                                              PFilterVelocityScale);
        PFilterVelocityScaleFunction = xml.getpar127(
            "velocity_sensing",
            PFilterVelocityScaleFunction);
        if(xml.enterbranch("FILTER")) {
            VoiceFilter->getfromXML(xml);
            xml.exitbranch();
        }

        PFilterEnvelopeEnabled = xml.getparbool("filter_envelope_enabled",
                                                 PFilterEnvelopeEnabled);
        if(xml.enterbranch("FILTER_ENVELOPE")) {
            FilterEnvelope->getfromXML(xml);
            xml.exitbranch();
        }

        PFilterLfoEnabled = xml.getparbool("filter_lfo_enabled",
                                            PFilterLfoEnabled);
        if(xml.enterbranch("FILTER_LFO")) {
            FilterLfo->getfromXML(xml);
            xml.exitbranch();
        }

        PFilterSeqEnabled = xml.getparbool("filter_seq_enabled",
                                            PFilterSeqEnabled);
        if(xml.enterbranch("FILTER_SEQ")) {
            FilterSeq->getfromXML(xml);
            xml.exitbranch();
        }
        xml.exitbranch(); // FILTER_PARAMETERS
    }

    if(xml.enterbranch("FM_PARAMETERS")) {
        const bool upgrade_3_0_3 = (xml.fileversion() < version_type(3,0,3)) ||
            (xml.getparreal("volume", -1) < 0);

        PFMVoice      = xml.getpar("input_voice", PFMVoice, -1, nvoice - 1);
        if (upgrade_3_0_3) {
            int Pvolume = xml.getpar127("volume", 0);
            FMvolume    = 100.0f * Pvolume / 127.0f;
        } else {
            FMvolume    = xml.getparreal("volume", FMvolume);
        }
        PFMVolumeDamp = xml.getpar127("volume_damp", PFMVolumeDamp);
        PFMVelocityScaleFunction = xml.getpar127("velocity_sensing",
                                                  PFMVelocityScaleFunction);

        PFMAmpEnvelopeEnabled = xml.getparbool("amp_envelope_enabled",
                                                PFMAmpEnvelopeEnabled);
        if(xml.enterbranch("AMPLITUDE_ENVELOPE")) {
            FMAmpEnvelope->getfromXML(xml);
            xml.exitbranch();
        }

        if(xml.enterbranch("MODULATOR")) {
            PFMDetune = xml.getpar("detune", PFMDetune, 0, 16383);
            PFMCoarseDetune = xml.getpar("coarse_detune",
                                          PFMCoarseDetune,
                                          0,
                                          16383);
            PFMDetuneType = xml.getpar127("detune_type", PFMDetuneType);

            PFMFreqEnvelopeEnabled = xml.getparbool("freq_envelope_enabled",
                                                     PFMFreqEnvelopeEnabled);
            PFMFixedFreq = xml.getparbool("fixed_freq",
                                                     PFMFixedFreq);
            if(xml.enterbranch("FREQUENCY_ENVELOPE")) {
                FMFreqEnvelope->getfromXML(xml);
                xml.exitbranch();
            }

            if(xml.enterbranch("OSCIL")) {
                FmGn->getfromXML(xml);
                xml.exitbranch();
            }

            xml.exitbranch(); // MODULATOR
        }
        xml.exitbranch(); // FM_PARAMETERS
    }
}

}
