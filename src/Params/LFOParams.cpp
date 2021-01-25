/*
  ZynAddSubFX - a software synthesizer

  LFOParams.cpp - Parameters for LFO
  Copyright (C) 2002-2005 Nasca Octavian Paul
  Author: Nasca Octavian Paul

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.
*/

#include <cmath>
#include <cstdio>
#include "../globals.h"
#include "../Misc/Util.h"
#include "../Misc/XMLwrapper.h"
#include "../Misc/Time.h"
#include "LFOParams.h"

#include <rtosc/port-sugar.h>
#include <rtosc/ports.h>
using namespace rtosc;

namespace zyn {

#define rObject LFOParams
#undef rChangeCb
#define rChangeCb if (obj->time) { obj->last_update_timestamp = obj->time->time(); }
#define rBegin [](const char *msg, rtosc::RtData &d) {
#define rEnd }
static const rtosc::Ports _ports = {
    rSelf(LFOParams),
    rPaste,
    rOption(loc, rProp(internal),
            rOptions(ad_global_amp, ad_global_freq, ad_global_filter,
                     ad_voice_amp, ad_voice_freq, ad_voice_filter, unspecified),
            "location of the filter"),
    rParamF(freq, rShort("freq"), rUnit(HZ), rLog(0.0775679,85.25),
            rDefaultDepends(loc),
            rPreset(ad_global_amp, 6.49), // 80
            rPreset(ad_global_freq, 3.71), // 70
            rPreset(ad_global_filter, 6.49),
            rPreset(ad_voice_amp, 11.25), // 90
            rPreset(ad_voice_freq, 1.19), // 50
            rPreset(ad_voice_filter, 1.19),
            "frequency of LFO\n"
            "lfo frequency = Pfreq * stretch\n"
            "true frequency is [0,85.25] Hz"),
    {"Pfreq::f", rShort("freq.") rLinear(0, 1.0) rDoc("frequency of LFO "
     "lfo frequency = Pfreq * stretch "
     "true frequency is [0,85.25] Hz"), NULL,
     [](const char *msg, RtData &d)
     {
         rObject *obj = (rObject *)d.obj;
         if (!rtosc_narguments(msg)) {
             d.reply(d.loc, "f", log2f(12.0f * obj->freq + 1.0f) / 10.0f);
         } else {
             obj->freq = (powf(2, 10.0f * rtosc_argument(msg, 0).f) - 1.0f) / 12.0f;
         }
     }},
    rParamZyn(Pintensity, rShort("depth"),
              rDefaultDepends(loc),
              rDefault(0), rPreset(ad_voice_amp, 32),
              rPreset(ad_voice_freq, 40), rPreset(ad_voice_filter, 20),
              "Intensity of LFO"),
    rParamZyn(Pstartphase, rShort("start"), rSpecial(random),
              rDefaultDepends(loc), rDefault(64), rPreset(ad_voice_freq, 0),
              "Starting Phase"),
    rParamZyn(Pcutoff, rShort("lp"), rDefault(127),
            "RND/SQR lp-filter freq"),
    rOption(PLFOtype, rShort("type"), rOptions(sine, triangle, square, up, down,
                exp1, exp2, random), rLinear(0,127), rDefault(sine), "Shape of LFO"),
    rParamZyn(Prandomness, rShort("a.r."), rSpecial(disable), rDefault(0),
            "Amplitude Randomness (calculated uniformly at each cycle)"),
    rParamZyn(Pfreqrand, rShort("f.r."), rSpecial(disable), rDefault(0),
            "Frequency Randomness (calculated uniformly at each cycle)"),
    rParamF(delay, rShort("delay"), rSpecial(disable), rUnit(S),
              rLinear(0.0, 4.0), rDefaultDepends(loc), rDefault(0),
              rPreset(ad_voice_amp, 0.94),
              "Delay before LFO start\n0..4 second delay"),
    {"Pdelay::i", rShort("delay") rLinear(0,127)
     rDoc("Delay before LFO start\n0..4 second delay"), NULL,
    [](const char *msg, RtData &d)
     {
         rObject *obj = (rObject *)d.obj;
         if (!rtosc_narguments(msg)) {
             d.reply(d.loc, "i", (int)roundf(127.0f * obj->delay / 4.0f));
         } else {
             obj->delay = 4.0f * rtosc_argument(msg, 0).i / 127.0f;
         }
     }},
    rToggle(Pcontinous, rShort("c"), rDefault(false),
            "Enable for global operation"),
    rToggle(Psync, rShort("s"), rDefault(false),
            "Enable for BPM sync"),
    rParamZyn(Pstretch, rShort("str"), rCentered, rDefault(64),
        "Note frequency stretch"),
// these are currently not yet implemented and must be hidden therefore
#ifdef DEAD_PORTS
    //Float valued aliases
    {"delay::f", rProp(parameter) rMap(units, ms) rLog(0,4000), 0,
        rBegin;

        rEnd},
#define rPseudoLog(a,b) rLog(a,b)
    {"period::f", rProp(parameter) rMap(units, ms) rPseudoLog(0.10, 1500.0), 0,
        rBegin;
        rEnd},
#endif
};
#undef rPseudoLog
#undef rBegin
#undef rEnd
#undef rChangeCb

const rtosc::Ports &LFOParams::ports = _ports;

void LFOParams::setup()
{
    switch(loc) {
        case loc_unspecified:
            fel = consumer_location_type_t::unspecified;
            break;
        case ad_global_freq:
        case ad_voice_freq:
            fel = consumer_location_type_t::freq;
            setpresettype("Plfofrequency");
            break;
        case ad_global_amp:
        case ad_voice_amp:
            fel = consumer_location_type_t::amp;
            setpresettype("Plfoamplitude");
            break;
        case ad_global_filter:
        case ad_voice_filter:
            fel = consumer_location_type_t::filter;
            setpresettype("Plfofilter");
            break;
        default:
            throw std::logic_error("Invalid envelope consumer location");
    }

    defaults();
}

// TODO: reuse
LFOParams::LFOParams(const AbsTime *time_) :
    LFOParams(2.65, 0, 0, 0, 0, 0, 0, 0, 0, loc_unspecified, time_)
{
}

LFOParams::LFOParams(float freq_,
                     char Pintensity_,
                     char Pstartphase_,
                     char Pcutoff_,
                     char PLFOtype_,
                     char Prandomness_,
                     float Pdelay_,
                     char Pcontinous_,
                     char Psync_,
                     consumer_location_t loc,
                     const AbsTime *time_) : loc(loc),
                                             time(time_),
                                             last_update_timestamp(0) {
    Dfreq       = freq_;
    Dintensity  = Pintensity_;
    Dstartphase = Pstartphase_;
    Dcutoff     = Pcutoff_;
    DLFOtype    = PLFOtype_;
    Drandomness = Prandomness_;
    Ddelay      = Pdelay_;
    Dcontinous  = Pcontinous_;
    Dsync       = Psync_;

    setup();
}

LFOParams::LFOParams(consumer_location_t loc,
                     const AbsTime *time_) : loc(loc),
                                             time(time_),
                                             last_update_timestamp(0) {

    auto init =
        [&](float freq_, char Pintensity_, char Pstartphase_, char Pcutoff_, char PLFOtype_,
            char Prandomness_, float delay_, char Pcontinous_, char Psync_)
    {
        Dfreq       = freq_;
        Dintensity  = Pintensity_;
        Dstartphase = Pstartphase_;
        Dcutoff     = Pcutoff_;
        DLFOtype    = PLFOtype_;
        Drandomness = Prandomness_;
        Ddelay      = delay_;
        Dcontinous  = Pcontinous_;
        Dsync       = Psync_;
    };

    switch(loc)
    {
        case ad_global_amp:    init(6.49, 0, 64, 0, 0, 0, 0, 0, 0); break;
        case ad_global_freq:   init(3.71, 0, 64, 0, 0, 0, 0, 0, 0); break;
        case ad_global_filter: init(6.49, 0, 64, 0, 0, 0, 0, 0, 0); break;
        case ad_voice_amp:     init(11.25, 32, 64, 0, 0, 0.94, 0, 0, 0); break;
        case ad_voice_freq:    init(1.19, 40,  0, 0, 0,  0, 0, 0, 0); break;
        case ad_voice_filter:  init(1.19, 20, 64, 0, 0,  0, 0, 0, 0); break;
        default: throw std::logic_error("Invalid LFO consumer location");
    }

    setup();
}

LFOParams::~LFOParams()
{}

void LFOParams::defaults()
{
    freq        = Dfreq;
    Pintensity  = Dintensity;
    Pstartphase = Dstartphase;
    Pcutoff     = Dcutoff;
    PLFOtype    = DLFOtype;
    Prandomness = Drandomness;
    delay       = Ddelay;
    Pcontinous  = Dcontinous;
    Psync       = Dsync;
    Pfreqrand   = 0;
    Pstretch    = 64;
}


void LFOParams::add2XML(XMLwrapper& xml)
{
    xml.addparreal("freq", freq);
    xml.addpar("intensity", Pintensity);
    xml.addpar("start_phase", Pstartphase);
    xml.addpar("cutoff", Pcutoff);
    xml.addpar("lfo_type", PLFOtype);
    xml.addpar("randomness_amplitude", Prandomness);
    xml.addpar("randomness_frequency", Pfreqrand);
    xml.addparreal("delay", delay);
    xml.addpar("stretch", Pstretch);
    xml.addparbool("continous", Pcontinous);
    xml.addparbool("sync", Psync);
}

void LFOParams::getfromXML(XMLwrapper& xml)
{
    if (xml.fileversion() < version_type(3, 0, 4)) {
        freq = (powf(2.0f, 10.0f * xml.getparreal("freq", freq, 0.0f, 1.0f)) -1) / 12.0;
    } else {
        freq       = xml.getparreal("freq", freq);
    }
    Pintensity  = xml.getpar127("intensity", Pintensity);
    Pstartphase = xml.getpar127("start_phase", Pstartphase);
    Pcutoff = xml.getpar127("cutoff", Pcutoff);
    PLFOtype    = xml.getpar127("lfo_type", PLFOtype);
    Prandomness = xml.getpar127("randomness_amplitude", Prandomness);
    Pfreqrand   = xml.getpar127("randomness_frequency", Pfreqrand);
    if (xml.hasparreal("delay")) {
        delay      = xml.getparreal("delay", delay);
    } else {
        delay      = 4.0f * xml.getpar127("delay", (int)delay *127.0f/4.0f)
                     / 127.0f;
    }
    Pstretch    = xml.getpar127("stretch", Pstretch);
    Pcontinous  = xml.getparbool("continous", Pcontinous);
    Psync  = xml.getparbool("sync", Psync);
}

#define COPY(y) this->y=x.y
void LFOParams::paste(LFOParams &x)
{
    COPY(freq);
    COPY(Pintensity);
    COPY(Pstartphase);
    COPY(Pcutoff);
    COPY(PLFOtype);
    COPY(Prandomness);
    COPY(Pfreqrand);
    COPY(delay);
    COPY(Pcontinous);
    COPY(Psync);
    COPY(Pstretch);

    if ( time ) {
        last_update_timestamp = time->time();
    }
}
#undef COPY

}
