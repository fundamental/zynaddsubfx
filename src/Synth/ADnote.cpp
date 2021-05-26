/*
  ZynAddSubFX - a software synthesizer

  ADnote.cpp - The "additive" synthesizer
  Copyright (C) 2002-2005 Nasca Octavian Paul
  Author: Nasca Octavian Paul

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.
*/

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cassert>
#include <stdint.h>

#include "../globals.h"
#include "../Misc/Util.h"
#include "../Misc/Allocator.h"
#include "../Params/ADnoteParameters.h"
#include "../Containers/ScratchString.h"
#include "../Containers/NotePool.h"
#include "ModFilter.h"
#include "OscilGen.h"
#include "ADnote.h"

#define LENGTHOF(x) ((int)(sizeof(x)/sizeof(x[0])))

namespace zyn {
ADnote::ADnote(ADnoteParameters *pars_, const SynthParams &spars,
        WatchManager *wm, const char *prefix)
    :SynthNote(spars), watch_be4_add(wm, prefix, "noteout/be4_mix"), watch_after_add(wm,prefix,"noteout/after_mix"),
    watch_punch(wm, prefix, "noteout/punch"), watch_legato(wm, prefix, "noteout/legato"), pars(*pars_)
{
    memory.beginTransaction();
    tmpwavel = memory.valloc<float>(synth.buffersize);
    tmpwaver = memory.valloc<float>(synth.buffersize);
    bypassl  = memory.valloc<float>(synth.buffersize);
    bypassr  = memory.valloc<float>(synth.buffersize);

    ADnoteParameters &pars = *pars_;
    portamento  = spars.portamento;
    note_log2_freq = spars.note_log2_freq;
    NoteEnabled = ON;
    velocity    = spars.velocity;
    initial_seed = spars.seed;
    current_prng_state = spars.seed;
    stereo = pars.GlobalPar.PStereo;

    NoteGlobalPar.Detune = getdetune(pars.GlobalPar.PDetuneType,
                                     pars.GlobalPar.PCoarseDetune,
                                     pars.GlobalPar.PDetune);
    bandwidthDetuneMultiplier = pars.getBandwidthDetuneMultiplier();

    if(pars.GlobalPar.PPanning == 0)
        NoteGlobalPar.Panning = getRandomFloat();
    else
        NoteGlobalPar.Panning = pars.GlobalPar.PPanning / 128.0f;

    NoteGlobalPar.Fadein_adjustment =
        pars.GlobalPar.Fadein_adjustment / (float)FADEIN_ADJUSTMENT_SCALE;
    NoteGlobalPar.Fadein_adjustment *= NoteGlobalPar.Fadein_adjustment;
    if(pars.GlobalPar.PPunchStrength != 0) {
        NoteGlobalPar.Punch.Enabled = 1;
        NoteGlobalPar.Punch.t = 1.0f; //start from 1.0f and to 0.0f
        NoteGlobalPar.Punch.initialvalue =
            ((powf(10, 1.5f * pars.GlobalPar.PPunchStrength / 127.0f) - 1.0f)
             * VelF(velocity,
                    pars.GlobalPar.PPunchVelocitySensing));
        float time =
            powf(10, 3.0f * pars.GlobalPar.PPunchTime / 127.0f) / 10000.0f;   //0.1f .. 100 ms
        float stretch = powf(440.0f / powf(2.0f, spars.note_log2_freq),
                             pars.GlobalPar.PPunchStretch / 64.0f);
        NoteGlobalPar.Punch.dt = 1.0f / (time * synth.samplerate_f * stretch);
    }
    else
        NoteGlobalPar.Punch.Enabled = 0;

    for(int nvoice = 0; nvoice < NUM_VOICES; ++nvoice)
        setupVoice(nvoice);

    max_unison = 1;
    for(int nvoice = 0; nvoice < NUM_VOICES; ++nvoice)
        if(NoteVoicePar[nvoice].unison_size > max_unison)
            max_unison = NoteVoicePar[nvoice].unison_size;


    tmpwave_unison = memory.valloc<float*>(max_unison);
    for(int k = 0; k < max_unison; ++k) {
        tmpwave_unison[k] = memory.valloc<float>(synth.buffersize);
        memset(tmpwave_unison[k], 0, synth.bufferbytes);
    }

    initparameters(wm, prefix);
    memory.endTransaction();
}

void ADnote::setupVoice(int nvoice)
{
    auto &param = pars.VoicePar[nvoice];
    auto &voice = NoteVoicePar[nvoice];


    for (int i = 0; i < 14; i++)
        voice.pinking[i] = 0.0;

    param.OscilGn->newrandseed(prng());
    voice.OscilSmp = NULL;
    voice.FMSmp    = NULL;
    voice.VoiceOut = NULL;

    voice.FMVoice = -1;
    voice.unison_size = 1;

    if(!pars.VoicePar[nvoice].Enabled) {
        voice.Enabled = OFF;
        return;   //the voice is disabled
    }
    NoteVoicePar[nvoice].AAEnabled =
                pars.VoicePar[nvoice].PAAEnabled;

    const int BendAdj = pars.VoicePar[nvoice].PBendAdjust - 64;
    if (BendAdj % 24 == 0)
        voice.BendAdjust = BendAdj / 24;
    else
        voice.BendAdjust = BendAdj / 24.0f;

    const float offset_val = (param.POffsetHz - 64)/64.0f;
    voice.OffsetHz   = 15.0f*(offset_val * sqrtf(fabsf(offset_val)));

    voice.unison_stereo_spread =
        pars.VoicePar[nvoice].Unison_stereo_spread / 127.0f;

    int unison = setupVoiceUnison(nvoice);


    voice.oscfreqhi   = memory.valloc<int>(unison);
    voice.oscfreqlo   = memory.valloc<float>(unison);

    voice.oscfreqhiFM = memory.valloc<unsigned int>(unison);
    voice.oscfreqloFM = memory.valloc<float>(unison);
    voice.oscposhi    = memory.valloc<int>(unison);
    voice.oscposlo    = memory.valloc<float>(unison);
    voice.oscposhiFM  = memory.valloc<unsigned int>(unison);
    voice.oscposloFM  = memory.valloc<float>(unison);

    voice.Enabled     = ON;
    voice.fixedfreq   = pars.VoicePar[nvoice].Pfixedfreq;
    voice.fixedfreqET = pars.VoicePar[nvoice].PfixedfreqET;

    setupVoiceDetune(nvoice);

    for(int k = 0; k < unison; ++k) {
        voice.oscposhi[k]   = 0;
        voice.oscposlo[k]   = 0.0f;
        voice.oscposhiFM[k] = 0;
        voice.oscposloFM[k] = 0.0f;
    }

    //the extra points contains the first point
    voice.OscilSmp =
        memory.valloc<float>(synth.oscilsize + OSCIL_SMP_EXTRA_SAMPLES);

    //Get the voice's oscil or external's voice oscil
    int vc = nvoice;
    if(pars.VoicePar[nvoice].Pextoscil != -1)
        vc = pars.VoicePar[nvoice].Pextoscil;
    if(!pars.GlobalPar.Hrandgrouping)
        pars.VoicePar[vc].OscilGn->newrandseed(prng());
    int oscposhi_start =
        pars.VoicePar[vc].OscilGn->get(NoteVoicePar[nvoice].OscilSmp,
                getvoicebasefreq(nvoice),
                pars.VoicePar[nvoice].Presonance);

    // This code was planned for biasing the carrier in MOD_RING
    // but that's on hold for the moment.  Disabled 'cos small
    // machines run this stuff too.
    //
    // //Find range of generated wave
    // float min = NoteVoicePar[nvoice].OscilSmp[0];
    // float max = min;
    // float *smpls = &(NoteVoicePar[nvoice].OscilSmp[1]);
    // for (int i = synth.oscilsize-1; i--; smpls++)
    //     if (*smpls > max)
    //         max = *smpls;
    //     else if (*smpls < min)
    //         min = *smpls;
    // NoteVoicePar[nvoice].OscilSmpMin = min;
    // NoteVoicePar[nvoice].OscilSmpMax = max;

    //I store the first elements to the last position for speedups
    for(int i = 0; i < OSCIL_SMP_EXTRA_SAMPLES; ++i)
        voice.OscilSmp[synth.oscilsize + i] = voice.OscilSmp[i];

    voice.phase_offset = (int)((pars.VoicePar[nvoice].Poscilphase
                    - 64.0f) / 128.0f * synth.oscilsize + synth.oscilsize * 4);
    oscposhi_start += NoteVoicePar[nvoice].phase_offset;

    int kth_start = oscposhi_start;
    for(int k = 0; k < unison; ++k) {
        voice.oscposhi[k] = kth_start % synth.oscilsize;
        //put random starting point for other subvoices
        kth_start      = oscposhi_start +
            (int)(RND * pars.VoicePar[nvoice].Unison_phase_randomness /
                    127.0f * (synth.oscilsize - 1));
    }

    voice.FreqLfo      = NULL;
    voice.FreqSeq      = NULL;
    voice.FreqEnvelope = NULL;

    voice.AmpLfo      = NULL;
    voice.AmpSeq      = NULL;
    voice.AmpEnvelope = NULL;

    voice.Filter         = NULL;
    voice.FilterEnvelope = NULL;
    voice.FilterLfo      = NULL;
    voice.FilterSeq      = NULL;

    voice.filterbypass = param.Pfilterbypass;

    setupVoiceMod(nvoice);

    voice.FMVoice = param.PFMVoice;
    voice.FMFreqEnvelope = NULL;
    voice.FMAmpEnvelope  = NULL;

    voice.FMoldsmp = memory.valloc<float>(unison);
    for(int k = 0; k < unison; ++k)
        voice.FMoldsmp[k] = 0.0f; //this is for FM (integration)

    voice.firsttick = 1;
    voice.DelayTicks =
        (int)((expf(param.PDelay / 127.0f * logf(50.0f))
                    - 1.0f) / synth.buffersize_f / 10.0f * synth.samplerate_f);
}

int ADnote::setupVoiceUnison(int nvoice)
{
    auto &voice = NoteVoicePar[nvoice];

    int unison = pars.VoicePar[nvoice].Unison_size;
    if(unison < 1)
        unison = 1;

    bool is_pwm = pars.VoicePar[nvoice].PFMEnabled == FMTYPE::PW_MOD;

    if (pars.VoicePar[nvoice].Type != 0) {
        // Since noise unison of greater than two is touch goofy...
        if (unison > 2)
            unison = 2;
    } else if (is_pwm) {
        /* Pulse width mod uses pairs of subvoices. */
        unison *= 2;
        // This many is likely to sound like noise anyhow.
        if (unison > 64)
            unison = 64;
    }

    //compute unison
    voice.unison_size = unison;

    voice.unison_base_freq_rap = memory.valloc<float>(unison);
    voice.unison_freq_rap      = memory.valloc<float>(unison);
    voice.unison_invert_phase  = memory.valloc<bool>(unison);
    const float unison_spread =
        pars.getUnisonFrequencySpreadCents(nvoice);
    const float unison_real_spread = powf(2.0f, (unison_spread * 0.5f) / 1200.0f);
    const float unison_vibratto_a  =
        pars.VoicePar[nvoice].Unison_vibratto / 127.0f; //0.0f .. 1.0f

    const int true_unison = unison / (is_pwm ? 2 : 1);
    switch(true_unison) {
        case 1:
            voice.unison_base_freq_rap[0] = 1.0f; //if the unison is not used, always make the only subvoice to have the default note
            break;
        case 2: //unison for 2 subvoices
            voice.unison_base_freq_rap[0] = 1.0f / unison_real_spread;
            voice.unison_base_freq_rap[1] = unison_real_spread;
            break;
        default: //unison for more than 2 subvoices
        {
            float unison_values[true_unison];
            float min = -1e-6f, max = 1e-6f;
            for(int k = 0; k < true_unison; ++k) {
                float step = (k / (float) (true_unison - 1)) * 2.0f - 1.0f; //this makes the unison spread more uniform
                float val  = step + (RND * 2.0f - 1.0f) / (true_unison - 1);
                unison_values[k] = val;
                if (min > val) {
                    min = val;
                }
                if (max < val) {
                    max = val;
                }
            }
            const float diff = max - min;
            for(int k = 0; k < true_unison; ++k) {
                unison_values[k] =
                    (unison_values[k] - (max + min) * 0.5f) / diff;         //the lowest value will be -1 and the highest will be 1
                voice.unison_base_freq_rap[k] =
                    powf(2.0f, (unison_spread * unison_values[k]) / 1200);
            }
        }
    }
    if (is_pwm)
        for (int i = true_unison - 1; i >= 0; i--) {
            voice.unison_base_freq_rap[2*i + 1] =
                voice.unison_base_freq_rap[i];
            voice.unison_base_freq_rap[2*i] =
                voice.unison_base_freq_rap[i];
        }

    //unison vibrattos
    if(unison > 2 || (!is_pwm && unison > 1))
        for(int k = 0; k < unison; ++k) //reduce the frequency difference for larger vibrattos
            voice.unison_base_freq_rap[k] = 1.0f
                + (voice.unison_base_freq_rap[k] - 1.0f)
                * (1.0f - unison_vibratto_a);
    voice.unison_vibratto.step      = memory.valloc<float>(unison);
    voice.unison_vibratto.position  = memory.valloc<float>(unison);
    voice.unison_vibratto.amplitude =
        (unison_real_spread - 1.0f) * unison_vibratto_a;

    const float increments_per_second = synth.samplerate_f / synth.buffersize_f;
    const float vib_speed = pars.VoicePar[nvoice].Unison_vibratto_speed / 127.0f;
    const float vibratto_base_period  = 0.25f * powf(2.0f, (1.0f - vib_speed) * 4.0f);
    for(int k = 0; k < unison; ++k) {
        voice.unison_vibratto.position[k] = RND * 1.8f - 0.9f;
        //make period to vary randomly from 50% to 200% vibratto base period
        const float vibratto_period = vibratto_base_period
            * powf(2.0f, RND * 2.0f - 1.0f);

        const float m = (RND < 0.5f ? -1.0f : 1.0f) *
            4.0f / (vibratto_period * increments_per_second);
        voice.unison_vibratto.step[k] = m;

        // Ugly, but the alternative is likely uglier.
        if (is_pwm)
            for (int i = 0; i < unison; i += 2) {
                voice.unison_vibratto.step[i+1] =
                    voice.unison_vibratto.step[i];
                voice.unison_vibratto.position[i+1] =
                    voice.unison_vibratto.position[i];
            }
    }

    if(unison <= 2) { //no vibratto for a single voice
        if (is_pwm) {
            voice.unison_vibratto.step[1]     = 0.0f;
            voice.unison_vibratto.position[1] = 0.0f;
        }
        if (is_pwm || unison == 1) {
            voice.unison_vibratto.step[0]     = 0.0f;
            voice.unison_vibratto.position[0] = 0.0f;
            voice.unison_vibratto.amplitude   = 0.0f;
        }
    }

    //phase invert for unison
    voice.unison_invert_phase[0] = false;
    if(unison != 1) {
        int inv = pars.VoicePar[nvoice].Unison_invert_phase;
        switch(inv) {
            case 0:
                for(int k = 0; k < unison; ++k)
                    voice.unison_invert_phase[k] = false;
                break;
            case 1:
                for(int k = 0; k < unison; ++k)
                    voice.unison_invert_phase[k] = (RND > 0.5f);
                break;
            default:
                for(int k = 0; k < unison; ++k)
                    voice.unison_invert_phase[k] =
                        (k % inv == 0) ? true : false;
                break;
        }
    }
    return unison;
}

void ADnote::setupVoiceDetune(int nvoice)
{
    //use the Globalpars.detunetype if the detunetype is 0
    if(pars.VoicePar[nvoice].PDetuneType != 0) {
        NoteVoicePar[nvoice].Detune = getdetune(
                pars.VoicePar[nvoice].PDetuneType,
                pars.VoicePar[nvoice].
                PCoarseDetune,
                8192); //coarse detune
        NoteVoicePar[nvoice].FineDetune = getdetune(
                pars.VoicePar[nvoice].PDetuneType,
                0,
                pars.VoicePar[nvoice].PDetune); //fine detune
    }
    else {
        NoteVoicePar[nvoice].Detune = getdetune(
                pars.GlobalPar.PDetuneType,
                pars.VoicePar[nvoice].
                PCoarseDetune,
                8192); //coarse detune
        NoteVoicePar[nvoice].FineDetune = getdetune(
                pars.GlobalPar.PDetuneType,
                0,
                pars.VoicePar[nvoice].PDetune); //fine detune
    }
    if(pars.VoicePar[nvoice].PFMDetuneType != 0)
        NoteVoicePar[nvoice].FMDetune = getdetune(
                pars.VoicePar[nvoice].PFMDetuneType,
                pars.VoicePar[nvoice].
                PFMCoarseDetune,
                pars.VoicePar[nvoice].PFMDetune);
    else
        NoteVoicePar[nvoice].FMDetune = getdetune(
                pars.GlobalPar.PDetuneType,
                pars.VoicePar[nvoice].
                PFMCoarseDetune,
                pars.VoicePar[nvoice].PFMDetune);
}

void ADnote::setupVoiceMod(int nvoice, bool first_run)
{
    auto &param = pars.VoicePar[nvoice];
    auto &voice = NoteVoicePar[nvoice];
    float FMVolume;

    if (param.Type != 0)
        voice.FMEnabled = FMTYPE::NONE;
    else
        voice.FMEnabled = param.PFMEnabled;

    voice.FMFreqFixed  = param.PFMFixedFreq;

    //Triggers when a user enables modulation on a running voice
    if(!first_run && voice.FMEnabled != FMTYPE::NONE && voice.FMSmp == NULL && voice.FMVoice < 0) {
        param.FmGn->newrandseed(prng());
        voice.FMSmp = memory.valloc<float>(synth.oscilsize + OSCIL_SMP_EXTRA_SAMPLES);
        memset(voice.FMSmp, 0, sizeof(float)*(synth.oscilsize + OSCIL_SMP_EXTRA_SAMPLES));
        int vc = nvoice;
        if(param.PextFMoscil != -1)
            vc = param.PextFMoscil;

        float tmp = 1.0f;
        if((pars.VoicePar[vc].FmGn->Padaptiveharmonics != 0)
                || (voice.FMEnabled == FMTYPE::MIX)
                || (voice.FMEnabled == FMTYPE::RING_MOD))
            tmp = getFMvoicebasefreq(nvoice);

        if(!pars.GlobalPar.Hrandgrouping)
            pars.VoicePar[vc].FmGn->newrandseed(prng());

        for(int k = 0; k < voice.unison_size; ++k)
            voice.oscposhiFM[k] = (voice.oscposhi[k]
                    + pars.VoicePar[vc].FmGn->get(
                        voice.FMSmp, tmp))
                % synth.oscilsize;

        for(int i = 0; i < OSCIL_SMP_EXTRA_SAMPLES; ++i)
            voice.FMSmp[synth.oscilsize + i] = voice.FMSmp[i];
        int oscposhiFM_add =
            (int)((param.PFMoscilphase
                        - 64.0f) / 128.0f * synth.oscilsize
                    + synth.oscilsize * 4);
        for(int k = 0; k < voice.unison_size; ++k) {
            voice.oscposhiFM[k] += oscposhiFM_add;
            voice.oscposhiFM[k] %= synth.oscilsize;
        }
    }


    //Compute the Voice's modulator volume (incl. damping)
    float fmvoldamp = powf(440.0f / getvoicebasefreq(nvoice),
            param.PFMVolumeDamp / 64.0f - 1.0f);
    const float fmvolume_ = param.FMvolume / 100.0f;
    switch(voice.FMEnabled) {
        case FMTYPE::PHASE_MOD:
        case FMTYPE::PW_MOD:
            fmvoldamp = powf(440.0f / getvoicebasefreq(nvoice),
                    param.PFMVolumeDamp / 64.0f);
            FMVolume = (expf(fmvolume_ * FM_AMP_MULTIPLIER) - 1.0f)
                * fmvoldamp * 4.0f;
            break;
        case FMTYPE::FREQ_MOD:
            FMVolume = (expf(fmvolume_ * FM_AMP_MULTIPLIER) - 1.0f)
                * fmvoldamp * 4.0f;
            break;
        default:
            if(fmvoldamp > 1.0f)
                fmvoldamp = 1.0f;
            FMVolume = fmvolume_ * fmvoldamp;
            break;
    }

    //Voice's modulator velocity sensing
    voice.FMVolume = FMVolume *
        VelF(velocity, pars.VoicePar[nvoice].PFMVelocityScaleFunction);
}

SynthNote *ADnote::cloneLegato(void)
{
    SynthParams sp{memory, ctl, synth, time, velocity,
                (bool)portamento, legato.param.note_log2_freq, true,
                initial_seed };
    return memory.alloc<ADnote>(&pars, sp);
}

// ADlegatonote: This function is (mostly) a copy of ADnote(...) and
// initparameters() stuck together with some lines removed so that it
// only alter the already playing note (to perform legato). It is
// possible I left stuff that is not required for this.
void ADnote::legatonote(const LegatoParams &lpars)
{
    //ADnoteParameters &pars = *partparams;
    // Manage legato stuff
    if(legato.update(lpars))
        return;

    portamento = lpars.portamento;
    note_log2_freq = lpars.note_log2_freq;
    initial_seed = lpars.seed;
    current_prng_state = lpars.seed;

    if(lpars.velocity > 1.0f)
        velocity = 1.0f;
    else
        velocity = lpars.velocity;

    const float basefreq = powf(2.0f, note_log2_freq);

    NoteGlobalPar.Detune = getdetune(pars.GlobalPar.PDetuneType,
                                     pars.GlobalPar.PCoarseDetune,
                                     pars.GlobalPar.PDetune);
    bandwidthDetuneMultiplier = pars.getBandwidthDetuneMultiplier();

    if(pars.GlobalPar.PPanning)
        NoteGlobalPar.Panning = pars.GlobalPar.PPanning / 128.0f;

    NoteGlobalPar.Filter->updateSense(velocity,
                pars.GlobalPar.PFilterVelocityScale,
                pars.GlobalPar.PFilterVelocityScaleFunction);


    for(int nvoice = 0; nvoice < NUM_VOICES; ++nvoice) {
        auto &voice = NoteVoicePar[nvoice];
        float FMVolume;

        if(voice.Enabled == OFF)
            continue;  //(gf) Stay the same as first note in legato.

        voice.fixedfreq   = pars.VoicePar[nvoice].Pfixedfreq;
        voice.fixedfreqET = pars.VoicePar[nvoice].PfixedfreqET;

        //use the Globalpars.detunetype if the detunetype is 0
        if(pars.VoicePar[nvoice].PDetuneType != 0) {
            voice.Detune = getdetune(
                pars.VoicePar[nvoice].PDetuneType,
                pars.VoicePar[nvoice].PCoarseDetune,
                8192); //coarse detune
            voice.FineDetune = getdetune(
                pars.VoicePar[nvoice].PDetuneType,
                0,
                pars.VoicePar[nvoice].PDetune); //fine detune
        }
        else {
            voice.Detune = getdetune(
                pars.GlobalPar.PDetuneType,
                pars.VoicePar[nvoice].PCoarseDetune,
                8192); //coarse detune
            voice.FineDetune = getdetune(
                pars.GlobalPar.PDetuneType,
                0,
                pars.VoicePar[nvoice].PDetune); //fine detune
        }
        if(pars.VoicePar[nvoice].PFMDetuneType != 0)
            voice.FMDetune = getdetune(
                pars.VoicePar[nvoice].PFMDetuneType,
                pars.VoicePar[nvoice].PFMCoarseDetune,
                pars.VoicePar[nvoice].PFMDetune);
        else
            voice.FMDetune = getdetune(
                pars.GlobalPar.PDetuneType,
                pars.VoicePar[nvoice].PFMCoarseDetune,
                pars.VoicePar[nvoice].PFMDetune);

        auto &voiceFilter = voice.Filter;
        if(voiceFilter) {
            const auto  &vce     = pars.VoicePar[nvoice];
            voiceFilter->updateSense(velocity, vce.PFilterVelocityScale,
                        vce.PFilterVelocityScaleFunction);
        }

        voice.filterbypass =
            pars.VoicePar[nvoice].Pfilterbypass;


        voice.FMVoice = pars.VoicePar[nvoice].PFMVoice;

        //Compute the Voice's modulator volume (incl. damping)
        float fmvoldamp = powf(440.0f / getvoicebasefreq(nvoice),
                               pars.VoicePar[nvoice].PFMVolumeDamp / 64.0f
                               - 1.0f);

        switch(voice.FMEnabled) {
            case FMTYPE::PHASE_MOD:
            case FMTYPE::PW_MOD:
                fmvoldamp =
                    powf(440.0f / getvoicebasefreq(
                             nvoice), pars.VoicePar[nvoice].PFMVolumeDamp
                         / 64.0f);
                FMVolume =
                    (expf(pars.VoicePar[nvoice].FMvolume / 100.0f
                          * FM_AMP_MULTIPLIER) - 1.0f) * fmvoldamp * 4.0f;
                break;
            case FMTYPE::FREQ_MOD:
                FMVolume =
                    (expf(pars.VoicePar[nvoice].FMvolume / 100.0f
                          * FM_AMP_MULTIPLIER) - 1.0f) * fmvoldamp * 4.0f;
                break;
            default:
                if(fmvoldamp > 1.0f)
                    fmvoldamp = 1.0f;
                FMVolume =
                    pars.VoicePar[nvoice].FMvolume
                    / 100.0f * fmvoldamp;
                break;
        }

        //Voice's modulator velocity sensing
        voice.FMVolume = FMVolume *
            VelF(velocity,
                 pars.VoicePar[nvoice].PFMVelocityScaleFunction);
    }
    ///    initparameters();

    ///////////////
    // Altered content of initparameters():

    int tmp[NUM_VOICES];

    NoteGlobalPar.Volume = dB2rap(pars.GlobalPar.Volume) //-60 dB .. 20 dB
                           * VelF(
                               velocity,
                               pars.GlobalPar.PAmpVelocityScaleFunction); //velocity sensing

    {
        auto        *filter  = NoteGlobalPar.Filter;
        filter->updateSense(velocity, pars.GlobalPar.PFilterVelocityScale,
                            pars.GlobalPar.PFilterVelocityScaleFunction);
        filter->updateNoteFreq(basefreq);
    }

    // Forbids the Modulation Voice to be greater or equal than voice
    for(int i = 0; i < NUM_VOICES; ++i)
        if(NoteVoicePar[i].FMVoice >= i)
            NoteVoicePar[i].FMVoice = -1;

    // Voice Parameter init
    for(unsigned nvoice = 0; nvoice < NUM_VOICES; ++nvoice) {
        Voice& vce = NoteVoicePar[nvoice];
        if(NoteVoicePar[nvoice].Enabled == 0)
            continue;

        NoteVoicePar[nvoice].noisetype = pars.VoicePar[nvoice].Type;
        /* Voice Amplitude Parameters Init */
        NoteVoicePar[nvoice].Volume =
            dB2rap(pars.VoicePar[nvoice].volume)             // -60 dB .. 0 dB
            * VelF(velocity,
                   pars.VoicePar[nvoice].PAmpVelocityScaleFunction); //velocity
        if(pars.VoicePar[nvoice].volume == -60.0)
            NoteVoicePar[nvoice].Volume = 0;

        if(pars.VoicePar[nvoice].PVolumeminus != 0)
            NoteVoicePar[nvoice].Volume = -NoteVoicePar[nvoice].Volume;

        NoteVoicePar[nvoice].AAEnabled =
                pars.VoicePar[nvoice].PAAEnabled;

        if(pars.VoicePar[nvoice].PPanning == 0) {
            NoteVoicePar[nvoice].Panning = getRandomFloat();
        } else
            NoteVoicePar[nvoice].Panning =
                pars.VoicePar[nvoice].PPanning / 128.0f;

        vce.newamplitude = 1.0f;
        if(pars.VoicePar[nvoice].PAmpEnvelopeEnabled
           && NoteVoicePar[nvoice].AmpEnvelope)
            vce.newamplitude *= NoteVoicePar[nvoice].AmpEnvelope->envout_dB();

        if(pars.VoicePar[nvoice].PAmpLfoEnabled && NoteVoicePar[nvoice].AmpLfo)
            vce.newamplitude *= NoteVoicePar[nvoice].AmpLfo->amplfoout();
            
        if(pars.VoicePar[nvoice].PAmpSeqEnabled && NoteVoicePar[nvoice].AmpSeq)
            vce.newamplitude *= NoteVoicePar[nvoice].AmpSeq->ampseqout();

        auto *voiceFilter = NoteVoicePar[nvoice].Filter;
        if(voiceFilter) {
            voiceFilter->updateSense(velocity, pars.VoicePar[nvoice].PFilterVelocityScale,
                                     pars.VoicePar[nvoice].PFilterVelocityScaleFunction);
            voiceFilter->updateNoteFreq(basefreq);
        }

        /* Voice Modulation Parameters Init */
        if((NoteVoicePar[nvoice].FMEnabled != FMTYPE::NONE)
           && (NoteVoicePar[nvoice].FMVoice < 0)) {
            pars.VoicePar[nvoice].FmGn->newrandseed(prng());

            //Perform Anti-aliasing only on MIX or RING MODULATION

            int vc = nvoice;
            if(pars.VoicePar[nvoice].PextFMoscil != -1)
                vc = pars.VoicePar[nvoice].PextFMoscil;

            if(!pars.GlobalPar.Hrandgrouping)
                pars.VoicePar[vc].FmGn->newrandseed(prng());

            for(int i = 0; i < OSCIL_SMP_EXTRA_SAMPLES; ++i)
                NoteVoicePar[nvoice].FMSmp[synth.oscilsize + i] =
                    NoteVoicePar[nvoice].FMSmp[i];
        }

        vce.FMnewamplitude = NoteVoicePar[nvoice].FMVolume
                                 * ctl.fmamp.relamp;

        if(pars.VoicePar[nvoice].PFMAmpEnvelopeEnabled
           && NoteVoicePar[nvoice].FMAmpEnvelope)
            vce.FMnewamplitude *=
                NoteVoicePar[nvoice].FMAmpEnvelope->envout_dB();
    }

    for(int nvoice = 0; nvoice < NUM_VOICES; ++nvoice) {
        for(unsigned i = nvoice + 1; i < NUM_VOICES; ++i)
            tmp[i] = 0;
        for(unsigned i = nvoice + 1; i < NUM_VOICES; ++i)
            if((NoteVoicePar[i].FMVoice == nvoice) && (tmp[i] == 0))
                tmp[i] = 1;
    }
}


/*
 * Kill a voice of ADnote
 */
void ADnote::KillVoice(int nvoice)
{
    auto &voice = NoteVoicePar[nvoice];

    memory.devalloc(voice.oscfreqhi);
    memory.devalloc(voice.oscfreqlo);
    memory.devalloc(voice.oscfreqhiFM);
    memory.devalloc(voice.oscfreqloFM);
    memory.devalloc(voice.oscposhi);
    memory.devalloc(voice.oscposlo);
    memory.devalloc(voice.oscposhiFM);
    memory.devalloc(voice.oscposloFM);

    memory.devalloc(voice.unison_base_freq_rap);
    memory.devalloc(voice.unison_freq_rap);
    memory.devalloc(voice.unison_invert_phase);
    memory.devalloc(voice.FMoldsmp);
    memory.devalloc(voice.unison_vibratto.step);
    memory.devalloc(voice.unison_vibratto.position);

    NoteVoicePar[nvoice].kill(memory, synth);
}

/*
 * Kill the note
 */
void ADnote::KillNote()
{
    for(unsigned nvoice = 0; nvoice < NUM_VOICES; ++nvoice) {
        if(NoteVoicePar[nvoice].Enabled == ON)
            KillVoice(nvoice);

        if(NoteVoicePar[nvoice].VoiceOut)
            memory.dealloc(NoteVoicePar[nvoice].VoiceOut);
    }

    NoteGlobalPar.kill(memory);

    NoteEnabled = OFF;
}

ADnote::~ADnote()
{
    if(NoteEnabled == ON)
        KillNote();
    memory.devalloc(tmpwavel);
    memory.devalloc(tmpwaver);
    memory.devalloc(bypassl);
    memory.devalloc(bypassr);
    for(int k = 0; k < max_unison; ++k)
        memory.devalloc(tmpwave_unison[k]);
    memory.devalloc(tmpwave_unison);
}


/*
 * Init the parameters
 */
void ADnote::initparameters(WatchManager *wm, const char *prefix)
{
    int tmp[NUM_VOICES];
    ScratchString pre = prefix;
    const float basefreq = powf(2.0f, note_log2_freq);

    // Global Parameters
    NoteGlobalPar.initparameters(pars.GlobalPar, synth,
                                 time,
                                 memory, basefreq, velocity,
                                 stereo, wm, prefix);

    NoteGlobalPar.AmpEnvelope->envout_dB(); //discard the first envelope output
    globalnewamplitude = NoteGlobalPar.Volume
                         * NoteGlobalPar.AmpEnvelope->envout_dB()
                         * NoteGlobalPar.AmpLfo->amplfoout()
                         * NoteGlobalPar.AmpSeq->ampseqout();

    // Forbids the Modulation Voice to be greater or equal than voice
    for(int i = 0; i < NUM_VOICES; ++i)
        if(NoteVoicePar[i].FMVoice >= i)
            NoteVoicePar[i].FMVoice = -1;

    // Voice Parameter init
    for(int nvoice = 0; nvoice < NUM_VOICES; ++nvoice) {
        Voice &vce = NoteVoicePar[nvoice];
        ADnoteVoiceParam &param = pars.VoicePar[nvoice];

        if(vce.Enabled == 0)
            continue;

        vce.noisetype = param.Type;
        /* Voice Amplitude Parameters Init */
        vce.Volume = dB2rap(param.volume) // -60dB..0dB
                     * VelF(velocity, param.PAmpVelocityScaleFunction);
        if(param.volume == -60.0f)
            vce.Volume = 0;

        if(param.PVolumeminus)
            vce.Volume = -vce.Volume;

        if(param.PPanning == 0) {
            vce.Panning = getRandomFloat();
        } else
            vce.Panning = param.PPanning / 128.0f;

        vce.newamplitude = 1.0f;
        if(param.PAmpEnvelopeEnabled) {
            vce.AmpEnvelope = memory.alloc<Envelope>(*param.AmpEnvelope,
                    basefreq, synth.dt(), wm,
                    (pre+"VoicePar"+nvoice+"/AmpEnvelope/").c_str);
            vce.AmpEnvelope->envout_dB(); //discard the first envelope sample
            vce.newamplitude *= vce.AmpEnvelope->envout_dB();
        }

        if(param.PAmpLfoEnabled) {
            vce.AmpLfo = memory.alloc<LFO>(*param.AmpLfo, basefreq, time, wm,
                    (pre+"VoicePar"+nvoice+"/AmpLfo/").c_str);
            vce.newamplitude *= vce.AmpLfo->amplfoout();
        }

        if(param.PAmpSeqEnabled) {
            vce.AmpSeq = memory.alloc<SEQ>(*param.AmpSeq, basefreq, time, wm,
                    (pre+"VoicePar"+nvoice+"/AmpSeq/").c_str);
            vce.newamplitude *= vce.AmpSeq->ampseqout();
        }

        /* Voice Frequency Parameters Init */
        if(param.PFreqEnvelopeEnabled)
            vce.FreqEnvelope = memory.alloc<Envelope>(*param.FreqEnvelope,
                    basefreq, synth.dt(), wm,
                    (pre+"VoicePar"+nvoice+"/FreqEnvelope/").c_str);

        if(param.PFreqLfoEnabled)
            vce.FreqLfo = memory.alloc<LFO>(*param.FreqLfo, basefreq, time, wm,
                    (pre+"VoicePar"+nvoice+"/FreqLfo/").c_str);

        if(param.PFreqSeqEnabled)
            vce.FreqSeq = memory.alloc<SEQ>(*param.FreqSeq, basefreq, time, wm,
                    (pre+"VoicePar"+nvoice+"/FreqSeq/").c_str);

        /* Voice Filter Parameters Init */
        if(param.PFilterEnabled) {
            vce.Filter = memory.alloc<ModFilter>(*param.VoiceFilter, synth, time, memory, stereo,
                      basefreq);
            vce.Filter->updateSense(velocity, param.PFilterVelocityScale,
                    param.PFilterVelocityScaleFunction);


            if(param.PFilterEnvelopeEnabled) {
                vce.FilterEnvelope =
                    memory.alloc<Envelope>(*param.FilterEnvelope,
                            basefreq, synth.dt(), wm,
                            (pre+"VoicePar"+nvoice+"/FilterEnvelope/").c_str);
                vce.Filter->addMod(*vce.FilterEnvelope);
            }

            if(param.PFilterLfoEnabled) {
                vce.FilterLfo = memory.alloc<LFO>(*param.FilterLfo, basefreq, time, wm,
                        (pre+"VoicePar"+nvoice+"/FilterLfo/").c_str);
                vce.Filter->addMod(*vce.FilterLfo);
            }

            if(param.PFilterSeqEnabled) {
                vce.FilterSeq = memory.alloc<SEQ>(*param.FilterSeq, basefreq, time, wm,
                        (pre+"VoicePar"+nvoice+"/FilterSeq/").c_str);
                vce.Filter->addMod(*vce.FilterSeq);
            }
        }

        /* Voice Modulation Parameters Init */
        if((vce.FMEnabled != FMTYPE::NONE) && (vce.FMVoice < 0)) {
            param.FmGn->newrandseed(prng());
            vce.FMSmp = memory.valloc<float>(synth.oscilsize + OSCIL_SMP_EXTRA_SAMPLES);

            //Perform Anti-aliasing only on MIX or RING MODULATION

            int vc = nvoice;
            if(param.PextFMoscil != -1)
                vc = param.PextFMoscil;

            float tmp = 1.0f;
            if((pars.VoicePar[vc].FmGn->Padaptiveharmonics != 0)
               || (vce.FMEnabled == FMTYPE::MIX)
               || (vce.FMEnabled == FMTYPE::RING_MOD))
                tmp = getFMvoicebasefreq(nvoice);

            if(!pars.GlobalPar.Hrandgrouping)
                pars.VoicePar[vc].FmGn->newrandseed(prng());

            for(int k = 0; k < vce.unison_size; ++k)
                vce.oscposhiFM[k] = (vce.oscposhi[k]
                                         + pars.VoicePar[vc].FmGn->get(
                                             vce.FMSmp, tmp))
                                        % synth.oscilsize;

            for(int i = 0; i < OSCIL_SMP_EXTRA_SAMPLES; ++i)
                vce.FMSmp[synth.oscilsize + i] = vce.FMSmp[i];
            int oscposhiFM_add =
                (int)((param.PFMoscilphase
                       - 64.0f) / 128.0f * synth.oscilsize
                      + synth.oscilsize * 4);
            for(int k = 0; k < vce.unison_size; ++k) {
                vce.oscposhiFM[k] += oscposhiFM_add;
                vce.oscposhiFM[k] %= synth.oscilsize;
            }
        }

        if(param.PFMFreqEnvelopeEnabled)
            vce.FMFreqEnvelope = memory.alloc<Envelope>(*param.FMFreqEnvelope,
                    basefreq, synth.dt(), wm,
                    (pre+"VoicePar"+nvoice+"/FMFreqEnvelope/").c_str);

        vce.FMnewamplitude = vce.FMVolume * ctl.fmamp.relamp;

        if(param.PFMAmpEnvelopeEnabled) {
            vce.FMAmpEnvelope =
                memory.alloc<Envelope>(*param.FMAmpEnvelope,
                        basefreq, synth.dt(), wm,
                        (pre+"VoicePar"+nvoice+"/FMAmpEnvelope/").c_str);
            vce.FMnewamplitude *= vce.FMAmpEnvelope->envout_dB();
        }
    }

    for(int nvoice = 0; nvoice < NUM_VOICES; ++nvoice) {
        for(int i = nvoice + 1; i < NUM_VOICES; ++i)
            tmp[i] = 0;
        for(int i = nvoice + 1; i < NUM_VOICES; ++i)
            if((NoteVoicePar[i].FMVoice == nvoice) && (tmp[i] == 0)) {
                NoteVoicePar[nvoice].VoiceOut =
                    memory.valloc<float>(synth.buffersize);
                tmp[i] = 1;
            }

        if(NoteVoicePar[nvoice].VoiceOut)
            memset(NoteVoicePar[nvoice].VoiceOut, 0, synth.bufferbytes);
    }
}


/*
 * Computes the relative frequency of each unison voice and it's vibratto
 * This must be called before setfreq* functions
 */
void ADnote::compute_unison_freq_rap(int nvoice) {
    Voice &vce = NoteVoicePar[nvoice];
    if(vce.unison_size == 1) { //no unison
        vce.unison_freq_rap[0] = 1.0f;
        return;
    }
    float relbw = ctl.bandwidth.relbw * bandwidthDetuneMultiplier;
    for(int k = 0; k < vce.unison_size; ++k) {
        float pos  = vce.unison_vibratto.position[k];
        float step = vce.unison_vibratto.step[k];
        pos += step;
        if(pos <= -1.0f) {
            pos  = -1.0f;
            step = -step;
        }
        if(pos >= 1.0f) {
            pos  = 1.0f;
            step = -step;
        }
        float vibratto_val = (pos - 0.333333333f * pos * pos * pos) * 1.5f; //make the vibratto lfo smoother
        vce.unison_freq_rap[k] = 1.0f
                                     + ((vce.unison_base_freq_rap[k]
                                         - 1.0f) + vibratto_val
                                        * vce.unison_vibratto.amplitude)
                                     * relbw;

        vce.unison_vibratto.position[k] = pos;
        step = vce.unison_vibratto.step[k] = step;
    }
}


/*
 * Computes the frequency of an oscillator
 */
void ADnote::setfreq(int nvoice, float in_freq)
{
    Voice &vce = NoteVoicePar[nvoice];
    for(int k = 0; k < vce.unison_size; ++k) {
        float freq  = fabsf(in_freq) * vce.unison_freq_rap[k];
        float speed = freq * synth.oscilsize_f / synth.samplerate_f;
        if(speed > synth.oscilsize_f)
            speed = synth.oscilsize_f;

        F2I(speed, vce.oscfreqhi[k]);
        vce.oscfreqlo[k] = speed - floorf(speed);

    }
}

/*
 * Computes the frequency of an modullator oscillator
 */
void ADnote::setfreqFM(int nvoice, float in_freq)
{
    Voice &vce = NoteVoicePar[nvoice];
    for(int k = 0; k < vce.unison_size; ++k) {
        float freq  = fabsf(in_freq) * vce.unison_freq_rap[k];
        float speed = freq * synth.oscilsize_f / synth.samplerate_f;
        if(speed > synth.samplerate_f)
            speed = synth.samplerate_f;

        F2I(speed, vce.oscfreqhiFM[k]);
        vce.oscfreqloFM[k] = speed - floorf(speed);
    }
}

/*
 * Get Voice base frequency
 */
float ADnote::getvoicebasefreq(int nvoice, float adjust_log2) const
{
    const float detune = NoteVoicePar[nvoice].Detune / 100.0f
                   + NoteVoicePar[nvoice].FineDetune / 100.0f
                   * ctl.bandwidth.relbw * bandwidthDetuneMultiplier
                   + NoteGlobalPar.Detune / 100.0f;

    if(NoteVoicePar[nvoice].fixedfreq == 0) {
        return powf(2.0f, note_log2_freq + detune / 12.0f + adjust_log2);
    }
    else { //the fixed freq is enabled
        const int fixedfreqET = NoteVoicePar[nvoice].fixedfreqET;
        float fixedfreq_log2 = log2f(440.0f);

        if(fixedfreqET != 0) { //if the frequency varies according the keyboard note
            float tmp_log2 = (note_log2_freq - fixedfreq_log2) *
                (powf(2.0f, (fixedfreqET - 1) / 63.0f) - 1.0f);
            if(fixedfreqET <= 64)
                fixedfreq_log2 += tmp_log2;
            else
                fixedfreq_log2 += tmp_log2 * log2f(3.0f);
        }
        return powf(2.0f, fixedfreq_log2 + detune / 12.0f + adjust_log2);
    }
}

/*
 * Get Voice's Modullator base frequency
 */
float ADnote::getFMvoicebasefreq(int nvoice) const
{
    return getvoicebasefreq(nvoice, NoteVoicePar[nvoice].FMDetune / 1200.0f);
}

/*
 * Computes all the parameters for each tick
 */
void ADnote::computecurrentparameters()
{
    const float relfreq = getFilterCutoffRelFreq();
    int   nvoice;
    float voicefreq, voicepitch, FMfreq,
          FMrelativepitch, globalpitch;

    globalpitch = 0.01f * (NoteGlobalPar.FreqEnvelope->envout()
                           + NoteGlobalPar.FreqLfo->lfoout() * ctl.modwheel.relmod
                           + NoteGlobalPar.FreqSeq->seqout() );

    globaloldamplitude = globalnewamplitude;
    globalnewamplitude = NoteGlobalPar.Volume
                         * NoteGlobalPar.AmpEnvelope->envout_dB()
                         * NoteGlobalPar.AmpLfo->amplfoout()
                         * NoteGlobalPar.AmpSeq->ampseqout();

    NoteGlobalPar.Filter->update(relfreq, ctl.filterq.relq);

    //compute the portamento, if it is used by this note
    float portamentofreqdelta_log2 = 0.0f;
    if(portamento != 0) { //this voice use portamento
        portamentofreqdelta_log2 = ctl.portamento.freqdelta_log2;
        if(ctl.portamento.used == 0) //the portamento has finished
            portamento = 0;  //this note is no longer "portamented"
    }

    //compute parameters for all voices
    for(nvoice = 0; nvoice < NUM_VOICES; ++nvoice) {
        Voice& vce = NoteVoicePar[nvoice];
        if(NoteVoicePar[nvoice].Enabled != ON)
            continue;
        NoteVoicePar[nvoice].DelayTicks -= 1;
        if(NoteVoicePar[nvoice].DelayTicks > 0)
            continue;

        compute_unison_freq_rap(nvoice);

        /*******************/
        /* Voice Amplitude */
        /*******************/
        vce.oldamplitude = vce.newamplitude;
        vce.newamplitude = 1.0f;

        if(NoteVoicePar[nvoice].AmpEnvelope)
            vce.newamplitude *= NoteVoicePar[nvoice].AmpEnvelope->envout_dB();

        if(NoteVoicePar[nvoice].AmpLfo)
            vce.newamplitude *= NoteVoicePar[nvoice].AmpLfo->amplfoout();

        if(NoteVoicePar[nvoice].AmpSeq)
            vce.newamplitude *= NoteVoicePar[nvoice].AmpSeq->ampseqout();

        /****************/
        /* Voice Filter */
        /****************/
        auto *voiceFilter = NoteVoicePar[nvoice].Filter;
        if(voiceFilter) {
            voiceFilter->update(relfreq, ctl.filterq.relq);
        }

        if(NoteVoicePar[nvoice].noisetype == 0) { //compute only if the voice isn't noise
            /*******************/
            /* Voice Frequency */
            /*******************/
            voicepitch = 0.0f;
            if(NoteVoicePar[nvoice].FreqLfo)
                voicepitch += NoteVoicePar[nvoice].FreqLfo->lfoout() / 100.0f
                              * ctl.bandwidth.relbw;

            if(NoteVoicePar[nvoice].FreqEnvelope)
                voicepitch += NoteVoicePar[nvoice].FreqEnvelope->envout()
                              / 100.0f;

            if(NoteVoicePar[nvoice].FreqSeq)
                voicepitch += NoteVoicePar[nvoice].FreqSeq->seqout()
                              / 100.0f;

            voicefreq = getvoicebasefreq(nvoice, portamentofreqdelta_log2 +
                (voicepitch + globalpitch) / 12.0f); //Hz frequency
            voicefreq *=
                powf(ctl.pitchwheel.relfreq, NoteVoicePar[nvoice].BendAdjust); //change the frequency by the controller
            setfreq(nvoice, voicefreq + NoteVoicePar[nvoice].OffsetHz);

            /***************/
            /*  Modulator */
            /***************/


            if(NoteVoicePar[nvoice].FMEnabled != FMTYPE::NONE) {
                FMrelativepitch = NoteVoicePar[nvoice].FMDetune / 100.0f;
                if(NoteVoicePar[nvoice].FMFreqEnvelope)
                    FMrelativepitch +=
                        NoteVoicePar[nvoice].FMFreqEnvelope->envout() / 100.0f;
                if (NoteVoicePar[nvoice].FMFreqFixed)
                    FMfreq = powf(2.0f, FMrelativepitch / 12.0f) * 440.0f;
                else
                    FMfreq = powf(2.0f, FMrelativepitch / 12.0f) * voicefreq;
                setfreqFM(nvoice, FMfreq);

                vce.FMoldamplitude = vce.FMnewamplitude;
                vce.FMnewamplitude = NoteVoicePar[nvoice].FMVolume
                                         * ctl.fmamp.relamp;
                if(NoteVoicePar[nvoice].FMAmpEnvelope)
                    vce.FMnewamplitude *=
                        NoteVoicePar[nvoice].FMAmpEnvelope->envout_dB();
            }
        }
    }
}


/*
 * Fadein in a way that removes clicks but keep sound "punchy"
 */
inline void ADnote::fadein(float *smps) const
{
    int zerocrossings = 0;
    for(int i = 1; i < synth.buffersize; ++i)
        if((smps[i - 1] < 0.0f) && (smps[i] > 0.0f))
            zerocrossings++;  //this is only the positive crossings

    float tmp = (synth.buffersize_f - 1.0f) / (zerocrossings + 1) / 3.0f;
    if(tmp < 8.0f)
        tmp = 8.0f;
    tmp *= NoteGlobalPar.Fadein_adjustment;

    int n;
    F2I(tmp, n); //how many samples is the fade-in
    if(n > synth.buffersize)
        n = synth.buffersize;
    for(int i = 0; i < n; ++i) { //fade-in
        float tmp = 0.5f - cosf((float)i / (float) n * PI) * 0.5f;
        smps[i] *= tmp;
    }
}

/*
 * Computes the Oscillator (Without Modulation) - LinearInterpolation
 */

/* As the code here is a bit odd due to optimization, here is what happens
 * First the current position and frequency are retrieved from the running
 * state. These are broken up into high and low portions to indicate how many
 * samples are skipped in one step and how many fractional samples are skipped.
 * Outside of this method the fractional samples are just handled with floating
 * point code, but that's a bit slower than it needs to be. In this code the low
 * portions are known to exist between 0.0 and 1.0 and it is known that they are
 * stored in single precision floating point IEEE numbers. This implies that
 * a maximum of 24 bits are significant. The below code does your standard
 * linear interpolation that you'll see throughout this codebase, but by
 * sticking to integers for tracking the overflow of the low portion, around 15%
 * of the execution time was shaved off in the ADnote test.
 */
inline void ADnote::ComputeVoiceOscillator_LinearInterpolation(int nvoice)
{
    Voice& vce = NoteVoicePar[nvoice];
    for(int k = 0; k < vce.unison_size; ++k) {
        int    poshi  = vce.oscposhi[k];
        // convert floating point fractional part (sample interval phase)
        // with range [0.0 ... 1.0] to fixed point with 1 digit is 2^-24
        // by multiplying with precalculated 2^24 and casting to integer:
        int    poslo  = (int)(vce.oscposlo[k] * 16777216.0f);
        int    freqhi = vce.oscfreqhi[k];
        // same for phase increment:
        int    freqlo = (int)(vce.oscfreqlo[k] * 16777216.0f);
        float *smps   = NoteVoicePar[nvoice].OscilSmp;
        float *tw     = tmpwave_unison[k];
        assert(vce.oscfreqlo[k] < 1.0f);
        for(int i = 0; i < synth.buffersize; ++i) {
            tw[i]  = (smps[poshi] * (0x01000000 - poslo) + smps[poshi + 1] * poslo)/(16777216.0f);
            poslo += freqlo;                // increment fractional part (sample interval phase)
            poshi += freqhi + (poslo>>24);  // add overflow over 24 bits in poslo to poshi
            poslo &= 0xffffff;              // remove overflow from poslo
            poshi &= synth.oscilsize - 1;   // remove overflow
        }
        vce.oscposhi[k] = poshi;
        vce.oscposlo[k] = poslo/(16777216.0f);
    }
}


/*
 * Computes the Oscillator (Without Modulation) - windowed sinc Interpolation
 */

/* As the code here is a bit odd due to optimization, here is what happens
 * First the current position and frequency are retrieved from the running
 * state. These are broken up into high and low portions to indicate how many
 * samples are skipped in one step and how many fractional samples are skipped.
 * Outside of this method the fractional samples are just handled with floating
 * point code, but that's a bit slower than it needs to be. In this code the low
 * portions are known to exist between 0.0 and 1.0 and it is known that they are
 * stored in single precision floating point IEEE numbers. This implies that
 * a maximum of 24 bits are significant. The below code does your standard
 * linear interpolation that you'll see throughout this codebase, but by
 * sticking to integers for tracking the overflow of the low portion, around 15%
 * of the execution time was shaved off in the ADnote test.
 */
inline void ADnote::ComputeVoiceOscillator_SincInterpolation(int nvoice)
{
    // windowed sinc kernel factor Fs*0.3, rejection 80dB
    const float_t kernel[] = {
        0.0010596256917418426f,
        0.004273442181254887f,
        0.0035466063043375785f,
        -0.014555483937137638f,
        -0.04789321342588484f,
        -0.050800020978553066f,
        0.04679847159974432f,
        0.2610646708018185f,
        0.4964802251145513f,
        0.6000513532962539f,
        0.4964802251145513f,
        0.2610646708018185f,
        0.04679847159974432f,
        -0.050800020978553066f,
        -0.04789321342588484f,
        -0.014555483937137638f,
        0.0035466063043375785f,
        0.004273442181254887f,
        0.0010596256917418426f
        };



    Voice& vce = NoteVoicePar[nvoice];
    for(int k = 0; k < vce.unison_size; ++k) {
        int    poshi  = vce.oscposhi[k];
        int    poslo  = (int)(vce.oscposlo[k] * (1<<24));
        int    freqhi = vce.oscfreqhi[k];
        int    freqlo = (int)(vce.oscfreqlo[k] * (1<<24));
        int    ovsmpfreqhi = vce.oscfreqhi[k] / 2;
        int    ovsmpfreqlo = (int)((vce.oscfreqlo[k] / 2) * (1<<24));

        int    ovsmpposlo;
        int    ovsmpposhi;
        int    uflow;
        float *smps   = NoteVoicePar[nvoice].OscilSmp;
        float *tw     = tmpwave_unison[k];
        assert(vce.oscfreqlo[k] < 1.0f);
        float out = 0;

        for(int i = 0; i < synth.buffersize; ++i) {
            ovsmpposlo  = poslo - (LENGTHOF(kernel)-1)/2 * ovsmpfreqlo;
            uflow = ovsmpposlo>>24;
            ovsmpposhi  = poshi - (LENGTHOF(kernel)-1)/2 * ovsmpfreqhi - ((0x00 - uflow) & 0xff);
            ovsmpposlo &= 0xffffff;
            ovsmpposhi &= synth.oscilsize - 1;
            out = 0;
            for (int l = 0; l<LENGTHOF(kernel); l++) {
                out += kernel[l] * (
                    smps[ovsmpposhi]     * ((1<<24) - ovsmpposlo) +
                    smps[ovsmpposhi + 1] * ovsmpposlo)/(1.0f*(1<<24));
                // advance to next kernel sample
                ovsmpposlo += ovsmpfreqlo;
                ovsmpposhi += ovsmpfreqhi + (ovsmpposlo>>24); // add the 24-bit overflow
                ovsmpposlo &= 0xffffff;
                ovsmpposhi &= synth.oscilsize - 1;

            }

            // advance to next sample
            poslo += freqlo;
            poshi += freqhi + (poslo>>24);
            poslo &= 0xffffff;
            poshi &= synth.oscilsize - 1;

            tw[i] = out;

        }
        vce.oscposhi[k] = poshi;
        vce.oscposlo[k] = poslo/(1.0f*(1<<24));
    }
}


/*
 * Computes the Oscillator (Mixing)
 */
inline void ADnote::ComputeVoiceOscillatorMix(int nvoice)
{
    ComputeVoiceOscillator_LinearInterpolation(nvoice);

    Voice& vce = NoteVoicePar[nvoice];
    if(vce.FMnewamplitude > 1.0f)
        vce.FMnewamplitude = 1.0f;
    if(vce.FMoldamplitude > 1.0f)
        vce.FMoldamplitude = 1.0f;

    if(NoteVoicePar[nvoice].FMVoice >= 0) {
        //if I use VoiceOut[] as modullator
        int FMVoice = NoteVoicePar[nvoice].FMVoice;
        for(int k = 0; k < vce.unison_size; ++k) {
            float *tw = tmpwave_unison[k];
            for(int i = 0; i < synth.buffersize; ++i) {
                float amp = INTERPOLATE_AMPLITUDE(vce.FMoldamplitude,
                                            vce.FMnewamplitude,
                                            i,
                                            synth.buffersize);
                tw[i] = tw[i]
                    * (1.0f - amp) + amp * NoteVoicePar[FMVoice].VoiceOut[i];
            }
        }
    }
    else
        for(int k = 0; k < vce.unison_size; ++k) {
            int    poshiFM  = vce.oscposhiFM[k];
            float  posloFM  = vce.oscposloFM[k];
            int    freqhiFM = vce.oscfreqhiFM[k];
            float  freqloFM = vce.oscfreqloFM[k];
            float *tw = tmpwave_unison[k];

            for(int i = 0; i < synth.buffersize; ++i) {
                float amp = INTERPOLATE_AMPLITUDE(vce.FMoldamplitude,
                                            vce.FMnewamplitude,
                                            i,
                                            synth.buffersize);
                tw[i] = tw[i] * (1.0f - amp) + amp
                        * (NoteVoicePar[nvoice].FMSmp[poshiFM] * (1 - posloFM)
                           + NoteVoicePar[nvoice].FMSmp[poshiFM + 1] * posloFM);
                posloFM += freqloFM;
                if(posloFM >= 1.0f) {
                    posloFM -= 1.0f;
                    poshiFM++;
                }
                poshiFM += freqhiFM;
                poshiFM &= synth.oscilsize - 1;
            }
            vce.oscposhiFM[k] = poshiFM;
            vce.oscposloFM[k] = posloFM;
        }
}

/*
 * Computes the Oscillator (Ring Modulation)
 */
inline void ADnote::ComputeVoiceOscillatorRingModulation(int nvoice)
{
    ComputeVoiceOscillator_LinearInterpolation(nvoice);

    Voice& vce = NoteVoicePar[nvoice];
    if(vce.FMnewamplitude > 1.0f)
        vce.FMnewamplitude = 1.0f;
    if(vce.FMoldamplitude > 1.0f)
        vce.FMoldamplitude = 1.0f;
    if(NoteVoicePar[nvoice].FMVoice >= 0)
        // if I use VoiceOut[] as modullator
        for(int k = 0; k < vce.unison_size; ++k) {
            float *tw = tmpwave_unison[k];
            for(int i = 0; i < synth.buffersize; ++i) {
                float amp = INTERPOLATE_AMPLITUDE(vce.FMoldamplitude,
                                            vce.FMnewamplitude,
                                            i,
                                            synth.buffersize);
                int FMVoice = NoteVoicePar[nvoice].FMVoice;
                tw[i] *= (1.0f - amp) + amp * NoteVoicePar[FMVoice].VoiceOut[i];
            }
        }
    else
        for(int k = 0; k < vce.unison_size; ++k) {
            int    poshiFM  = vce.oscposhiFM[k];
            float  posloFM  = vce.oscposloFM[k];
            int    freqhiFM = vce.oscfreqhiFM[k];
            float  freqloFM = vce.oscfreqloFM[k];
            float *tw = tmpwave_unison[k];

            for(int i = 0; i < synth.buffersize; ++i) {
                float amp = INTERPOLATE_AMPLITUDE(vce.FMoldamplitude,
                                            vce.FMnewamplitude,
                                            i,
                                            synth.buffersize);
                tw[i] *= (NoteVoicePar[nvoice].FMSmp[poshiFM] * (1.0f - posloFM)
                          + NoteVoicePar[nvoice].FMSmp[poshiFM
                                                       + 1] * posloFM) * amp
                         + (1.0f - amp);
                posloFM += freqloFM;
                if(posloFM >= 1.0f) {
                    posloFM -= 1.0f;
                    poshiFM++;
                }
                poshiFM += freqhiFM;
                poshiFM &= synth.oscilsize - 1;
            }
            vce.oscposhiFM[k] = poshiFM;
            vce.oscposloFM[k] = posloFM;
        }
}

/*
 * Computes the Oscillator (Phase Modulation or Frequency Modulation)
 */
inline void ADnote::ComputeVoiceOscillatorFrequencyModulation(int nvoice,
                                                              FMTYPE FMmode)
{
    Voice& vce = NoteVoicePar[nvoice];
    if(NoteVoicePar[nvoice].FMVoice >= 0) {
        //if I use VoiceOut[] as modulator
        for(int k = 0; k < vce.unison_size; ++k) {
            float *tw = tmpwave_unison[k];
            const float *smps = NoteVoicePar[NoteVoicePar[nvoice].FMVoice].VoiceOut;
            if (FMmode == FMTYPE::PW_MOD && (k & 1))
                for (int i = 0; i < synth.buffersize; ++i)
                    tw[i] = -smps[i];
            else
                memcpy(tw, smps, synth.bufferbytes);
        }
    } else {
        //Compute the modulator and store it in tmpwave_unison[][]
        for(int k = 0; k < vce.unison_size; ++k) {
            int    poshiFM  = vce.oscposhiFM[k];
            int    posloFM  = (int)(vce.oscposloFM[k]  * (1<<24));
            int    freqhiFM = vce.oscfreqhiFM[k];
            int    freqloFM = (int)(vce.oscfreqloFM[k] * (1<<24));
            float *tw = tmpwave_unison[k];
            const float *smps = NoteVoicePar[nvoice].FMSmp;

            for(int i = 0; i < synth.buffersize; ++i) {
                tw[i] = (smps[poshiFM] * ((1<<24) - posloFM)
                         + smps[poshiFM + 1] * posloFM) / (1.0f*(1<<24));
                if (FMmode == FMTYPE::PW_MOD && (k & 1))
                    tw[i] = -tw[i];

                posloFM += freqloFM;
                if(posloFM >= (1<<24)) {
                    posloFM &= 0xffffff;//fmod(posloFM, 1.0f);
                    poshiFM++;
                }
                poshiFM += freqhiFM;
                poshiFM &= synth.oscilsize - 1;
            }
            vce.oscposhiFM[k] = poshiFM;
            vce.oscposloFM[k] = posloFM/((1<<24)*1.0f);
        }
    }
    // Amplitude interpolation
    if(ABOVE_AMPLITUDE_THRESHOLD(vce.FMoldamplitude,
                                 vce.FMnewamplitude)) {
        for(int k = 0; k < vce.unison_size; ++k) {
            float *tw = tmpwave_unison[k];
            for(int i = 0; i < synth.buffersize; ++i)
                tw[i] *= INTERPOLATE_AMPLITUDE(vce.FMoldamplitude,
                                               vce.FMnewamplitude,
                                               i,
                                               synth.buffersize);
        }
    } else {
        for(int k = 0; k < vce.unison_size; ++k) {
            float *tw = tmpwave_unison[k];
            for(int i = 0; i < synth.buffersize; ++i)
                tw[i] *= vce.FMnewamplitude;
        }
    }


    //normalize: makes all sample-rates, oscil_sizes to produce same sound
    if(FMmode == FMTYPE::FREQ_MOD) { //Frequency modulation
        const float normalize = synth.oscilsize_f / 262144.0f * 44100.0f
                          / synth.samplerate_f;
        for(int k = 0; k < vce.unison_size; ++k) {
            float *tw    = tmpwave_unison[k];
            float  fmold = vce.FMoldsmp[k];
            for(int i = 0; i < synth.buffersize; ++i) {
                fmold = fmodf(fmold + tw[i] * normalize, synth.oscilsize);
                tw[i] = fmold;
            }
            vce.FMoldsmp[k] = fmold;
        }
    }
    else {  //Phase or PWM modulation
        const float normalize = synth.oscilsize_f / 262144.0f;
        for(int k = 0; k < vce.unison_size; ++k) {
            float *tw = tmpwave_unison[k];
            for(int i = 0; i < synth.buffersize; ++i)
                tw[i] *= normalize;
        }
    }

    //do the modulation
    for(int k = 0; k < vce.unison_size; ++k) {
        float *smps   = NoteVoicePar[nvoice].OscilSmp;
        float *tw     = tmpwave_unison[k];
        int    poshi  = vce.oscposhi[k];
        int    poslo  = (int)(vce.oscposlo[k] * (1<<24));
        int    freqhi = vce.oscfreqhi[k];
        int    freqlo = (int)(vce.oscfreqlo[k] * (1<<24));

        for(int i = 0; i < synth.buffersize; ++i) {
            int FMmodfreqhi = 0;
            F2I(tw[i], FMmodfreqhi);
            float FMmodfreqlo = tw[i]-FMmodfreqhi;//fmod(tw[i] /*+ 0.0000000001f*/, 1.0f);
            if(FMmodfreqhi < 0)
                FMmodfreqlo++;

            //carrier
            int carposhi = poshi + FMmodfreqhi;
            int carposlo = (int)(poslo + FMmodfreqlo);
            if (FMmode == FMTYPE::PW_MOD && (k & 1))
                carposhi += NoteVoicePar[nvoice].phase_offset;

            if(carposlo >= (1<<24)) {
                carposhi++;
                carposlo &= 0xffffff;//fmod(carposlo, 1.0f);
            }
            carposhi &= (synth.oscilsize - 1);

            tw[i] = (smps[carposhi] * ((1<<24) - carposlo)
                    + smps[carposhi + 1] * carposlo)/(1.0f*(1<<24));

            poslo += freqlo;
            if(poslo >= (1<<24)) {
                poslo &= 0xffffff;//fmod(poslo, 1.0f);
                poshi++;
            }

            poshi += freqhi;
            poshi &= synth.oscilsize - 1;
        }
        vce.oscposhi[k] = poshi;
        vce.oscposlo[k] = (poslo)/((1<<24)*1.0f);
    }
}


/*
 * Computes the Noise
 */
inline void ADnote::ComputeVoiceWhiteNoise(int nvoice)
{
    for(int k = 0; k < NoteVoicePar[nvoice].unison_size; ++k) {
        float *tw = tmpwave_unison[k];
        for(int i = 0; i < synth.buffersize; ++i)
            tw[i] = RND * 2.0f - 1.0f;
    }
}

inline void ADnote::ComputeVoicePinkNoise(int nvoice)
{
    Voice& vce = NoteVoicePar[nvoice];
    for(int k = 0; k < vce.unison_size; ++k) {
        float *tw = tmpwave_unison[k];
        float *f = &vce.pinking[k > 0 ? 7 : 0];
        for(int i = 0; i < synth.buffersize; ++i) {
            float white = (RND-0.5f)/4.0f;
            f[0] = 0.99886f*f[0]+white*0.0555179f;
            f[1] = 0.99332f*f[1]+white*0.0750759f;
            f[2] = 0.96900f*f[2]+white*0.1538520f;
            f[3] = 0.86650f*f[3]+white*0.3104856f;
            f[4] = 0.55000f*f[4]+white*0.5329522f;
            f[5] = -0.7616f*f[5]-white*0.0168980f;
            tw[i] = f[0]+f[1]+f[2]+f[3]+f[4]+f[5]+f[6]+white*0.5362f;
            f[6] = white*0.115926f;
        }
    }
}

inline void ADnote::ComputeVoiceDC(int nvoice)
{
    for(int k = 0; k < NoteVoicePar[nvoice].unison_size; ++k) {
        float *tw = tmpwave_unison[k];
        for(int i = 0; i < synth.buffersize; ++i)
            tw[i] = 1.0f;
    }
}



/*
 * Compute the ADnote samples
 * Returns 0 if the note is finished
 */
int ADnote::noteout(float *outl, float *outr)
{
    memcpy(outl, synth.denormalkillbuf, synth.bufferbytes);
    memcpy(outr, synth.denormalkillbuf, synth.bufferbytes);

    if(NoteEnabled == OFF)
        return 0;

    memset(bypassl, 0, synth.bufferbytes);
    memset(bypassr, 0, synth.bufferbytes);

    //Update Changed Parameters From UI
    for(unsigned nvoice = 0; nvoice < NUM_VOICES; ++nvoice) {
        if((NoteVoicePar[nvoice].Enabled != ON)
           || (NoteVoicePar[nvoice].DelayTicks > 0))
            continue;
        setupVoiceDetune(nvoice);
        setupVoiceMod(nvoice, false);
    }

    computecurrentparameters();

    for(unsigned nvoice = 0; nvoice < NUM_VOICES; ++nvoice) {
        if((NoteVoicePar[nvoice].Enabled != ON)
           || (NoteVoicePar[nvoice].DelayTicks > 0))
            continue;
        switch (NoteVoicePar[nvoice].noisetype) {
            case 0: //voice mode=sound
                switch(NoteVoicePar[nvoice].FMEnabled) {
                    case FMTYPE::MIX:
                        ComputeVoiceOscillatorMix(nvoice);
                        break;
                    case FMTYPE::RING_MOD:
                        ComputeVoiceOscillatorRingModulation(nvoice);
                        break;
                    case FMTYPE::FREQ_MOD:
                    case FMTYPE::PHASE_MOD:
                    case FMTYPE::PW_MOD:
                        ComputeVoiceOscillatorFrequencyModulation(nvoice,
                                                                  NoteVoicePar[nvoice].FMEnabled);
                        break;
                    default:
                        if(NoteVoicePar[nvoice].AAEnabled) ComputeVoiceOscillator_SincInterpolation(nvoice);
                        else ComputeVoiceOscillator_LinearInterpolation(nvoice);
                        //if (config.cfg.Interpolation) ComputeVoiceOscillator_CubicInterpolation(nvoice);
                }
                break;
            case 1:
                ComputeVoiceWhiteNoise(nvoice);
                break;
            case 2:
                ComputeVoicePinkNoise(nvoice);
                break;
            default:
                ComputeVoiceDC(nvoice);
                break;
        }
        // Voice Processing

        Voice& vce = NoteVoicePar[nvoice];
        //mix subvoices into voice
        memset(tmpwavel, 0, synth.bufferbytes);
        if(stereo)
            memset(tmpwaver, 0, synth.bufferbytes);
        for(int k = 0; k < vce.unison_size; ++k) {
            float *tw = tmpwave_unison[k];
            if(stereo) {
                float stereo_pos = 0;
                bool is_pwm = NoteVoicePar[nvoice].FMEnabled == FMTYPE::PW_MOD;
                if (is_pwm) {
                    if(vce.unison_size > 2)
                        stereo_pos = k/2
                            / (float)(vce.unison_size/2
                                      - 1) * 2.0f - 1.0f;
                } else if(vce.unison_size > 1) {
                    stereo_pos = k
                        / (float)(vce.unison_size
                                  - 1) * 2.0f - 1.0f;
                }
                float stereo_spread = vce.unison_stereo_spread * 2.0f; //between 0 and 2.0f
                if(stereo_spread > 1.0f) {
                    float stereo_pos_1 = (stereo_pos >= 0.0f) ? 1.0f : -1.0f;
                    stereo_pos =
                        (2.0f
                         - stereo_spread) * stereo_pos
                        + (stereo_spread - 1.0f) * stereo_pos_1;
                }
                else
                    stereo_pos *= stereo_spread;

                if(vce.unison_size == 1 ||
                   (is_pwm && vce.unison_size == 2))
                    stereo_pos = 0.0f;
                float panning = (stereo_pos + 1.0f) * 0.5f;


                float lvol = (1.0f - panning) * 2.0f;
                if(lvol > 1.0f)
                    lvol = 1.0f;

                float rvol = panning * 2.0f;
                if(rvol > 1.0f)
                    rvol = 1.0f;

                if(vce.unison_invert_phase[k]) {
                    lvol = -lvol;
                    rvol = -rvol;
                }

                for(int i = 0; i < synth.buffersize; ++i)
                    tmpwavel[i] += tw[i] * lvol;
                for(int i = 0; i < synth.buffersize; ++i)
                    tmpwaver[i] += tw[i] * rvol;
            }
            else
                for(int i = 0; i < synth.buffersize; ++i)
                    tmpwavel[i] += tw[i];
            if(nvoice == 0)
                watch_be4_add(tmpwavel,synth.buffersize);
        }

        float unison_amplitude = 1.0f / sqrtf(vce.unison_size); //reduce the amplitude for large unison sizes
        // Amplitude
        float oldam = vce.oldamplitude * unison_amplitude;
        float newam = vce.newamplitude * unison_amplitude;

        if(ABOVE_AMPLITUDE_THRESHOLD(oldam, newam)) {
            int rest = synth.buffersize;
            //test if the amplitude if raising and the difference is high
            if((newam > oldam) && ((newam - oldam) > 0.25f)) {
                rest = 10;
                if(rest > synth.buffersize)
                    rest = synth.buffersize;
                for(int i = 0; i < synth.buffersize - rest; ++i)
                    tmpwavel[i] *= oldam;
                if(stereo)
                    for(int i = 0; i < synth.buffersize - rest; ++i)
                        tmpwaver[i] *= oldam;
            }
            // Amplitude interpolation
            for(int i = 0; i < rest; ++i) {
                float amp = INTERPOLATE_AMPLITUDE(oldam, newam, i, rest);
                tmpwavel[i + (synth.buffersize - rest)] *= amp;
                if(stereo)
                    tmpwaver[i + (synth.buffersize - rest)] *= amp;
            }
        }
        else {
            for(int i = 0; i < synth.buffersize; ++i)
                tmpwavel[i] *= newam;
            if(stereo)
                for(int i = 0; i < synth.buffersize; ++i)
                    tmpwaver[i] *= newam;
        }

        // Fade in
        if(vce.firsttick != 0) {
            fadein(&tmpwavel[0]);
            if(stereo)
                fadein(&tmpwaver[0]);
            vce.firsttick = 0;
        }

        // Filter
        if(NoteVoicePar[nvoice].Filter) {
            if(stereo)
                NoteVoicePar[nvoice].Filter->filter(tmpwavel, tmpwaver);
            else
                NoteVoicePar[nvoice].Filter->filter(tmpwavel, 0);
        }

        //check if the amplitude envelope is finished, if yes, the voice will be fadeout
        if(NoteVoicePar[nvoice].AmpEnvelope)
            if(NoteVoicePar[nvoice].AmpEnvelope->finished()) {
                for(int i = 0; i < synth.buffersize; ++i)
                    tmpwavel[i] *= 1.0f - (float)i / synth.buffersize_f;
                if(stereo)
                    for(int i = 0; i < synth.buffersize; ++i)
                        tmpwaver[i] *= 1.0f - (float)i / synth.buffersize_f;
            }
        //the voice is killed later


        // Put the ADnote samples in VoiceOut (without applying Global volume, because I wish to use this voice as a modullator)
        if(NoteVoicePar[nvoice].VoiceOut) {
            if(stereo)
                for(int i = 0; i < synth.buffersize; ++i)
                    NoteVoicePar[nvoice].VoiceOut[i] = tmpwavel[i]
                                                       + tmpwaver[i];
            else   //mono
                for(int i = 0; i < synth.buffersize; ++i)
                    NoteVoicePar[nvoice].VoiceOut[i] = tmpwavel[i];
        }


        // Add the voice that do not bypass the filter to out
        if(NoteVoicePar[nvoice].filterbypass == 0) { //no bypass
            if(stereo)
                for(int i = 0; i < synth.buffersize; ++i) { //stereo
                    outl[i] += tmpwavel[i] * NoteVoicePar[nvoice].Volume
                               * NoteVoicePar[nvoice].Panning * 2.0f;
                    outr[i] += tmpwaver[i] * NoteVoicePar[nvoice].Volume
                               * (1.0f - NoteVoicePar[nvoice].Panning) * 2.0f;
                }
            else
                for(int i = 0; i < synth.buffersize; ++i) //mono
                    outl[i] += tmpwavel[i] * NoteVoicePar[nvoice].Volume;
        }
        else {  //bypass the filter
            if(stereo)
                for(int i = 0; i < synth.buffersize; ++i) { //stereo
                    bypassl[i] += tmpwavel[i] * NoteVoicePar[nvoice].Volume
                                  * NoteVoicePar[nvoice].Panning * 2.0f;
                    bypassr[i] += tmpwaver[i] * NoteVoicePar[nvoice].Volume
                                  * (1.0f
                                     - NoteVoicePar[nvoice].Panning) * 2.0f;
                }
            else
                for(int i = 0; i < synth.buffersize; ++i) //mono
                    bypassl[i] += tmpwavel[i] * NoteVoicePar[nvoice].Volume;
        }
        // check if there is necessary to process the voice longer (if the Amplitude envelope isn't finished)
        if(NoteVoicePar[nvoice].AmpEnvelope)
            if(NoteVoicePar[nvoice].AmpEnvelope->finished())
                KillVoice(nvoice);
    }

    //Processing Global parameters
    if(stereo) {
        NoteGlobalPar.Filter->filter(outl, outr);
    } else { //set the right channel=left channel
        NoteGlobalPar.Filter->filter(outl, 0);
        memcpy(outr, outl, synth.bufferbytes);
        memcpy(bypassr, bypassl, synth.bufferbytes);
    }

    for(int i = 0; i < synth.buffersize; ++i) {
        outl[i] += bypassl[i];
        outr[i] += bypassr[i];
    }

    if(ABOVE_AMPLITUDE_THRESHOLD(globaloldamplitude, globalnewamplitude))
        // Amplitude Interpolation
        for(int i = 0; i < synth.buffersize; ++i) {
            float tmpvol = INTERPOLATE_AMPLITUDE(globaloldamplitude,
                                                 globalnewamplitude,
                                                 i,
                                                 synth.buffersize);
            outl[i] *= tmpvol * NoteGlobalPar.Panning;
            outr[i] *= tmpvol * (1.0f - NoteGlobalPar.Panning);
        }
    else
        for(int i = 0; i < synth.buffersize; ++i) {
            outl[i] *= globalnewamplitude * NoteGlobalPar.Panning;
            outr[i] *= globalnewamplitude * (1.0f - NoteGlobalPar.Panning);
        }

    //Apply the punch
    if(NoteGlobalPar.Punch.Enabled != 0)
        for(int i = 0; i < synth.buffersize; ++i) {
            float punchamp = NoteGlobalPar.Punch.initialvalue
                             * NoteGlobalPar.Punch.t + 1.0f;
            outl[i] *= punchamp;
            outr[i] *= punchamp;
            NoteGlobalPar.Punch.t -= NoteGlobalPar.Punch.dt;
            if(NoteGlobalPar.Punch.t < 0.0f) {
                NoteGlobalPar.Punch.Enabled = 0;
                break;
            }
        }

    watch_punch(outl, synth.buffersize);
    watch_after_add(outl,synth.buffersize);

    // Apply legato-specific sound signal modifications
    legato.apply(*this, outl, outr);

    watch_legato(outl, synth.buffersize);

    // Check if the global amplitude is finished.
    // If it does, disable the note
    if(NoteGlobalPar.AmpEnvelope->finished()) {
        for(int i = 0; i < synth.buffersize; ++i) { //fade-out
            float tmp = 1.0f - (float)i / synth.buffersize_f;
            outl[i] *= tmp;
            outr[i] *= tmp;
        }
        KillNote();
    }
    return 1;
}


/*
 * Release the key (NoteOff)
 */
void ADnote::releasekey()
{
    for(int nvoice = 0; nvoice < NUM_VOICES; ++nvoice)
        NoteVoicePar[nvoice].releasekey();
    NoteGlobalPar.FreqEnvelope->releasekey();
    NoteGlobalPar.FilterEnvelope->releasekey();
    NoteGlobalPar.AmpEnvelope->releasekey();
}

/*
 * Check if the note is finished
 */
bool ADnote::finished() const
{
    if(NoteEnabled == ON)
        return 0;
    else
        return 1;
}

void ADnote::entomb(void)
{
    NoteGlobalPar.AmpEnvelope->forceFinish();
}

void ADnote::Voice::releasekey()
{
    if(!Enabled)
        return;
    if(AmpEnvelope)
        AmpEnvelope->releasekey();
    if(FreqEnvelope)
        FreqEnvelope->releasekey();
    if(FilterEnvelope)
        FilterEnvelope->releasekey();
    if(FMFreqEnvelope)
        FMFreqEnvelope->releasekey();
    if(FMAmpEnvelope)
        FMAmpEnvelope->releasekey();
}

void ADnote::Voice::kill(Allocator &memory, const SYNTH_T &synth)
{
    memory.devalloc(OscilSmp);
    memory.dealloc(FreqEnvelope);
    memory.dealloc(FreqLfo);
    memory.dealloc(FreqSeq);
    memory.dealloc(AmpEnvelope);
    memory.dealloc(AmpLfo);
    memory.dealloc(AmpSeq);
    memory.dealloc(Filter);
    memory.dealloc(FilterEnvelope);
    memory.dealloc(FilterLfo);
    memory.dealloc(FilterSeq);
    memory.dealloc(FMFreqEnvelope);
    memory.dealloc(FMAmpEnvelope);

    if((FMEnabled != FMTYPE::NONE) && (FMVoice < 0))
        memory.devalloc(FMSmp);

    if(VoiceOut)
        memset(VoiceOut, 0, synth.bufferbytes);
    //the buffer can't be safely deleted as it may be
    //an input to another voice

    Enabled = OFF;
}

void ADnote::Global::kill(Allocator &memory)
{
    memory.dealloc(FreqEnvelope);
    memory.dealloc(FreqLfo);
    memory.dealloc(FreqSeq);
    memory.dealloc(AmpEnvelope);
    memory.dealloc(AmpLfo);
    memory.dealloc(AmpSeq);
    memory.dealloc(Filter);
    memory.dealloc(FilterEnvelope);
    memory.dealloc(FilterLfo);
    memory.dealloc(FilterSeq);
}

void ADnote::Global::initparameters(const ADnoteGlobalParam &param,
                                    const SYNTH_T &synth,
                                    const AbsTime &time,
                                    class Allocator &memory,
                                    float basefreq, float velocity,
                                    bool stereo,
                                    WatchManager *wm,
                                    const char *prefix)
{
    ScratchString pre = prefix;
    FreqEnvelope = memory.alloc<Envelope>(*param.FreqEnvelope, basefreq,
            synth.dt(), wm, (pre+"GlobalPar/FreqEnvelope/").c_str);
    FreqLfo      = memory.alloc<LFO>(*param.FreqLfo, basefreq, time, wm,
                   (pre+"GlobalPar/FreqLfo/").c_str);
    FreqSeq      = memory.alloc<SEQ>(*param.FreqSeq, basefreq, time, wm,
                   (pre+"GlobalPar/FreqSeq/").c_str);

    AmpEnvelope = memory.alloc<Envelope>(*param.AmpEnvelope, basefreq,
            synth.dt(), wm, (pre+"GlobalPar/AmpEnvelope/").c_str);
    AmpLfo      = memory.alloc<LFO>(*param.AmpLfo, basefreq, time, wm,
                   (pre+"GlobalPar/AmpLfo/").c_str);
    AmpSeq      = memory.alloc<SEQ>(*param.AmpSeq, basefreq, time, wm,
                   (pre+"GlobalPar/AmpSeq/").c_str);

    Volume = dB2rap(param.Volume)
             * VelF(velocity, param.PAmpVelocityScaleFunction);     //sensing

    Filter = memory.alloc<ModFilter>(*param.GlobalFilter, synth, time, memory,
            stereo, basefreq);

    FilterEnvelope = memory.alloc<Envelope>(*param.FilterEnvelope, basefreq,
            synth.dt(), wm, (pre+"GlobalPar/FilterEnvelope/").c_str);
    FilterLfo      = memory.alloc<LFO>(*param.FilterLfo, basefreq, time, wm,
                   (pre+"GlobalPar/FilterLfo/").c_str);
    FilterSeq      = memory.alloc<SEQ>(*param.FilterSeq, basefreq, time, wm,
                   (pre+"GlobalPar/FilterSeq/").c_str);

    Filter->addMod(*FilterEnvelope);
    Filter->addMod(*FilterLfo);
    Filter->addMod(*FilterSeq);

    {
        Filter->updateSense(velocity, param.PFilterVelocityScale,
                param.PFilterVelocityScaleFunction);
    }
}

}
