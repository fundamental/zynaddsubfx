/*
  ZynAddSubFX - a software synthesizer

  Envelope.cpp - Envelope implementation
  Copyright (C) 2002-2005 Nasca Octavian Paul
  Author: Nasca Octavian Paul

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.
*/

#include <cmath>
#include "Envelope.h"
#include "../Params/EnvelopeParams.h"

namespace zyn {

Envelope::Envelope(EnvelopeParams &pars, float basefreq, float bufferdt,
        WatchManager *m, const char *watch_prefix)
    :watchOut(m, watch_prefix, "out")
{
    envpoints = pars.Penvpoints;
    if(envpoints > MAX_ENVELOPE_POINTS)
        envpoints = MAX_ENVELOPE_POINTS;
    envsustain     = (pars.Penvsustain == 0) ? -1 : pars.Penvsustain;
    forcedrelease   = pars.Pforcedrelease;
    envstretch     = powf(440.0f / basefreq, pars.Penvstretch / 64.0f);
    linearenvelope = pars.Plinearenvelope;
    repeating = pars.Prepeating;

    if(!pars.Pfreemode)
        pars.converttofree();

    mode = pars.Envmode;

    //for amplitude envelopes
    if((mode == 1) && !linearenvelope)
        mode = 2;                              //change to log envelope
    if((mode == 2) && linearenvelope)
        mode = 1;                              //change to linear

    for(int i = 0; i < MAX_ENVELOPE_POINTS; ++i) {
        const float dtstretched = pars.getdt(i) * envstretch;
        if(dtstretched > bufferdt)
            envdt[i] = bufferdt / dtstretched;
        else
            envdt[i] = 2.0f;  //any value larger than 1

        switch(mode) {
            case 2:
                envval[i] = (1.0f - pars.Penvval[i] / 127.0f) * -40;
                break;
            case 3:
                envval[i] =
                    (powf(2, 6.0f
                          * fabsf(pars.Penvval[i]
                                 - 64.0f) / 64.0f) - 1.0f) * 100.0f;
                if(pars.Penvval[i] < 64)
                    envval[i] = -envval[i];
                break;
            case 4:
                envval[i] = (pars.Penvval[i] - 64.0f) / 64.0f * 6.0f; //6 octaves (filtru)
                break;
            case 5:
                envval[i] = (pars.Penvval[i] - 64.0f) / 64.0f * 10;
                break;
            default:
                envval[i] = pars.Penvval[i] / 127.0f;
        }
    }

    envdt[0] = 1.0f;

    currentpoint = 1; //the envelope starts from 1
    keyreleased  = false;
    t = 0.0f;
    envfinish = false;
    inct      = envdt[1];
    envoutval = 0.0f;
}

Envelope::~Envelope()
{}


/*
 * Release the key (note envelope)
 */
void Envelope::releasekey()
{
    if(keyreleased)
        return;
    keyreleased = true;
    if(forcedrelease)
        t = 0.0f;
}

void Envelope::forceFinish(void)
{
    envfinish = true;
}

void Envelope::watch(float time, float value)
{
    float pos[2];
    float factor1;
    float factor2;
    pos[0] = time;
    switch(mode) {
        case 2:
            pos[1] = 1 - value / -40.f;
            watchOut(pos, 2);
            break;
        case 3:
            factor1 = log(value/100. + 1.) / (6. * log(2));
            factor2 = log(1. - value/100.) / (6. * log(2));
            pos[1] = ((0.5 * factor1) >= 0) ? (0.5 * factor1 + 0.5) : (0.5  - factor2 * 0.5);
            watchOut(pos, 2);
            break;
        case 4:
            pos[1] = (value + 6.) / 12.f;
            watchOut(pos, 2);
            break;
        case 5:
            pos[1] = (value + 10.) / 20.f;
            watchOut(pos, 2);
            break;
        default:
            pos[1] = value;
            watchOut(pos, 2);
    }
}

inline float lerp(float a, float b, float t) 
{
    return a + (b-a)*t;
}

inline float bezier(float a, float bRel, float c, float w2) 
{
    // square bezier direct form p = (1-t)^2 *P0 + 2*(1-t)*t*P1 + t*t*P2
    // formula in zest visualization:
    // y[i]+y[i-1])/2 + ((y[i]-y[i-1]).abs+0.1)*4*c[i]
    float w1 = 1.0f - w2;
    float b = (a+c) + (fabs(c-a)+0.1)*8.0f*bRel;
    return w1*w1*a + w1*w2*b + w2*w2*c;
}

/*
 * Envelope Output
 */
float Envelope::envout(bool doWatch)
{
    float out;
    if(envfinish) { //if the envelope is finished
        envoutval = envval[envpoints - 1];
        if(doWatch) {
            watch(envpoints - 1, envoutval);
        }
        return envoutval;
    }

    if((currentpoint == envsustain + 1) && !keyreleased) { //if it is sustaining now
        envoutval = envval[envsustain];
        bool zerorelease = true;
        for (auto i = envsustain; i<envpoints; i++)
            if (envval[i] != -40.0f) zerorelease = false;
        if (zerorelease &&                             //if sustaining at zero with zero until env ends
            (mode == ADSR_lin || mode == ADSR_dB)) {   // and its an amp envelope
            envfinish = true;   // finish voice to free ressources
        }
        if(doWatch) {
            watch(envsustain, envoutval);
        }
        return envoutval;
    }

    if(keyreleased && forcedrelease) { //do the forced release
        int releaseindex = (envsustain < 0) ? (envpoints - 1) : (envsustain + 1); //if there is no sustain point, use the last point for release

        if(envdt[releaseindex] < 0.00000001f)
            out = envval[releaseindex];
        else
            out = envoutval + (envval[releaseindex] - envoutval) * t; // linear interpolation envoutval and envval[releaseindex]

        t += envdt[releaseindex];

        if(t >= 1.0f) { // move to the next segment
            currentpoint = envsustain + 2;
            forcedrelease = 0;
            t    = 0.0f;
            inct = envdt[currentpoint];
            if((currentpoint >= envpoints) || (envsustain < 0))
                envfinish = true;
        }

        if(doWatch) {
            watch(releaseindex + t, envoutval);
        }

        return out;
    }
    if(inct >= 1.0f)
        out = envval[currentpoint];
    else // interpolation
        out = bezier(envval[currentpoint - 1], envcp[currentpoint-1], envval[currentpoint], t);

    t += inct;

    if(t >= 1.0f) {
        if(currentpoint >= envpoints - 1) // if last point reached
            envfinish = true;
        // but if reached sustain point, repeating activated and key still pressed or sustained
        else if (repeating && currentpoint == envsustain && !keyreleased) {
            // set first value to sustain value to prevent jump
            envval[0] = envval[currentpoint];
            // reset current point
            currentpoint = 1;
        }
        // otherwise proceed to the next segment
        else currentpoint++;

        t    = 0.0f;
        inct = envdt[currentpoint];
    }

    envoutval = out;

    if(doWatch) {
        watch(currentpoint + t, envoutval);
    }
    return out;
}

/*
 * Envelope Output (dB)
 */
float Envelope::envout_dB()
{
    float out;
    if(linearenvelope)
        return envout(true);

    if((currentpoint == 1) && (!keyreleased || !forcedrelease)) { //first point is always lineary interpolated <- seems to have odd effects
        float v1 = EnvelopeParams::env_dB2rap(envval[0]);
        float v2 = EnvelopeParams::env_dB2rap(envval[1]);
        out = v1 + (v2 - v1) * t;

        t += inct;

        if(t >= 1.0f) {
            t    = 0.0f;
            inct = envdt[2];
            currentpoint++;
            out = v2;
        }

        if(out > 0.001f)
            envoutval = EnvelopeParams::env_rap2dB(out);
        else
            envoutval = MIN_ENVELOPE_DB;
        out = envoutval;
    } else
        out = envout(false);

    watch(currentpoint + t, out);
    return EnvelopeParams::env_dB2rap(out);

}

bool Envelope::finished() const
{
    return envfinish;
}

}
