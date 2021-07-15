/*
  ZynAddSubFX - a software synthesizer

  LFO.cpp - LFO implementation
  Copyright (C) 2002-2005 Nasca Octavian Paul
  Author: Nasca Octavian Paul

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.
*/

#include "LFO.h"
#include "../Params/LFOParams.h"
#include "../Misc/Util.h"

#include <cstdlib>
#include <cstdio>
#include <cmath>

namespace zyn {

LFO::LFO(const LFOParams &lfopars_, float basefreq_, const AbsTime &t, WatchManager *m,
        const char *watch_prefix)
    :first_half(-1),
    time(t),
    delayTime(t, lfopars_.delay), //0..4 sec
    deterministic(!lfopars_.Pfreqrand),
    dt(t.dt()),
    lfopars(lfopars_), 
    basefreq(basefreq_),
    watchOut(m, watch_prefix, "out")
{
    updatePars();
    
    if(!lfopars.Pcontinous) {
        if(!lfopars.Pstartphase)
            phase = RND;
        else
            phase = 0.0f;
    }
    else {
        phase = fmod((float)t.time() * phaseInc, 1.0f);
    }

    lfornd = limit(lfopars.Prandomness / 127.0f, 0.0f, 1.0f);
    lfofreqrnd = powf(lfopars.Pfreqrand / 127.0f, 2.0f) * 4.0f;

    switch(lfopars.fel) {
        case consumer_location_type_t::amp:
            lfointensity = lfopars.Pintensity / 127.0f;
            break;
        case consumer_location_type_t::filter:
            lfointensity = lfopars.Pintensity / 127.0f * 4.0f;
            break; //in octave
        case consumer_location_type_t::freq:
        case consumer_location_type_t::unspecified:
            lfointensity = powf(2, lfopars.Pintensity / 127.0f * 11.0f) - 1.0f; //in centi
            phase -= 0.25f; //chance the starting phase
            break;
    }

    lfo_state=lfo_state_type::delaying;
    
    rampUp = 0.0f;
    rampDown = 1.0f;

    amp1     = (1 - lfornd) + lfornd * RND;
    amp2     = (1 - lfornd) + lfornd * RND;
    incrnd   = nextincrnd = 1.0f;
    computeNextFreqRnd();
    computeNextFreqRnd(); //twice because I want incrnd & nextincrnd to be random
    z1 = 0.0;
    z2 = 0.0;
}

LFO::~LFO()
{}

void LFO::updatePars()
{
    waveShape = lfopars.PLFOtype;
    int stretch = lfopars.Pstretch;
    if(stretch == 0)
        stretch = 1;

    // stretch max 2x/octave
    const float lfostretch = powf(basefreq / 440.0f, (stretch - 64.0f) / 63.0f);

    float lfofreq;
    if (!lfopars.numerator || !lfopars.denominator) {
        lfofreq = lfopars.freq * lfostretch;   
    } else {
        tempo = time.tempo;
        lfofreq = float(tempo) * float(lfopars.denominator)/(240.0f * float(lfopars.numerator));
    }
    phaseInc = fabsf(lfofreq) * dt;
    
    //Limit the Frequency(or else...)
    if(phaseInc > 0.49999999f)
        phaseInc = 0.499999999f;
}

float LFO::baseOut(const char waveShape, const float phase)
{
    float lfo_out;
    switch(waveShape) {
        case LFO_TRIANGLE:
            if(phase >= 0.0f && phase < 0.25f)
                return 4.0f * phase;
            else if(phase > 0.25f && phase < 0.75f)
                return 2 - 4 * phase;
            else
                return 4.0f * phase - 4.0f;
            break;
        case LFO_SQUARE:
            if(phase < 0.5f)
                lfo_out = -1;
            else
                lfo_out = 1;

            return biquad(lfo_out);
            break;
        case LFO_RAMPUP:    return (phase - 0.5f) * 2.0f;
        case LFO_RAMPDOWN:  return (0.5f - phase) * 2.0f;
        case LFO_EXP_DOWN1: return powf(0.05f, phase) * 2.0f - 1.0f;
        case LFO_EXP_DOWN2: return powf(0.001f, phase) * 2.0f - 1.0f;
        case LFO_RANDOM:
            if ((phase < 0.5) != first_half) {
                first_half = phase < 0.5;
                last_random = 2*RND-1;
            }
            return biquad(last_random);
            break;
        default:
            return cosf(phase * 2.0f * PI); //LFO_SINE
    }
}

float LFO::biquad(float input)
{
    float output;
    if (lfopars.Pcutoff!=cutoff ) // calculate coeffs only if cutoff changed
    {
        cutoff = lfopars.Pcutoff;
        if (cutoff != 127) // at cutoff 127 we bypass filter, no coeffs needed
        {
            // calculate biquad coefficients
            FcAbs = (cutoff + 7.0f)*(cutoff + 7.0f)/ 450.56f; // max value < 40
            K = tan(PI * limit(FcAbs * dt,0.001f,0.4f)); // FcRel * dt_ max 40 * 0.01 = 0.4,
            // LIMIT in case of LFO sampling frequency lower than 100 Hz

            norm = 1 / (1 + K / 0.7071f + K * K);
            a0 = K * K * norm;
            a1 = 2 * a0;
            a2 = a0;
            b1 = 2 * (K * K - 1) * norm;
            b2 = (1 - K / 0.7071f + K * K) * norm;
        }
    }
    if (cutoff != 127) // at cutoff 127 we bypass filter, nothing to do
    {
        // lp filter the (s&h) random LFO
        output = limit(input * a0 + z1, -1.0f, 1.0f);
        z1 = input * a1 + z2 - b1 * output;
        z2 = input * a2 - b2 * output;
    }
    return (cutoff==127) ? input : output; // at cutoff 127 bypass filter
}

void LFO::releasekey()
{
    if (lfopars.fadeout==10.0f) { // deactivated
        fadeOutDuration = 0.0f;
        return;
    }
    // store current ramp value in case of release while fading in
    rampOnRelease = rampUp;
    // burn in the current amount of outStartValue    
    // therefor multiply its current reverse ramping factor 
    // and divide by current ramp factor. It will be multiplied during fading out.
    outStartValue *= (1.0f - rampOnRelease);
    // store current time
    releaseTimestamp = lfopars.time->time();
    // calculate fade out duration in frames
    fadeOutDuration = lfopars.fadeout * lfopars.time->framesPerSec();
    // set fadeout state
    lfo_state = lfo_state_type::fadingOut;
}
float LFO::lfoout()
{
    //update internals
    if ( ! lfopars.time || lfopars.last_update_timestamp == lfopars.time->time())
    {
        updatePars();
        switch(lfopars.fel) {
        
            case consumer_location_type_t::amp:
                lfointensity = lfopars.Pintensity / 127.0f; // [0...1]
                break;
            case consumer_location_type_t::filter:
                lfointensity = lfopars.Pintensity / 127.0f * 4.0f; // [0...4] octaves
                break; 
            case consumer_location_type_t::freq:
            case consumer_location_type_t::unspecified:
                lfointensity = powf(2, lfopars.Pintensity / 127.0f * 11.0f) - 1.0f; // [0...2047] cent
                break;
        }
    }
    
    // refresh freq if tempo has changed
    if (lfopars.numerator && lfopars.denominator && tempo != time.tempo) {
        tempo = time.tempo;
        float lfofreq = float(tempo) * float(lfopars.denominator)/(240.0f * float(lfopars.numerator));
        phaseInc = fabsf(lfofreq) * dt;
    }
    float phaseWithStartphase = fmod(phase + (lfopars.Pstartphase + 63.0f) / 127.0f, 1.0f);
    float out = baseOut(waveShape, phaseWithStartphase);
    if(waveShape == LFO_SINE || waveShape == LFO_TRIANGLE)
        out *= lfointensity * (amp1 + phaseWithStartphase * (amp2 - amp1));
    else
        out *= lfointensity * amp2;
    
    
    // handle lfo state (delay, fade in, fade out)
    switch(lfo_state) {
        case delaying:
        
            outStartValue = out; // keep start value to prevent jump
            if (delayTime.inFuture()) {
                return out;
            }else{
                fadeInTimestamp = lfopars.time->time();
                fadeInDuration = lfopars.fadein * lfopars.time->framesPerSec();
                lfo_state = lfo_state_type::fadingIn;
            }
            
            break;

        case fadingIn:
            
            if (fadeInDuration && rampUp < 1.0) {
                rampUp = ((float)(lfopars.time->time() - fadeInTimestamp) / (float)fadeInDuration);
                rampUp *= rampUp; // square for soft start
                    
            } 
            else {
                rampUp = 1.0f;
                lfo_state = lfo_state_type::running;
            } 
            
            out *= rampUp;
            out += outStartValue * (1.0f-rampUp);
            
            break;
        
        case fadingOut:
            if(fadeOutDuration && rampDown) {// no division by zero, please
                rampDown = 1.0f - ( (float)(lfopars.time->time() - releaseTimestamp) / (float)fadeOutDuration );
                rampDown *= rampDown; // square for soft end
            }
            else // no ramp down
                rampDown = 0.0f; 
            
            
            out *= rampOnRelease * rampDown;
            out += outStartValue*rampDown;
            
            break;
        
        case running:
        default:
            break;
    }

    //Start oscillating
    if(deterministic)
        phase += phaseInc;
    else {
        const float tmp = (incrnd * (1.0f - phase) + nextincrnd * phase);
        phase += phaseInc * limit(tmp, 0.0f, 1.0f);
    }
    if(phase >= 1) {
        phase    = fmod(phase, 1.0f);
        amp1 = amp2;
        amp2 = (1 - lfornd) + lfornd * RND;

        computeNextFreqRnd();
    }
            
    float watch_data[2] = {phaseWithStartphase, out};
    watchOut(watch_data, 2);

    return out;
}

/*
 * LFO out (for amplitude)
 */
float LFO::amplfoout()
{
    return limit(1.0f - lfointensity + lfoout(), -1.0f, 1.0f);
}


void LFO::computeNextFreqRnd()
{
    if(deterministic)
        return;
    incrnd     = nextincrnd;
    nextincrnd = powf(0.5f, lfofreqrnd) + RND * (powf(2.0f, lfofreqrnd) - 1.0f);
}

}
