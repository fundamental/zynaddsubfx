/*
  ZynAddSubFX - a software synthesizer

  SEQ.cpp - SEQ implementation
  Copyright (C) 2002-2005 Nasca Octavian Paul
  Author: Nasca Octavian Paul

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.
*/

#include "SEQ.h"
#include "../Params/SEQParams.h"
#include "../Misc/Util.h"

#include <cstdlib>
#include <cstdio>
#include <cmath>

namespace zyn {

SEQ::SEQ(const SEQParams &seqpars, float basefreq, const AbsTime &t, WatchManager *m,
        const char *watch_prefix)
    :time(t),
    delayTime(t, seqpars.delay), //0..4 sec
    tStart(t.time()),
    dt_(t.dt()),
    seqpars_(seqpars), 
    basefreq_(basefreq),
    watchOut(m, watch_prefix, "out")
{

    //max 2x/octave
    if(!seqpars.speedratio) {
        if (seqpars.freq < 0.001) {
            duration = 1000.0f;
        }
        else
            duration = 1.0/seqpars.freq; 
    } else {
        duration = 60.0f / limit(((float(time.tempo)) * seqpars.speedratio), 0.001f, 1000.0f);
    }
    
    if (seqpars.continous)
        tRef = 0;
    else
        tRef = tStart;

    z1 = 0.0;
    z2 = 0.0;
    
    if(seqpars_.fel == consumer_location_type_t::amp)
        neutralOut = 1.0f;
    else
        neutralOut = 0.0f;

}

SEQ::~SEQ()
{}

float SEQ::baseOut()
{
    float fraction = (float)(time.time() - tRef) * time.dt() / duration;
    int index = (int)fraction;
    phase = fraction - (float)index;
    //~ printf("time.time(): %d, tRef: %d\n", time.time(), tRef);
    currentSegment = index % limit((int)seqpars_.steps, 1, 128);
    //~ printf("index: %d, currentSegment: %d\n", index, currentSegment);
    float seq_out = seqpars_.sequence[currentSegment];
    //~ printf("seq_out: %f\n", seq_out);
    return biquad(seq_out);
}


float SEQ::biquad(float input)
{
    float output;
    if (seqpars_.cutoff!=cutoff ) // calculate coeffs only if cutoff changed
    {
        cutoff = seqpars_.cutoff;
        if (cutoff!= MAX_CUTOFF) // at cutoff 0.0f we bypass filter, no coeffs needed
        {
            // calculate biquad coefficients
            FcAbs = limit(cutoff, 1.0f, 40.0f);
            K = tan(PI * limit(FcAbs * dt_,0.001f,0.4f)); // FcRel * dt_ max 40 * 0.01 = 0.4,
            // LIMIT in case of SEQ sampling frequency lower than 100 Hz

            norm = 1.0f / (1.0f + K / 0.7071f + K * K);
            a0 = K * K * norm;
            a1 = 2.0f * a0;
            a2 = a0;
            b1 = 2.0f * (K * K - 1.0f) * norm;
            b2 = (1.0f - K / 0.7071f + K * K) * norm;
        }
    }
    if (cutoff != MAX_CUTOFF) // at cutoff 127 we bypass filter, nothing to do
    {
        output = limit(input * a0 + z1, -1.0f, 1.0f);
        z1 = input * a1 + z2 - b1 * output;
        z2 = input * a2 - b2 * output;
    }
    return (cutoff==MAX_CUTOFF) ? input : output; // at cutoff 127 bypass filter
}


float SEQ::seqout()
{
    if(delayTime.inFuture() || seqpars_.steps == 0) 
    {
        return neutralOut;
    }
    //update internals
    if ( ! seqpars_.time || seqpars_.last_update_timestamp == seqpars_.time->time())
    {
        if((seqpars_.speedratio)<0.001) { 
            if (seqpars_.freq < 0.001) {
                return neutralOut;
            }  
            duration = 1.0/seqpars_.freq;
        } else {
            duration = 60.0f / limit(((float(time.tempo)) * seqpars_.speedratio), 0.001f, 10000.0f);
        }
    }

    float out = baseOut();
    
    switch(seqpars_.fel) {
        case consumer_location_type_t::amp:
            seqintensity = seqpars_.intensity / 100.0f;
            break;
        case consumer_location_type_t::filter:
            seqintensity = seqpars_.intensity / 25.0f;
            break; 
        case consumer_location_type_t::freq:
        case consumer_location_type_t::unspecified:
            seqintensity = seqpars_.intensity * 19.0f;
            break;
        }
    

    out *= seqintensity;
    
    float watch_data[2] = {((float)currentSegment+phase)/(float)NUM_SEQ_STEPS, out};
    watchOut(watch_data, 2);

    return out;
}

/*
 * SEQ out (for amplitude)
 */
float SEQ::ampseqout()
{
    return limit(1.0f - seqintensity + seqout(), -1.0f, 1.0f);
}

}
