/*
  ZynAddSubFX - a software synthesizer

  SEQ.h - SEQ implementation
  Copyright (C) 2002-2005 Nasca Octavian Paul
  Author: Nasca Octavian Paul

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.
*/

#ifndef SEQ_H
#define SEQ_H

#include "../globals.h"
#include "../Misc/Time.h"
#include "WatchPoint.h"

namespace zyn {

/**Class for creating Step Sequencers*/
class SEQ
{
    public:
        /**Constructor
         *
         * @param seqpars pointer to a SEQParams object
         * @param basefreq base frequency of SEQ
         */
        SEQ(const SEQParams &seqpars, float basefreq, const AbsTime &t, WatchManager *m=0,
                const char *watch_prefix=0);
        ~SEQ();

        float seqout();
        float ampseqout();
    private:
        float baseOut();
        float biquad(float input);

        //step duration
        float duration;
        
        // biquad filter state
        float z1, z2;

        // Ref to AbsTime object for time.bpm
        const AbsTime &time;
        //Delay before starting
        RelTime delayTime;
        //ReferenceTime
        int64_t tStart;
        int64_t tRef;
        // store the constant out value before oscillating starts
        float outConst;
        // step index
        int currentSegment;
        float phase;

        //If After initialization there are no calls to random number gen.
        bool  deterministic;

        const float     dt_;
        const SEQParams &seqpars_;
        const float basefreq_;

        float seqintensity;

        float FcAbs, K, norm;

        //biquad coefficients for lp filtering in noise-SEQ
        float a0 = 0.0007508914611009499;
        float a1 = 0.0015017829222018998;
        float a2 = 0.0007508914611009499;
        float b1 = -1.519121359805288;
        float b2 =  0.5221249256496917;

        float cutoff = 0.5f;
        
        float neutralOut = 0.0f;
        
        VecWatchPoint watchOut;
};

}

#endif
