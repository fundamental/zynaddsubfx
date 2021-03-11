/*
  ZynAddSubFX - a software synthesizer

  LFO.h - LFO implementation
  Copyright (C) 2002-2005 Nasca Octavian Paul
  Author: Nasca Octavian Paul

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.
*/

#ifndef LFO_H
#define LFO_H

#include "../globals.h"
#include "../Misc/Time.h"
#include "WatchPoint.h"

namespace zyn {

/**Class for creating Low Frequency Oscillators*/
class LFO
{
    public:
        /**Constructor
         *
         * @param lfopars pointer to a LFOParams object
         * @param basefreq base frequency of LFO
         */
        LFO(const LFOParams &lfopars, float basefreq, const AbsTime &t, WatchManager *m=0,
                const char *watch_prefix=0);
        ~LFO();

        float lfoout();
        float amplfoout();
    private:
        float baseOut(const char waveShape, const float phase);
        float biquad(float input);
        //Phase of Oscillator
        float phase;
        //Phase Increment Per Frame
        float phaseInc;
        //Frequency Randomness
        float incrnd, nextincrnd;
        //Amplitude Randomness
        float amp1, amp2;

        // RND mode
        int first_half;
        float last_random;
        float z1, z2;

        //Intensity of the wave
        float lfointensity;
        //Amount Randomness
        float lfornd, lfofreqrnd;
        // Ref to AbsTime object for time.bpm
        const AbsTime& time;
        //Delay before starting
        RelTime delayTime;

        char  waveShape;

        //If After initialization there are no calls to random number gen.
        bool  deterministic;

        const float     dt_;
        const LFOParams &lfopars_;
        const float basefreq_;

        float FcAbs, K, norm;

        //biquad coefficients for lp filtering in noise-LFO
        float a0 = 0.0007508914611009499;
        float a1 = 0.0015017829222018998;
        float a2 = 0.0007508914611009499;
        float b1 = -1.519121359805288;
        float b2 =  0.5221249256496917;

        char cutoff = 127;

        VecWatchPoint watchOut;

        void computeNextFreqRnd(void);
};

}

#endif
