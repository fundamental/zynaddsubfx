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
        LFO(const LFOParams &lfopars_, float basefreq_, const AbsTime &t, WatchManager *m=0,
                const char *watch_prefix=0);
        ~LFO();

        float lfoout();
        float amplfoout();
        void releasekey();
    private:
        typedef enum lfo_state_type{
            delaying, 
            fadingIn,
            running,
            fadingOut
        } lfo_state_type;
    
        float baseOut(const char waveShape, const float phase);
        float biquad(float input);
        void updatePars();
        
        lfo_state_type lfo_state;
        
        //tempo stored to detect changes
        unsigned int tempo;
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
        // Ref to AbsTime object for time.tempo
        const AbsTime &time;
        //Delay before starting
        RelTime delayTime;

        int64_t fadeInDuration;
        //Timestamp of begin fadein
        int64_t fadeInTimestamp;
        //Timestamp of noteoff
        int64_t releaseTimestamp;
        //Time to ramp out
        
        int64_t fadeOutDuration;
        float rampUp, rampDown, rampOnRelease;
        // store the constant out value before oscillating starts
        float outStartValue = 0.0;

        char  waveShape;

        //If After initialization there are no calls to random number gen.
        bool  deterministic;

        const float     dt;
        const LFOParams &lfopars;
        const float basefreq;

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
