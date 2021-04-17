/*
  ZynAddSubFX - a software synthesizer

  AdNoteTest.h - CxxTest for Synth/ADnote
  Copyright (C) 2009-2011 Mark McCurry
  Copyright (C) 2009 Harald Hvaal
  Authors: Mark McCurry, Harald Hvaal

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.
*/


#include <cxxtest/TestSuite.h>
#include <iostream>
#include <fstream>
#include <ctime>
#include <string>
#include "../Misc/Util.h"
#include "../Misc/Allocator.h"
#include "../Synth/ADnote.h"
#include "../Synth/OscilGen.h"
#include "../Params/Presets.h"
#include "../Params/FilterParams.h"
#include "../DSP/FFTwrapper.h"
#include "../globals.h"
using namespace std;
using namespace zyn;

SYNTH_T *synth;

#define BUF 256
class AdNoteTest:public CxxTest::TestSuite
{
    public:

        ADnote       *note;
        FFTwrapper   *fft;
        Controller   *controller;
        float test_freq_log2;
        ADnoteParameters *params;
        AbsTime  *time;
        Alloc     memory;


        float outR[BUF], outL[BUF];

        void setUp() {
            //First the sensible settings and variables that have to be set:
            synth = new SYNTH_T;
            synth->buffersize = BUF;
            synth->alias();
            time = new AbsTime(*synth);

            memset(outL,0,sizeof(outL));
            memset(outR,0,sizeof(outR));

            fft = new FFTwrapper(BUF);
            //prepare the default settings
            params = new ADnoteParameters(*synth, fft, time, false);

            //sawtooth to make things a bit more interesting
            params->VoicePar[0].OscilGn->Pcurrentbasefunc = 3;

            params->GlobalPar.PFilterVelocityScale = 64;
            params->GlobalPar.GlobalFilter->basefreq = 5076.203125;

            controller = new Controller(*synth, time);

            //lets go with.... 50! as a nice note
            test_freq_log2 = log2f(440.0f) + (50.0 - 69.0f) / 12.0f;
        }

        void tearDown() {
            delete note;
            delete controller;
            delete fft;
            FFT_cleanup();
            delete time;
            delete synth;
            delete params;
        }

        void run_test(int a, int b, int c, int d, int e, int f, float values[4])
        {
            sprng(0);

            const int ADnote_unison_sizes[] =
            {1, 2, 3, 4, 5, 6, 8, 10, 12, 15, 20, 25, 30, 40, 50, 0};
            params->VoicePar[0].Unison_size = ADnote_unison_sizes[a];
            params->VoicePar[0].Unison_frequency_spread = b;
            params->VoicePar[0].Unison_stereo_spread    = c;
            params->VoicePar[0].Unison_vibratto         = d;
            params->VoicePar[0].Unison_vibratto_speed   = e;
            params->VoicePar[0].Unison_invert_phase     = f;

            SynthParams pars{memory, *controller, *synth, *time, 120, 0, test_freq_log2, false, prng()};
            note = new ADnote(params, pars);
            note->noteout(outL, outR);
            TS_ASSERT_DELTA(outL[80], values[0], 1.9e-5);
            printf("{%f,", outL[80]);
            note->noteout(outL, outR);
            TS_ASSERT_DELTA(outR[90], values[1], 1.9e-5);
            printf("%f,", outR[90]);
            note->noteout(outL, outR);
            TS_ASSERT_DELTA(outL[20], values[2], 1.9e-5);
            printf("%f,", outL[20]);
            note->noteout(outL, outR);
            TS_ASSERT_DELTA(outR[200], values[3], 1.9e-5);
            printf("%f},\n", outR[200]);
        }

        void testUnison() {
            sprng(0xbeef);

            float data[][4] = {
                {0.125972,0.029887,0.000000,0.138013},
                {-0.095414,-0.083965,-0.000000,0.009048},
                {-0.077587,-0.001760,-0.021463,-0.013995},
                {0.041240,-0.008561,-0.000000,-0.099298},
                {-0.098969,-0.048030,-0.000052,-0.087053},
                {0.104913,-0.081452,-0.017700,0.000978},
                {0.041270,0.003788,0.006064,0.002436},
                {-0.030791,-0.036072,-0.007964,-0.015141},
                {0.009218,0.015333,-0.007500,0.083076},
                {0.058909,0.064450,-0.002517,0.041595},
                {-0.007731,-0.009040,-0.068033,-0.016573},
                {-0.047286,-0.002355,-0.049196,0.016222},
                {0.014014,-0.002635,0.006542,0.050710},
                {-0.054877,-0.027135,0.040211,0.031927},
                {-0.048367,0.022010,0.018224,0.032846},
            };

            int freq_spread[15];
            int stereo_spread[15];
            int vibrato[15];
            int vibrato_speed[15];
            int inv_phase[15];
            for(int i=0; i<15; ++i)
            {
                freq_spread[i]    = prng()%0x7f;
                stereo_spread[i]  = prng()%0x7f;
                vibrato[i]        = prng()%0x7f;
                vibrato_speed[i]  = prng()%0x7f;
                inv_phase[i]      = prng()%5;
            }

            printf("\n");
            for(int i=0; i<15; ++i)
            {
                run_test(i, freq_spread[i], stereo_spread[i],
                        vibrato[i], vibrato_speed[i], inv_phase[i], data[i]);
            }
#if 0
            int sampleCount = 0;

            sampleCount += synth->buffersize;

            TS_ASSERT_DELTA(outL[255], 0.254609f, 0.0001f);

            note->noteout(outL, outR);
            sampleCount += synth->buffersize;
            TS_ASSERT_DELTA(outL[255], -0.102197f, 0.0001f);

            note->noteout(outL, outR);
            sampleCount += synth->buffersize;
            TS_ASSERT_DELTA(outL[255], -0.111422f, 0.0001f);

            note->noteout(outL, outR);
            sampleCount += synth->buffersize;
            TS_ASSERT_DELTA(outL[255], -0.021375f, 0.0001f);

            note->noteout(outL, outR);
            sampleCount += synth->buffersize;
            TS_ASSERT_DELTA(outL[255], 0.149882f, 0.0001f);
#endif
        }
};
