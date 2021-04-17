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
#include "../Misc/Master.h"
#include "../Misc/Util.h"
#include "../Misc/Allocator.h"
#include "../Synth/ADnote.h"
#include "../Params/Presets.h"
#include "../DSP/FFTwrapper.h"
#include "../Synth/LFO.h"
#include "../Params/LFOParams.h"
#include "../globals.h"
#include <rtosc/thread-link.h>

using namespace std;
using namespace zyn;

SYNTH_T *synth;


class AdNoteTest:public CxxTest::TestSuite
{
    public:
        rtosc::ThreadLink *tr;
        ADnote       *note;
        AbsTime      *time;
        FFTwrapper   *fft;
        ADnoteParameters *defaultPreset;
        Controller   *controller;
        Alloc         memory;
        float test_freq_log2;
        WatchManager *w;
        float *outR, *outL;

        LFO          *lfo;
        LFOParams    *lfop;
        int randval(int min, int max)
        {
            int ret = rand()%(1+max-min)+min;
            //printf("ret = %d (%d..%d)\n",ret, min,max);
            return ret;
        }

        void randomize_params(void) {
            lfop->Pintensity  = randval(0,255);
            lfop->Pstartphase = randval(0,255);
            lfop->Pcutoff     = randval(0,127);
            lfop->PLFOtype    = randval(0,6);
            lfop->Prandomness = randval(0,255);
            lfop->Pfreqrand   = randval(0,255);
            lfop->Pcontinous  = randval(0,1);
            lfop->Pstretch    = randval(0,255);
            lfop->fel         = (consumer_location_type_t) randval(1,2);
            
        }

        void run_lfo_randomtest(void)
        {
            lfo  = new LFO(*lfop, 440.0f, *time);
            for(int i=0; i<100; ++i) {
                float out = lfo->lfoout();
                switch(lfop->fel)
                {
                    case consumer_location_type_t::amp:
                        TS_ASSERT((-2.0f < out && out < 2.0f));
                        break;
                    case consumer_location_type_t::filter:
                        TS_ASSERT((-8.0f < out && out < 8.0f));
                        break;
                    case consumer_location_type_t::freq:
                    case consumer_location_type_t::unspecified:
                    default:
                        break;
                }
            }
        }



        void setUp() {
            //First the sensible settings and variables that have to be set:
            synth = new SYNTH_T;
            synth->buffersize = 256;
            //synth->alias();
            time  = new AbsTime(*synth);

            outL = new float[synth->buffersize];
            for(int i = 0; i < synth->buffersize; ++i)
                *(outL + i) = 0;
            outR = new float[synth->buffersize];
            for(int i = 0; i < synth->buffersize; ++i)
                *(outR + i) = 0;

            tr  = new rtosc::ThreadLink(1024,3);
            w   = new WatchManager(tr);

            fft = new FFTwrapper(synth->oscilsize);
            //prepare the default settings
            defaultPreset = new ADnoteParameters(*synth, fft, time, false);

            //Assert defaults
            TS_ASSERT(!defaultPreset->VoicePar[1].Enabled);

            XMLwrapper wrap;
            cout << string(SOURCE_DIR) + string("/guitar-adnote.xmz")
                 << endl;
            wrap.loadXMLfile(string(SOURCE_DIR)
                              + string("/guitar-adnote.xmz"));
            TS_ASSERT(wrap.enterbranch("MASTER"));
            TS_ASSERT(wrap.enterbranch("PART", 0));
            TS_ASSERT(wrap.enterbranch("INSTRUMENT"));
            TS_ASSERT(wrap.enterbranch("INSTRUMENT_KIT"));
            TS_ASSERT(wrap.enterbranch("INSTRUMENT_KIT_ITEM", 0));
            TS_ASSERT(wrap.enterbranch("ADD_SYNTH_PARAMETERS"));
            defaultPreset->getfromXML(wrap);
            //defaultPreset->defaults();

            //verify xml was loaded
            TS_ASSERT(defaultPreset->VoicePar[1].Enabled);



            controller = new Controller(*synth, time);

            //lets go with.... 50! as a nice note
            test_freq_log2 = log2f(440.0f) + (50.0 - 69.0f) / 12.0f;
            SynthParams pars{memory, *controller, *synth, *time, 120, 0, test_freq_log2, false, prng()};

            note = new ADnote(defaultPreset, pars,w);

        }

        void tearDown() {
            delete note;
            delete controller;
            delete defaultPreset;
            delete fft;
            delete [] outL;
            delete [] outR;
            FFT_cleanup();
            delete synth;
        }

        void testDefaults() {
            int sampleCount = 0;

//#define WRITE_OUTPUT

#ifdef WRITE_OUTPUT
            ofstream file("adnoteout", ios::out);
#endif
            note->noteout(outL, outR);
#ifdef WRITE_OUTPUT
            for(int i = 0; i < synth->buffersize; ++i)
                file << outL[i] << std::endl;

#endif
            sampleCount += synth->buffersize;
            TS_ASSERT_DELTA(outL[255], 0.1924f, 0.0001f);
            note->releasekey();

            TS_ASSERT(!tr->hasNext());
            w->add_watch("noteout/be4_mix");
            note->noteout(outL, outR);
            sampleCount += synth->buffersize;
            TS_ASSERT_DELTA(outL[255], -0.4717f, 0.0001f);
            w->tick();
            TS_ASSERT(tr->hasNext());

            note->noteout(outL, outR);
            sampleCount += synth->buffersize;
            w->tick();
            TS_ASSERT_DELTA(outL[255], 0.0646f, 0.0001f);

            note->noteout(outL, outR);
            sampleCount += synth->buffersize;
            TS_ASSERT_DELTA(outL[255], 0.1183f, 0.0001f);
            w->tick();
            note->noteout(outL, outR);
            sampleCount += synth->buffersize;
            TS_ASSERT_DELTA(outL[255], -0.1169f, 0.0001f);
            w->tick();

            TS_ASSERT(tr->hasNext());
            TS_ASSERT_EQUALS(string("noteout/be4_mix"), tr->read());
            TS_ASSERT(!tr->hasNext());

            note->noteout(outL, outR);
            sampleCount += synth->buffersize;


            while(!note->finished()) {
                note->noteout(outL, outR);
#ifdef WRITE_OUTPUT
                for(int i = 0; i < synth->buffersize; ++i)
                    file << outL[i] << std::endl;

#endif
                sampleCount += synth->buffersize;
            }
#ifdef WRITE_OUTPUT
            file.close();
#endif

            TS_ASSERT_EQUALS(sampleCount, 30208);

              lfop = new LFOParams();
            lfop->fel  = zyn::consumer_location_type_t::amp;
            lfop->freq = 2.0f;
            lfop->delay = 0.0f;
            for(int i=0; i<10000; ++i) {
                randomize_params();
                run_lfo_randomtest();
            }


        }

#define OUTPUT_PROFILE
#ifdef OUTPUT_PROFILE
        void testSpeed() {
            const int samps = 15000;

            int t_on = clock(); // timer before calling func
            for(int i = 0; i < samps; ++i)
                note->noteout(outL, outR);
            int t_off = clock(); // timer when func returns

            printf("AdNoteTest: %f seconds for %d Samples to be generated.\n",
                   (static_cast<float>(t_off - t_on)) / CLOCKS_PER_SEC, samps);
        }
#endif
};
