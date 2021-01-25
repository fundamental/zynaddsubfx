/*
  ZynAddSubFX - a software synthesizer

  Note.h - Abstract Base Class for synthesizers
  Copyright (C) 2010-2010 Mark McCurry
  Author: Mark McCurry

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.
*/
#ifndef SYNTH_NOTE_H
#define SYNTH_NOTE_H
#include "../globals.h"
#include "../Misc/Util.h"
#include "../Containers/NotePool.h"

namespace zyn {

class Allocator;
class Controller;
struct SynthParams
{
    Allocator &memory;   //Memory Allocator for the Note to use
    const Controller &ctl;
    const SYNTH_T    &synth;
    const AbsTime    &time;
    float     velocity;  //Velocity of the Note
    bool      portamento;//True if portamento is used for this note
    float     note_log2_freq; //Floating point value of the note
    bool      quiet;     //Initial output condition for legato notes
    prng_t    seed;      //Random seed
    
};

struct LegatoParams
{
    float velocity;
    bool portamento;
    float note_log2_freq; //Floating point value of the note
    bool externcall;
    prng_t seed;
};

class SynthNote
{
    public:
        SynthNote(const SynthParams &pars);
        virtual ~SynthNote() {}

        /**Compute Output Samples
         * @return 0 if note is finished*/
        virtual int noteout(float *outl, float *outr) = 0;

        //TODO fix this spelling error [noisey commit]
        /**Release the key for the note and start release portion of envelopes.*/
        virtual void releasekey() = 0;

        /**Return if note is finished.
         * @return finished=1 unfinished=0*/
        virtual bool finished() const = 0;

        /**Make a note die off next buffer compute*/
        virtual void entomb(void) = 0;

        virtual void legatonote(const LegatoParams &pars) = 0;

        virtual SynthNote *cloneLegato(void) = 0;

        /* For polyphonic aftertouch needed */
        void setVelocity(float velocity_);

        /* For per-note pitch */
        void setPitch(float log2_freq_);

        /* For per-note filter cutoff */
        void setFilterCutoff(float);
        float getFilterCutoffRelFreq(void);

        /* Random numbers with own seed */
        float getRandomFloat();
        prng_t getRandomUint();
        
        /* For midi synced LFOs*/
        void setBPM(int bpm_, int phase_);

        //Realtime Safe Memory Allocator For notes
        class Allocator  &memory;
    protected:
        // Legato transitions
        class Legato
        {
            public:
                Legato(const SYNTH_T &synth_, float vel, int port,
                       float note_log2_freq, bool quiet, prng_t seed);

                void apply(SynthNote &note, float *outl, float *outr);
                int update(const LegatoParams &pars);

            private:
                bool      silent;
                float     lastfreq_log2;
                LegatoMsg msg;
                int       decounter;
                struct { // Fade In/Out vars
                    int   length;
                    float m, step;
                } fade;
            public:
                struct { // Note parameters
                    float  freq, vel;
                    bool   portamento;
                    float  note_log2_freq;
                    prng_t seed;
                } param;
                const SYNTH_T &synth;

            public: /* Some get routines for legatonote calls (aftertouch feature)*/
                float getFreq() {return param.freq; }
                float getVelocity() {return param.vel; }
                bool  getPortamento() {return param.portamento; }
                float getNoteLog2Freq() {return param.note_log2_freq; }
                prng_t getSeed() {return param.seed;}
                void setSilent(bool silent_) {silent = silent_; }
                void setDecounter(int decounter_) {decounter = decounter_; }
        } legato;

        prng_t initial_seed;
        prng_t current_prng_state;
        const Controller &ctl;
        const SYNTH_T    &synth;
        const AbsTime    &time;
        WatchManager     *wm;
        smooth_float     filtercutoff_relfreq;
        int bpm, phase;
};

}

#endif
