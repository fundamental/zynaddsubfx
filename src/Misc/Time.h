/*
  ZynAddSubFX - a software synthesizer

  Time.h - Frame Tracker
  Copyright (C) 2016 Mark McCurry

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.
*/
#pragma once
#include <stdint.h>
#include "../globals.h"

namespace zyn {

class AbsTime
{
    public:
        AbsTime(const SYNTH_T &synth)
            :frames(0),
            s(synth){};
        void operator++(){++frames;};
        void operator++(int){frames++;};
        int64_t time() const {return frames;};
        unsigned long tStamp; // time stamp of the current buffer
        unsigned long tRef; // time stamp midi message
        int bpm;
        int64_t time2nanos(int64_t t) {return (t-frames)*framesPerSec()*1000000000.0f;}; // referenceTime
        float dt() const { return s.dt(); }
        float framesPerSec() const {return 1.0f/s.dt();}
        int   framesPerBuffer() const {return s.buffersize;}
        float bufferDuration() const {return s.buffersize*s.dt();}
    private:
        int64_t frames;
        const SYNTH_T &s;
};

//Marker for an event relative to some position of the absolute timer
class RelTime
{
    public:
        RelTime(const AbsTime &t_, float sec)
            :t(t_)
        {
            //Calculate time of event
            double deltaFrames = sec*t.framesPerSec();
            int64_t tmp = (int64_t)deltaFrames;
            frame = t.time() + tmp;
            sample = (int32_t)(t.framesPerBuffer()*(deltaFrames-tmp));
        }
        bool inThisFrame() {return t.time() == frame;};
        bool inPast() {return t.time() > frame;}
        bool inFuture() {return t.time() < frame;}
    private:
        int64_t frame;
        int32_t sample;
        const AbsTime &t;
};

}
