/*
  ZynAddSubFX - a software synthesizer

  InMgr.h - MIDI Input Manager
  Copyright (C) 2016 Mark McCurry

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.
*/
#ifndef INMGR_H
#define INMGR_H

#include <string>
#include "ZynSema.h"
#include "SafeQueue.h"

namespace zyn {

enum midi_type {
    M_NOTE = 1, // for note
    M_CONTROLLER = 2, // for controller
    M_PGMCHANGE  = 3, // for program change
    M_PRESSURE   = 4, // for polyphonic aftertouch
    M_FLOAT_NOTE = 5, // for floating point note
    M_FLOAT_CTRL = 6, // for floating point controller
    M_TC         = 7, // for MIDI Timecode
    M_SPP        = 8, // for Song Position Pointer
    M_CLOCK      = 9, // for MIDI Clock Pulse
    M_START      = 10,// for MIDI Transport Start
    M_STOP       = 11,// for MIDI Transport Stop
    M_CONTINUE   = 12,// for MIDI Transport Continue
    M_TIME_SIG   = 13 // for MIDI Time Signature
};

struct MidiEvent {
    MidiEvent();
    int channel; //the midi channel for the event
    int type;    //type=1 for note, type=2 for controller
    int num;     //note, controller or program number
    int value;   //velocity or controller value
    unsigned long time;    //absolute frame of the midi event
    unsigned long nanos;    //nanoseconds of the midi event
    unsigned long taudio;    //nanoseconds of the midi event
    float log2_freq;   //type=5,6 for logarithmic representation of note/parameter
};

//super simple class to manage the inputs
class InMgr
{
    public:
        static InMgr &getInstance();
        ~InMgr();

        void putEvent(MidiEvent ev);

        /**Flush the Midi Queue*/
        void flush(unsigned frameStart, unsigned frameStop);

        bool empty() const;

        bool setSource(std::string name);

        std::string getSource() const;

        void setMaster(class Master *master);

        friend class EngineMgr;
    private:
        InMgr();
        class MidiIn *getIn(std::string name);
        SafeQueue<MidiEvent> queue;
        mutable ZynSema work;
        class MidiIn * current;

        /**the link to the rest of zyn*/
        class Master *master;
        bool isPlaying = false;
};

}

#endif
