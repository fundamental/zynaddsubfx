/*
  ZynAddSubFX - a software synthesizer

  MidiIn.h - This class is inherited by all the Midi input classes
  Copyright (C) 2002-2005 Nasca Octavian Paul
  Copyright (C) 2009-2010 Mark McCurry
  Author: Nasca Octavian Paula
          Mark McCurry

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.
*/

#ifndef MIDI_IN_H
#define MIDI_IN_H

#include "Engine.h"

namespace zyn {

/**This class is inherited by all the Midi input classes*/
class MidiIn:public virtual Engine
{
    public:
        MidiIn();

        /**Enables or disables driver based upon value*/
        virtual void setMidiEn(bool nval) = 0;
        /**Returns if driver is initialized*/
        virtual bool getMidiEn() const = 0;
        void midiProcess(unsigned char head,
                         unsigned char num,
                         unsigned char value, 
                         unsigned long time = 0, 
                         unsigned long nanos = 0, 
                         unsigned long nframes = 0);

    private:
        uint8_t midiSysEx(unsigned char data);
        uint8_t sysex_offset;
        uint8_t sysex_data[64];
        int frameLo, frame, secondLo, second, minuteLo, minute, hourLo, hour;
        
        
};

}

#endif
