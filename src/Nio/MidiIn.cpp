/*
  ZynAddSubFX - a software synthesizer

  MidiIn.C - This class is inherited by all the Midi input classes
  Copyright (C) 2002-2005 Nasca Octavian Paul
  Author: Nasca Octavian Paul

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.
*/

#include "MidiIn.h"
#include "../globals.h"
#include "InMgr.h"
#include <string.h>

namespace zyn {

MidiIn::MidiIn()
{
    sysex_offset = 0;
    memset(sysex_data, 0, sizeof(sysex_data));
}

uint8_t MidiIn::midiSysEx(unsigned char data)
{
    if (data & 0x80) {
        if (data == 0xF0) {
            sysex_offset = 0; /* begin */
        } else if (data == 0xF7) {
                return (2); /* end */
        } else {
                return (1); /* error */
        }
    } else if (sysex_offset >= sizeof(sysex_data)) {
        return (1); /* error */
    }
    sysex_data[sysex_offset++] = data;
    return (0);
}

void MidiIn::midiProcess(unsigned char head,
                         unsigned char num,
                         unsigned char value, 
                         unsigned long time, 
                         unsigned long nanos,
                         unsigned long nframes)
{

    unsigned char chan = head & 0x0f;
    unsigned char type = head & 0xf0;
    
    MidiEvent ev;
    ev.channel = 0;
    ev.type = 0;
    ev.num = 0;
    ev.value = 0;
    ev.time = time;
    ev.nanos = nanos;
    //~ printf("midiProcess(h:%d, n:%d, v:%d, m:%lu)\n", head, num, value, nanos);
    /* SYSEX handling */
    if (head == 0xF0 || sysex_offset != 0) {
        uint8_t status = 0;

        status |= midiSysEx(type);
        status |= midiSysEx(num);
        status |= midiSysEx(value);

        if (status & 1) {
            /* error parsing SYSEX */
            sysex_offset = 0;
        } else if (status & 2) {
            /* message complete */

            if (sysex_offset >= 10 &&
                sysex_data[1] == 0x0A &&
                sysex_data[2] == 0x55) {
                ev.channel = sysex_data[3] & 0x0F;
                ev.log2_freq = sysex_data[6] +
                  (sysex_data[7] / (128.0f)) +
                  (sysex_data[8] / (128.0f * 128.0f)) +
                  (sysex_data[9] / (128.0f * 128.0f * 128.0f));

                switch (sysex_data[3] >> 4) {
                case 0: /* Note ON */
                    ev.type = M_FLOAT_NOTE;
                    ev.num = sysex_data[4];
                    ev.value = sysex_data[5];
                    ev.log2_freq /= 12.0f;
                    break;
                case 1: /* Pressure, Aftertouch */
                    ev.type = M_FLOAT_CTRL;
                    ev.num = C_aftertouch;
                    ev.value = sysex_data[4];
                    break;
                case 2: /* Controller */
                    ev.type = M_FLOAT_CTRL;
                    ev.num = sysex_data[5];
                    ev.value = sysex_data[4];
                    break;
                case 3: /* Absolute pitch */
                    ev.type = M_FLOAT_CTRL;
                    ev.num = C_pitch;
                    ev.value = sysex_data[4];
                    ev.log2_freq /= 12.0f;
                    break;
                default:
                    return;
                }
                InMgr::getInstance().putEvent(ev);
            }
            
            if (sysex_offset >= 7 &&
                sysex_data[1] == 0x7F &&
                sysex_data[3] == 0x03) {
                    ev.type = M_TIME_SIG;
                    ev.num = sysex_data[5]; // numerator
                    ev.value = sysex_data[6]; // denominator
                    InMgr::getInstance().putEvent(ev);
                }
            
            return; /* message complete */
        } else {
            return; /* wait for more data */
        }
    }
    //
    switch(type) {
        case 0x80: //Note Off
            ev.type    = M_NOTE;
            ev.channel = chan;
            ev.num     = num;
            ev.value   = 0;
            InMgr::getInstance().putEvent(ev);
            break;
        case 0x90: //Note On
            ev.type    = M_NOTE;
            ev.channel = chan;
            ev.num     = num;
            ev.value   = value;
            InMgr::getInstance().putEvent(ev);
            break;
        case 0xA0: /* pressure, aftertouch */
            ev.type    = M_PRESSURE;
            ev.channel = chan;
            ev.num     = num;
            ev.value   = value;
            InMgr::getInstance().putEvent(ev);
            break;
        case 0xb0: //Controller
            ev.type    = M_CONTROLLER;
            ev.channel = chan;
            ev.num     = num;
            ev.value   = value;
            InMgr::getInstance().putEvent(ev);
            break;
        case 0xc0: //Program Change
            ev.type    = M_PGMCHANGE;
            ev.channel = chan;
            ev.num     = num;
            ev.value   = 0;
            InMgr::getInstance().putEvent(ev);
            break;
        case 0xe0: //Pitch Wheel
            ev.type    = M_CONTROLLER;
            ev.channel = chan;
            ev.num     = C_pitchwheel;
            ev.value   = (num + value * (int) 128) - 8192;
            InMgr::getInstance().putEvent(ev);
            break;
        case 0x79: // reset to power-on state
            ev.type    = M_CONTROLLER;
            ev.channel = chan;
            ev.num     = C_resetallcontrollers;
            ev.value   = 0;
            InMgr::getInstance().putEvent(ev);
                break;
    }
    
    /* MTC handling */
    int subtype = (value >> 4) & 0xf;
    int val  = (value & 0xf);
    ev.channel = 0;

    switch(head) {
        case 0xF1: /* M_TC - MIDI Time Code Quarter Frame */
            switch (subtype)
            {
            case 0:
            {
                frameLo = val;
                break;
            }
            case 1:
            {
                frame = frameLo + (val * 16);
                break;
            }
            case 2:
            {
                secondLo = val;
                break;
            }
            case 3:
            {
                second = secondLo + (val * 16);
                break;
            }
            case 4:
            {
                minuteLo = val;
                break;
            }
            case 5:
            {
                minute = minuteLo + (val * 16);
                break;
            }
            case 6:
            {
                hourLo = (val & 31);
                break;
            }
            case 7:
            {
                hour = hourLo + ((val * 16) & 31);
                if(frame==23) {
                    ev.type  = M_TC;
                    ev.num = 0;
                    ev.value = second+60*minute+3600*hour;
                    InMgr::getInstance().putEvent(ev);
                }
                break;
            }
            }
            break;

        case 0xF2: /* M_SPP - Song Position Pointer*/
            ev.type  = M_SPP;
            ev.value = num + 128 * value;
            InMgr::getInstance().putEvent(ev);
            break;
        case 0xF8: /* M_CLOCK */
            ev.type  = M_CLOCK;
            InMgr::getInstance().putEvent(ev);
            break;
        case 0xFA: /* M_START */
            ev.type  = M_START;
            InMgr::getInstance().putEvent(ev);
            break;
        case 0xFB: /* M_CONTINUE */
            ev.type  = M_CONTINUE;
            InMgr::getInstance().putEvent(ev);
            break;
        case 0xFC: /* M_STOP */
            ev.type  = M_STOP;
            InMgr::getInstance().putEvent(ev);
            break;
    }

}

}
