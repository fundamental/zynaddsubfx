/*
  ZynAddSubFX - a software synthesizer

  InMgr.cpp - MIDI Input Manager
  Copyright (C) 2016 Mark McCurry

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.
*/
#include "InMgr.h"
#include "MidiIn.h"
#include "EngineMgr.h"
#include "../Misc/Master.h"
#include "../Misc/Part.h"
#include "../Misc/MiddleWare.h"
#include <rtosc/thread-link.h>
#include <iostream>
using namespace std;

extern zyn::MiddleWare *middleware;

namespace zyn {

ostream &operator<<(ostream &out, const MidiEvent &ev)
{
    switch(ev.type) {
        case M_NOTE:
            out << "MidiNote: note(" << ev.num << ")\n"
            << "          channel(" << ev.channel << ")\n"
            << "          velocity(" << ev.value << ")";
            break;

        case M_FLOAT_NOTE:
            out << "MidiNote: note(" << ev.num << ")\n"
            << "          channel(" << ev.channel << ")\n"
            << "          velocity(" << ev.value << ")\n"
            << "          log2_freq(" << ev.log2_freq << ")";
            break;

        case M_CONTROLLER:
            out << "MidiCtl: controller(" << ev.num << ")\n"
            << "         channel(" << ev.channel << ")\n"
            << "         value(" << ev.value << ")";
            break;

        case M_FLOAT_CTRL:
            out << "MidiNote: controller(" << ev.num << ")\n"
            << "          channel(" << ev.channel << ")\n"
            << "          note(" << ev.value << ")\n"
            << "          log2_value(" << ev.log2_freq << ")";
            break;

        case M_PGMCHANGE:
            out << "PgmChange: program(" << ev.num << ")\n"
            << "           channel(" << ev.channel << ")";
            break;

        case M_CLOCK:
            out << "M_CLOCK\n";
            break;
    }

    return out;
}

MidiEvent::MidiEvent()
    :channel(0), type(0), num(0), value(0), time(0)
{}

InMgr &InMgr::getInstance()
{
    static InMgr instance;
    return instance;
}

InMgr::InMgr()
    :queue(256), master(NULL)
{
    current = NULL;
    work.init(PTHREAD_PROCESS_PRIVATE, 0);
}

InMgr::~InMgr()
{
    //lets stop the consumer thread
}

void InMgr::putEvent(MidiEvent ev)
{
    if(queue.push(ev)) //check for error
        cerr << "ERROR: MIDI ringbuffer is FULL" << endl;
    else
        work.post();
}

void InMgr::flush(unsigned frameStart, unsigned frameStop)
{
    MidiEvent ev;
    while(!work.trywait()) {
        queue.peak(ev);
        if(ev.time < (int)frameStart || ev.time > (int)frameStop) {
            //Back out of transaction
            work.post();
            //~ printf("%lu vs [%d..%d]\n", ev.time, frameStart, frameStop);
            break;
        }
        queue.pop(ev);

        switch(ev.type) {
            case M_NOTE:
                master->noteOn(ev.channel, ev.num, ev.value);
                break;

            case M_FLOAT_NOTE:
                master->noteOn(ev.channel, ev.num, ev.value, ev.log2_freq);
                break;

            case M_CONTROLLER:
                if(ev.num == C_bankselectmsb) {        // Change current bank
                    middleware->spawnMaster()->bToU->write("/forward", "");
                    middleware->spawnMaster()->bToU->write("/bank/msb", "i", ev.value);
                    middleware->spawnMaster()->bToU->write("/bank/bank_select", "i", ev.value);

                } else if(ev.num == C_bankselectlsb)  {// Change current bank (LSB)
                    middleware->spawnMaster()->bToU->write("/forward", "");
                    middleware->spawnMaster()->bToU->write("/bank/lsb", "i", ev.value);
                } else
                    master->setController(ev.channel, ev.num, ev.value);
                break;

            case M_FLOAT_CTRL:
                master->setController(ev.channel, ev.num, ev.value, ev.log2_freq);
                break;

            case M_PGMCHANGE:
                for(int i=0; i < NUM_MIDI_PARTS; ++i) {
                    //set the program of the parts assigned to the midi channel
                    if(master->part[i]->Prcvchn == ev.channel) {
                        middleware->pendingSetProgram(i, ev.num);
                    }
                }
                break;

            case M_PRESSURE:
                master->polyphonicAftertouch(ev.channel, ev.num, ev.value);
                break;
            
            case M_CLOCK: // needed to determine bpm
                master->midiClock(ev.nanos);
                break;
            
            case M_TC: // suited to determine reference time for calculating phase
                master->midiTcSync(ev.nanos, ev.value);
                break;

            case M_START: // suited to determine reference time for calculating phase
                printf("M_START\n");
                if (!isPlaying) {
                    master->midiSppSync(ev.nanos, 0);
                    isPlaying = true;
                }
                break;

            case M_CONTINUE: // suited to determine reference time for calculating phase
                printf("M_CONTINUE\n");
                if (!isPlaying) {
                    //~ master->midiSppSync(ev.nanos, 0);
                    isPlaying = true;
                }
                break;

            case M_STOP: 
                isPlaying = false;
                printf("M_STOP\n");
                break;

            case M_TIME_SIG: // suited to determine reference time for calculating phase
                master->setSignature(ev.num, ev.value);
                break;

            case M_SPP: // suited to determine reference time for calculating phase
                printf("M_SPP\n");
                master->midiSppSync(ev.nanos, ev.value);
                break;
        }
    }
}

bool InMgr::empty(void) const
{
    int semvalue = work.getvalue();
    return semvalue <= 0;
}

bool InMgr::setSource(string name)
{
    MidiIn *src = getIn(name);

    if(!src)
        return false;

    if(current)
        current->setMidiEn(false);
    current = src;
    current->setMidiEn(true);

    bool success = current->getMidiEn();

    //Keep system in a valid state (aka with a running driver)
    if(!success)
        (current = getIn("NULL"))->setMidiEn(true);

    return success;
}

string InMgr::getSource() const
{
    if(current)
        return current->name;
    else
        return "ERROR";
}

MidiIn *InMgr::getIn(string name)
{
    EngineMgr &eng = EngineMgr::getInstance(NULL);
    return dynamic_cast<MidiIn *>(eng.getEng(name));
}

void InMgr::setMaster(Master *master_)
{
    master = master_;
}

}
