/*
  ZynAddSubFX - a software synthesizer

  JackEngine.cpp - Jack Driver
  Copyright (C) 2009 Alan Calvert
  Copyright (C) 2014 Mark McCurry

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.
*/

#include <iostream>

#include <jack/midiport.h>
#ifdef JACK_HAS_METADATA_API
# include <jack/metadata.h>
#endif // JACK_HAS_METADATA_API
#include "jack_osc.h"
#include <fcntl.h>
#include <sys/stat.h>
#include <cassert>
#include <cstring>
#include <chrono>
#include <unistd.h> // access()
#include <fstream> // std::istream


#include "Nio.h"
#include "OutMgr.h"
#include "InMgr.h"
#include "../Misc/Util.h"

#include "JackEngine.h"

extern char *instance_name;

namespace zyn {

using namespace std;

JackEngine::JackEngine(const SYNTH_T &synth)
    :AudioOut(synth), jackClient(NULL)
{
    name = "JACK";
    audio.jackSamplerate = 0;
    audio.jackNframes    = 0;
    for(int i = 0; i < 2; ++i) {
        audio.ports[i]     = NULL;
        audio.portBuffs[i] = NULL;
    }
    midi.inport = NULL;
    midi.jack_sync = false;
    osc.oscport = NULL;
    bufferindex = 0;
}

bool JackEngine::connectServer(string server)
{
    bool autostart_jack = true;
    if(jackClient)
        return true;

    string clientname = "zynaddsubfx";
    string postfix    = Nio::getPostfix();
    if(!postfix.empty())
        clientname += "_" + postfix;
    if(Nio::pidInClientName)
        clientname += "_" + os_pid_as_padded_string();

    jack_status_t jackstatus;
    bool use_server_name = server.size() && server.compare("default") != 0;
    jack_options_t jopts = (jack_options_t)
                           (((!instance_name
                              && use_server_name) ? JackServerName :
                             JackNullOption)
                            | ((autostart_jack) ? JackNullOption :
                               JackNoStartServer));

    if(instance_name)
        jackClient = jack_client_open(instance_name, jopts, &jackstatus);
    else {
        if(use_server_name)
            jackClient = jack_client_open(
                clientname.c_str(), jopts, &jackstatus,
                server.c_str());
        else
            jackClient = jack_client_open(
                clientname.c_str(), jopts, &jackstatus);
    }


    if(NULL != jackClient)
        return true;
    else
        cerr << "Error, failed to open jack client on server: " << server
             << " status " << jackstatus << endl;
    return false;
}

bool JackEngine::connectJack()
{
    connectServer("");
    if(NULL != jackClient) {
        setBufferSize(jack_get_buffer_size(jackClient));
        jack_set_error_function(_errorCallback);
        jack_set_info_function(_infoCallback);
        if(jack_set_buffer_size_callback(jackClient, _bufferSizeCallback, this))
            cerr << "Error setting the bufferSize callback" << endl;
        if((jack_set_xrun_callback(jackClient, _xrunCallback, this)))
            cerr << "Error setting jack xrun callback" << endl;
        if(jack_set_process_callback(jackClient, _processCallback, this)) {
            cerr << "Error, JackEngine failed to set process callback" << endl;
            return false;
        }
        if(jack_activate(jackClient)) {
            cerr << "Error, failed to activate jack client" << endl;
            return false;
        }

        return true;
    }
    else
        cerr << "Error, NULL jackClient through Start()" << endl;
    return false;
}

void JackEngine::disconnectJack()
{
    if(jackClient) {
        cout << "Deactivating and closing JACK client" << endl;

        jack_deactivate(jackClient);
        jack_client_close(jackClient);
        jackClient = NULL;
    }
}

bool JackEngine::Start()
{
    return openMidi() && openAudio();
}

void JackEngine::Stop()
{
    stopMidi();
    stopAudio();
}

void JackEngine::setMidiEn(bool nval)
{
    if(nval)
        openMidi();
    else
        stopMidi();
}

bool JackEngine::getMidiEn() const
{
    return midi.inport;
}

void JackEngine::setAudioEn(bool nval)
{
    if(nval)
        openAudio();
    else
        stopAudio();
}

bool JackEngine::getAudioEn() const
{
    return audio.ports[0];
}

bool JackEngine::openAudio()
{
    if(getAudioEn())
        return true;

    if(!getMidiEn())
        if(!connectJack())
            return false;


    const char *portnames[] = { "out_1", "out_2" };
    for(int port = 0; port < 2; ++port)
        audio.ports[port] = jack_port_register(
            jackClient,
            portnames[port],
            JACK_DEFAULT_AUDIO_TYPE,
            JackPortIsOutput
            | JackPortIsTerminal,
            0);
    if((NULL != audio.ports[0]) && (NULL != audio.ports[1])) {
        audio.jackSamplerate = jack_get_sample_rate(jackClient);
        audio.jackNframes    = jack_get_buffer_size(jackClient);
        samplerate = audio.jackSamplerate;
        bufferSize = audio.jackNframes;


        //Attempt to autoConnect when specified
        if(Nio::autoConnect) {
            const char **outPorts = jack_get_ports(
                jackClient,
                NULL,
                NULL,
                JackPortIsPhysical
                | JackPortIsInput);
            if(outPorts != NULL) {
                //Verify that stereo is available
                assert(outPorts[0]);
                assert(outPorts[1]);

                //Connect to physical outputs
                jack_connect(jackClient, jack_port_name(
                                 audio.ports[0]), outPorts[0]);
                jack_connect(jackClient, jack_port_name(
                                 audio.ports[1]), outPorts[1]);
            }
            else
                cerr << "Warning, No outputs to autoconnect to" << endl;
        }
        midi.jack_sync = true;
        osc.oscport = jack_port_register(jackClient, "osc",
                JACK_DEFAULT_OSC_TYPE, JackPortIsInput, 0);
#ifdef JACK_HAS_METADATA_API
        jack_uuid_t uuid = jack_port_uuid(osc.oscport);
        jack_set_property(jackClient, uuid, "http://jackaudio.org/metadata/event-types", JACK_EVENT_TYPE__OSC, "text/plain");
#endif // JACK_HAS_METADATA_API
        return true;
    }
    else
        cerr << "Error, failed to register jack audio ports" << endl;
    midi.jack_sync = false;
    return false;
}

void JackEngine::stopAudio()
{
    for(int i = 0; i < 2; ++i) {
        jack_port_t *port = audio.ports[i];
        audio.ports[i] = NULL;
        if(jackClient != NULL && NULL != port)
            jack_port_unregister(jackClient, port);
    }
    midi.jack_sync = false;
    if(osc.oscport) {
       if (jackClient != NULL) {
               jack_port_unregister(jackClient, osc.oscport);
#ifdef JACK_HAS_METADATA_API
               jack_uuid_t uuid = jack_port_uuid(osc.oscport);
               jack_remove_property(jackClient, uuid, "http://jackaudio.org/metadata/event-types");
#endif // JACK_HAS_METADATA_API
       }
    }
    if(!getMidiEn())
        disconnectJack();
}

bool JackEngine::openMidi()
{
    if(getMidiEn())
        return true;
    if(!getAudioEn())
        if(!connectJack())
            return false;

    midi.inport = jack_port_register(jackClient, "midi_input",
                                     JACK_DEFAULT_MIDI_TYPE,
                                     JackPortIsInput | JackPortIsTerminal, 0);
    return midi.inport;
}

void JackEngine::stopMidi()
{
    jack_port_t *port = midi.inport;
    midi.inport = NULL;
    if(port)
        jack_port_unregister(jackClient, port);

    if(!getAudioEn())
        disconnectJack();
}

int JackEngine::clientId()
{
    if(NULL != jackClient)
        return (long)jack_client_thread_id(jackClient);
    else
        return -1;
}

string JackEngine::clientName()
{
    if(NULL != jackClient)
        return string(jack_get_client_name(jackClient));
    else
        cerr << "Error, clientName() with null jackClient" << endl;
    return string("Oh, yoshimi :-(");
}

int JackEngine::_processCallback(jack_nframes_t nframes, void *arg)
{
    return static_cast<JackEngine *>(arg)->processCallback(nframes);
}

int JackEngine::processCallback(jack_nframes_t nframes)
{
    bool okaudio = true;
    bufferindex++;
    handleMidi(nframes);
    if((NULL != audio.ports[0]) && (NULL != audio.ports[1]))
        okaudio = processAudio(nframes);
    return okaudio ? 0 : -1;
}

bool JackEngine::processAudio(jack_nframes_t nframes)
{
    //handle rt osc events first
    void   *oscport    = jack_port_get_buffer(osc.oscport, nframes);
    size_t osc_packets = jack_osc_get_event_count(oscport);

    for(size_t i = 0; i < osc_packets; ++i) {
        jack_osc_event_t event;
        if(jack_osc_event_get(&event, oscport, i))
            continue;
        if(*event.buffer!='/') //Bundles are unhandled
            continue;
        //TODO validate message length
        OutMgr::getInstance().applyOscEventRt((char*)event.buffer);
    }

    for(int port = 0; port < 2; ++port) {
        audio.portBuffs[port] =
            (jsample_t *)jack_port_get_buffer(audio.ports[port], nframes);
        if(NULL == audio.portBuffs[port]) {
            cerr << "Error, failed to get jack audio port buffer: "
                 << port << endl;
            return false;
        }
    }
    jack_nframes_t frame = jack_last_frame_time(jackClient);
    tstampaudio = jack_frames_to_time (jackClient,frame)*1000;

    Stereo<float *> smp = getNext();

    //Assumes size of smp.l == nframes
    memcpy(audio.portBuffs[0], smp.l, bufferSize * sizeof(float));
    memcpy(audio.portBuffs[1], smp.r, bufferSize * sizeof(float));
    return true;
}

int JackEngine::_xrunCallback(void *)
{
    cerr << "Jack reports xrun" << endl;
    return 0;
}

void JackEngine::_errorCallback(const char *msg)
{
    cerr << "Jack reports error: " << msg << endl;
}

void JackEngine::_infoCallback(const char *msg)
{
    cerr << "Jack info message: " << msg << endl;
}

int JackEngine::_bufferSizeCallback(jack_nframes_t nframes, void *arg)
{
    return static_cast<JackEngine *>(arg)->bufferSizeCallback(nframes);
}

int JackEngine::bufferSizeCallback(jack_nframes_t nframes)
{
    cerr << "Jack buffer resized" << endl;
    setBufferSize(nframes);
    return 0;
}

void JackEngine::handleMidi(unsigned long frames)
{
    if(!midi.inport)
        return;
    void *midi_buf = jack_port_get_buffer(midi.inport, frames);
    jack_midi_event_t jack_midi_event;
    jack_nframes_t    event_index = 0;
    unsigned char     buf[3];

    while(jack_midi_event_get(&jack_midi_event, midi_buf,
                              event_index++) == 0) 
    {
        //~ memset(buf, 0, sizeof(buf));
        //~ memcpy(buf, jack_midi_event.buffer,
            //~ std::min(sizeof(buf), jack_midi_event.size));

        
        jack_nframes_t time = jack_midi_event.time;    
        jack_nframes_t nframes = jack_last_frame_time(jackClient);
        jack_time_t micros = jack_frames_to_time (jackClient,nframes+time);
        
        for (size_t x = 0; x < jack_midi_event.size; x += 3) 
        {
            size_t y = jack_midi_event.size - x;
            if (y >= 3) {
                memcpy(buf, (uint8_t *)jack_midi_event.buffer + x, 3);
            } else {
                memset(buf, 0, sizeof(buf));
                memcpy(buf, (uint8_t *)jack_midi_event.buffer + x, y);
            }
            /* make sure the values are within range */
            buf[1] &= 0x7F;
            buf[2] &= 0x7F;
            midiProcess(buf[0], buf[1], buf[2], time, micros*1000);
        }
    }
}

}
