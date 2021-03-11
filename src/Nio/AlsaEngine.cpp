/*
  ZynAddSubFX - a software synthesizer

  AlsaEngine.cpp - ALSA Driver
  Copyright (C) 2009 Alan Calvert
  Copyright (C) 2014 Mark McCurry

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.
*/

#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <chrono>
#include <poll.h>

#include "../Misc/Util.h"
#include "../Misc/Config.h"
#include "InMgr.h"
#include "AlsaEngine.h"
#include "Compressor.h"
#include "Nio.h"

using namespace std;

namespace zyn {

AlsaEngine::AlsaEngine(const SYNTH_T &synth)
    :AudioOut(synth)
{
    audio.buffer = new short[synth.buffersize * 2];
    name = "ALSA";
    audio.handle = NULL;
    audio.peaks[0] = 0;

    midi.handle  = NULL;
    midi.alsaId  = -1;
    midi.pThread = 0;
    
    framecount = 0;
}

AlsaEngine::~AlsaEngine()
{
    Stop();
    delete[] audio.buffer;
}

void *AlsaEngine::_AudioThread(void *arg)
{
    return (static_cast<AlsaEngine *>(arg))->AudioThread();
}

void *AlsaEngine::AudioThread()
{
    set_realtime();
    return processAudio();
}

bool AlsaEngine::Start()
{
    return openAudio() && openMidi();
}

void AlsaEngine::Stop()
{
    if(getMidiEn())
        setMidiEn(false);
    if(getAudioEn())
        setAudioEn(false);
    snd_config_update_free_global();
}

void AlsaEngine::setMidiEn(bool nval)
{
    if(nval)
        openMidi();
    else
        stopMidi();
}

bool AlsaEngine::getMidiEn() const
{
    return midi.handle;
}

void AlsaEngine::setAudioEn(bool nval)
{
    if(nval)
        openAudio();
    else
        stopAudio();
}

bool AlsaEngine::getAudioEn() const
{
    return audio.handle;
}

void *AlsaEngine::_MidiThread(void *arg)
{
    return static_cast<AlsaEngine *>(arg)->MidiThread();
}


void *AlsaEngine::MidiThread(void)
{
    snd_seq_event_t *event; 
    //~ unsigned char     buf[3];
    MidiEvent ev = {};
    struct pollfd pfd[4 /* XXX 1 should be enough */];
    int error;

    set_realtime();
    while(1) {
        if(midi.exiting)
            break;
        error = snd_seq_event_input(midi.handle, &event);
        if (error < 0) {
            if(error != -EAGAIN && error != -EINTR)
                break;
            error = snd_seq_poll_descriptors(midi.handle, pfd, 4, POLLIN);
            if(error <= 0)
                break;
            error = poll(pfd, error, 1000 /* ms */);
            if(error < 0 &&
               errno != EAGAIN && errno != EINTR)
            break;
            continue;
        }
        //ensure ev is empty
        ev.channel = 0;
        ev.num     = 0;
        ev.value   = 0;
        ev.type    = 0;
        ev.time    = 0;
        
        tmidi = (unsigned long)event->time.time.tv_sec * 1000000UL + (unsigned long)event->time.time.tv_nsec/1000UL;
        
        
        ev.nanos = tmidi + tstamptrigger;
        //~ ev.taudio = htstampaudio.tv_sec * 1000000000 + htstampaudio.tv_nsec; // constant
        
        printf("\nevent->type =  %d\n", event->type);
        printf("tmidi + ttrigger %lu + %lu = (ev.nanos) %lu\n", tmidi, tstamptrigger, ev.nanos);

        if(!event)
            continue;
            
       //~ snd_midi_event_t* midiParser;
        
        //~ printf("event->type: %d\n", event->type);
        switch(event->type) {
            case SND_SEQ_EVENT_NOTEON:
                ev.type    = M_NOTE;
                ev.channel = event->data.note.channel;
                ev.num     = event->data.note.note;
                ev.value   = event->data.note.velocity;
                InMgr::getInstance().putEvent(ev);
                break;

            case SND_SEQ_EVENT_NOTEOFF:
                ev.type    = M_NOTE;
                ev.channel = event->data.note.channel;
                ev.num     = event->data.note.note;
                ev.value   = 0;
                InMgr::getInstance().putEvent(ev);
                break;

            case SND_SEQ_EVENT_KEYPRESS:
                ev.type    = M_PRESSURE;
                ev.channel = event->data.note.channel;
                ev.num     = event->data.note.note;
                ev.value   = event->data.note.velocity;
                InMgr::getInstance().putEvent(ev);
                break;

            case SND_SEQ_EVENT_PITCHBEND:
                ev.type    = M_CONTROLLER;
                ev.channel = event->data.control.channel;
                ev.num     = C_pitchwheel;
                ev.value   = event->data.control.value;
                InMgr::getInstance().putEvent(ev);
                break;

            case SND_SEQ_EVENT_CONTROLLER:
                ev.type    = M_CONTROLLER;
                ev.channel = event->data.control.channel;
                ev.num     = event->data.control.param;
                ev.value   = event->data.control.value;
                InMgr::getInstance().putEvent(ev);
                break;

            case SND_SEQ_EVENT_PGMCHANGE:
                ev.type    = M_PGMCHANGE;
                ev.channel = event->data.control.channel;
                ev.num     = event->data.control.value;
                InMgr::getInstance().putEvent(ev);
                break;

            case SND_SEQ_EVENT_RESET: // reset to power-on state
                ev.type    = M_CONTROLLER;
                ev.channel = event->data.control.channel;
                ev.num     = C_resetallcontrollers;
                ev.value   = 0;
                InMgr::getInstance().putEvent(ev);
                break;

            case SND_SEQ_EVENT_PORT_SUBSCRIBED: // ports connected
                if(true)
                    cout << "Info, alsa midi port connected" << endl;
                break;

            case SND_SEQ_EVENT_PORT_UNSUBSCRIBED: // ports disconnected
                if(true)
                    cout << "Info, alsa midi port disconnected" << endl;
                break;
            case SND_SEQ_EVENT_SENSING: // midi device still there
                break;
                
                
            case SND_SEQ_EVENT_SYSEX:   // system exclusive
                for (unsigned int x = 0; x < event->data.ext.len; x += 3) {
                    uint8_t buf[3];
                    int y = event->data.ext.len - x;
                    if (y >= 3) {
                        memcpy(buf, (uint8_t *)event->data.ext.ptr + x, 3);
                    } else {
                        memset(buf, 0, sizeof(buf));
                        memcpy(buf, (uint8_t *)event->data.ext.ptr + x, y);
                    }
                    midiProcess(buf[0], buf[1], buf[2], ev.time, ev.nanos);
                    
                }
                break;
            case SND_SEQ_EVENT_QFRAME :   // Midi Timecode
                midiProcess(0xf1, event->data.control.param, event->data.control.value, ev.time, ev.nanos);
                //~ ev.type  = M_TC;
                //~ ev.value = event->data.control.value;
                //~ InMgr::getInstance().putEvent(ev);
                break;
            case SND_SEQ_EVENT_SONGPOS:   // Song Position Pointer
                ev.type  = M_SPP;
                ev.value = event->data.control.value;
                InMgr::getInstance().putEvent(ev);
                break;
            case SND_SEQ_EVENT_CLOCK:   // Midi Clock
                ev.type  = M_CLOCK;
                InMgr::getInstance().putEvent(ev);
                break;
            case SND_SEQ_EVENT_START: /* M_START */
                midiProcess(0xFA, event->data.control.param, event->data.control.value, ev.time, ev.nanos);
                break;
            case SND_SEQ_EVENT_CONTINUE: /* M_CONTINUE */
                midiProcess(0xFB, event->data.control.param, event->data.control.value, ev.time, ev.nanos);
                break;
            case SND_SEQ_EVENT_STOP: /* M_STOP */
                midiProcess(0xFC, event->data.control.param, event->data.control.value, ev.time, ev.nanos);
                break;
            default:
                if(true)
                    cout << "Info, other non-handled midi event, type: "
                         << (int)event->type << endl;
                break;
        }
        snd_seq_free_event(event);
    }
    return NULL;
}

bool AlsaEngine::openMidi()
{
    if(getMidiEn())
        return true;

    //~ int alsaport, queue, client;
    midi.handle = NULL;

    if(snd_seq_open(&midi.handle, "default",
       SND_SEQ_OPEN_INPUT, SND_SEQ_NONBLOCK) != 0)
        return false;

    string clientname = "ZynAddSubFX";
    string postfix = Nio::getPostfix();
    if (!postfix.empty())
        clientname += "_" + postfix;
    if(Nio::pidInClientName)
        clientname += "_" + os_pid_as_padded_string();
    snd_seq_set_client_name(midi.handle, clientname.c_str());


	snd_seq_port_info_alloca(&midi.pinfo);
	snd_seq_port_info_set_name(midi.pinfo, "Input");
	snd_seq_port_info_set_type(midi.pinfo, SND_SEQ_PORT_TYPE_SYNTH);
	snd_seq_port_info_set_capability(midi.pinfo, SND_SEQ_PORT_CAP_WRITE | SND_SEQ_PORT_CAP_SUBS_WRITE);

	/* Enable timestamping for events sent by external subscribers. */
	snd_seq_port_info_set_timestamping(midi.pinfo, 1);
	snd_seq_port_info_set_timestamp_real(midi.pinfo, 1);
    int client, queue, port;
    if ((client = snd_seq_client_id(midi.handle))<0) {
		fprintf(stderr, "Cannot determine client number: %s\n", snd_strerror(client));
		return false;
	}
	printf("Client ID = %i\n", client);
	if ((queue = snd_seq_alloc_queue(midi.handle))<0) {
		fprintf(stderr, "Cannot allocate queue: %s\n", snd_strerror(queue));
		return false;
	}
    
    port = snd_seq_port_info_get_port(midi.pinfo);
    
    snd_seq_port_subscribe_t* sub;
	snd_seq_addr_t seq_addr;

	snd_seq_port_subscribe_alloca (&sub);
	seq_addr.client = client;
	seq_addr.port = port;
	snd_seq_port_subscribe_set_sender (sub, &seq_addr);
	//~ seq_addr.client = client;
	//~ seq_addr.port = port;
	//~ snd_seq_port_subscribe_set_dest (sub, &seq_addr);

	snd_seq_port_subscribe_set_time_update (sub, 1);
	snd_seq_port_subscribe_set_queue (sub, queue);
	snd_seq_port_subscribe_set_time_real (sub, 1);

	if ((snd_seq_create_port(midi.handle, midi.pinfo)) < 0) 
		return false;

    midi.exiting = false;
    pthread_attr_t attr;
    
	struct sched_param param; 
	pthread_attr_init(&attr);
	pthread_attr_setschedpolicy(&attr, SCHED_FIFO);
	pthread_attr_setinheritsched(&attr, PTHREAD_EXPLICIT_SCHED);
	param.sched_priority = 83;
	pthread_attr_setschedparam(&attr, &param); 

    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    pthread_create(&midi.pThread, &attr, _MidiThread, this);
    return true;
}

void AlsaEngine::stopMidi()
{
    if(!getMidiEn())
        return;

    snd_seq_t *handle = midi.handle;
    if((NULL != midi.handle) && midi.pThread) {
        midi.exiting = true;
        pthread_join(midi.pThread, 0);
    }
    midi.handle = NULL;
    if(handle)
        snd_seq_close(handle);
}

short *AlsaEngine::interleave(const Stereo<float *> &smps)
{
    /**\todo TODO fix repeated allocation*/
    short *shortInterleaved = audio.buffer;
    memset(shortInterleaved, 0, bufferSize * 2 * sizeof(short));
    int    idx = 0; //possible off by one error here
    double scaled;
    for(int frame = 0; frame < bufferSize; ++frame) { // with a nod to libsamplerate ...
        float l = smps.l[frame];
        float r = smps.r[frame];
        stereoCompressor(synth.samplerate, audio.peaks[0], l, r);

        scaled = l * (8.0f * 0x10000000);
        shortInterleaved[idx++] = (short int)(lrint(scaled) >> 16);
        scaled = r * (8.0f * 0x10000000);
        shortInterleaved[idx++] = (short int)(lrint(scaled) >> 16);
    }
    return shortInterleaved;
}

bool AlsaEngine::openAudio()
{
    if(getAudioEn())
        return true;

    int rc = 0;
    /* Open PCM device for playback. */
    audio.handle = NULL;

    const char *device = getenv("ALSA_DEVICE");
    if(device == 0)
        device = "hw:0";

    rc = snd_pcm_open(&audio.handle, device,
                      SND_PCM_STREAM_PLAYBACK, 0);
    if(rc < 0) {
        fprintf(stderr,
                "unable to open pcm device: %s\n",
                snd_strerror(rc));
        return false;
    }

    /* Allocate a hardware parameters object. */
    snd_pcm_hw_params_alloca(&audio.params);

    /* Fill it in with default values. */
    snd_pcm_hw_params_any(audio.handle, audio.params);

    /* Set the desired hardware parameters. */

    /* Interleaved mode */
    snd_pcm_hw_params_set_access(audio.handle, audio.params,
                                 SND_PCM_ACCESS_RW_INTERLEAVED);

    /* Signed 16-bit little-endian format */
    snd_pcm_hw_params_set_format(audio.handle, audio.params,
                                 SND_PCM_FORMAT_S16_LE);

    /* Two channels (stereo) */
    snd_pcm_hw_params_set_channels(audio.handle, audio.params, 2);

    audio.sampleRate = synth.samplerate;
    snd_pcm_hw_params_set_rate_near(audio.handle, audio.params,
                                    &audio.sampleRate, NULL);

    audio.frames = 512;
    snd_pcm_hw_params_set_period_size_near(audio.handle,
                                           audio.params, &audio.frames, NULL);

    audio.periods = 4;
    snd_pcm_hw_params_set_periods_near(audio.handle,
                                       audio.params, &audio.periods, NULL);

    /* Set buffer size (in frames). The resulting latency is given by */
    /* latency = periodsize * periods / (rate * bytes_per_frame)      */
    snd_pcm_uframes_t alsa_buffersize = synth.buffersize;
    rc = snd_pcm_hw_params_set_buffer_size_near(audio.handle,
                                               audio.params,
                                               &alsa_buffersize);

    /* At this place, ALSA's and zyn's buffer sizes may differ. */
    /* This should not be a problem.                            */
    if((int)alsa_buffersize != synth.buffersize)
     cerr << "ALSA buffer size: " << alsa_buffersize << endl;

    /* Write the parameters to the driver */
    rc = snd_pcm_hw_params(audio.handle, audio.params);
    if(rc < 0) {
        fprintf(stderr,
                "unable to set hw parameters: %s\n",
                snd_strerror(rc));
        return false;
    }
    
    /* Allocate a hardware parameters object. */
    snd_pcm_sw_params_alloca(&audio.swparams);
    
    rc = snd_pcm_sw_params_current(audio.handle, audio.swparams);
    if( rc < 0 )
    {
        fprintf (stderr, "cannot initialize software parameters structure (%s)\n", snd_strerror(rc) );
        return false;
    }
    
    /* Enable Timestamping */
    rc = snd_pcm_sw_params_set_tstamp_mode(audio.handle, audio.swparams, SND_PCM_TSTAMP_ENABLE);
    if (rc < 0) {
        fprintf(stderr, "Unable to set timestamp mode: %s\n", snd_strerror(rc));
        return false;
    }
    /* Set Timestamp type to Wallclock */
    rc = snd_pcm_sw_params_set_tstamp_type(audio.handle, audio.swparams, SND_PCM_TSTAMP_TYPE_GETTIMEOFDAY);
    if (rc < 0) {
        printf("Unable to set tstamp type : %s\n", snd_strerror(rc));
        return false;
    }
    
    /* Write the parameters to the driver */
    rc = snd_pcm_sw_params( audio.handle, audio.swparams );
        if(rc < 0) {
        fprintf(stderr,
                "unable to set sw parameters: %s\n",
                snd_strerror(rc));
        return false;
    }
    
    snd_pcm_hw_params_t *hwparams_c;
    snd_pcm_hw_params_alloca(&hwparams_c);
		/* get the current hwparams */
		rc = snd_pcm_hw_params_current(audio.handle, hwparams_c);
		if (rc < 0) {
			printf("Unable to determine current hwparams_c: %s\n", snd_strerror(rc));
			return false;
		}

		if (snd_pcm_hw_params_supports_audio_ts_type(hwparams_c, 0))
			printf("Capture supports audio compat timestamps\n");
		if (snd_pcm_hw_params_supports_audio_ts_type(hwparams_c, 1))
			printf("Capture supports audio default timestamps\n");
		if (snd_pcm_hw_params_supports_audio_ts_type(hwparams_c, 2))
			printf("Capture supports audio link timestamps\n");
		if (snd_pcm_hw_params_supports_audio_ts_type(hwparams_c, 3))
			printf("Capture supports audio link absolute timestamps\n");
		if (snd_pcm_hw_params_supports_audio_ts_type(hwparams_c, 4))
			printf("Capture supports audio link estimated timestamps\n");
		if (snd_pcm_hw_params_supports_audio_ts_type(hwparams_c, 5))
			printf("Capture supports audio link synchronized timestamps\n");
		if (snd_pcm_hw_params_supports_audio_wallclock_ts(hwparams_c))
			printf("Capture supports wallclock timestamps\n");
            
    snd_pcm_status_alloca(&audio.status);
        audio_tstamp_config.type_requested = 4;
        audio_tstamp_config.report_delay = 1;
        
    snd_pcm_status_set_audio_htstamp_config(audio.status, &audio_tstamp_config);
        if ((rc = snd_pcm_status(audio.handle, audio.status)) < 0) {
            printf("Stream status error: %s\n", snd_strerror(rc));
            exit(0);
        }
    snd_pcm_status_get_trigger_htstamp(audio.status, &htstamptrigger);
    tstamptrigger = htstamptrigger.tv_sec * 1000000ULL + htstamptrigger.tv_nsec/1000ULL;
    
    
    printf("htstamptrigger.tv_sec: %lu tv_nsec: %lu \n", htstamptrigger.tv_sec, htstamptrigger.tv_nsec);
            

    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
       
    pthread_create(&audio.pThread, &attr, _AudioThread, this);
    
    framecount = 0;
    return true;
}

void AlsaEngine::stopAudio()
{
    if(!getAudioEn())
        return;

    snd_pcm_t *handle = audio.handle;
    audio.handle = NULL;
    pthread_join(audio.pThread, NULL);
    snd_pcm_drain(handle);
    if(snd_pcm_close(handle))
        cout << "Error: in snd_pcm_close " << __LINE__ << ' ' << __FILE__
             << endl;
}

void *AlsaEngine::processAudio()
{
    while(audio.handle) {
        int rc;
        snd_pcm_t *handle = audio.handle;
        
        
        
        
        
        snd_pcm_status_get_htstamp(audio.status, &htstampaudio);
        
        // TBD: tstamp needs to be known in Master
        tstampaudio = htstampaudio.tv_sec * 1000000ULL + htstampaudio.tv_nsec/1000ULL;
        
        
        audio.buffer = interleave(getNext());
        
        rc = snd_pcm_writei(handle, audio.buffer, synth.buffersize);
        if(rc == -EPIPE) {
            /* EPIPE means underrun */
            cerr << "underrun occurred" << endl;
            snd_pcm_prepare(handle);
        }
        else
        if(rc < 0) {
            cerr << "AlsaEngine: Recovering connection..." << endl;
            rc = snd_pcm_recover(handle, rc, 0);
            if(rc < 0)
             throw "Could not recover ALSA connection";
        }
        

        framecount+=synth.buffersize;
    }
    return NULL;
}

}
