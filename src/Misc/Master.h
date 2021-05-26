/*
  ZynAddSubFX - a software synthesizer

  Master.h - It sends Midi Messages to Parts, receives samples from parts,
             process them with system/insertion effects and mix them
  Copyright (C) 2002-2005 Nasca Octavian Paul
  Author: Nasca Octavian Paul

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.
*/

#ifndef MASTER_H
#define MASTER_H
#include "../globals.h"
#include "Microtonal.h"
#include <atomic>
#include <rtosc/automations.h>
#include <rtosc/miditable.h>
#include <rtosc/savefile.h>

#include "Time.h"
#include "Bank.h"
#include "Recorder.h"

#include "../Params/Controller.h"
#include "../Synth/WatchPoint.h"
#include "../DSP/Value_Smoothing_Filter.h"

namespace zyn {

class Allocator;

struct vuData {
    vuData(void);
    float outpeakl, outpeakr, maxoutpeakl, maxoutpeakr,
          rmspeakl, rmspeakr;
    int clipped;
};


/** It sends Midi Messages to Parts, receives samples from parts,
 *  process them with system/insertion effects and mix them */
class Master
{
    public:
        Master(const Master& other) = delete;
        Master(Master&& other) = delete;

        /** Constructor TODO make private*/
        Master(const SYNTH_T &synth, class Config *config);
        /** Destructor*/
        ~Master();

        char last_xmz[XMZ_PATH_MAX];

        //applyOscEvent overlays
        bool applyOscEvent(const char *event, float *outl, float *outr,
                           bool offline, bool nio = true, int msg_id = -1);
        bool applyOscEvent(const char *event, bool nio = true, int msg_id = -1);

        /**Saves all settings to a XML file
         * @return 0 for ok or <0 if there is an error*/
        int saveXML(const char *filename);

        /**This adds the parameters to the XML data*/
        void add2XML(XMLwrapper& xml);

        static void saveAutomation(XMLwrapper &xml, const rtosc::AutomationMgr &midi);
        static void loadAutomation(XMLwrapper &xml,       rtosc::AutomationMgr &midi);

        void defaults();

        /**loads all settings from a XML file
         * @return 0 for ok or -1 if there is an error*/
        int loadXML(const char *filename);

        /**Append all settings to an OSC savefile (as specified by RT OSC)*/
        std::string saveOSC(std::string savefile);
        /**loads all settings from an OSC file (as specified by RT OSC)
         * @param dispatcher Message dispatcher and modifier
         * @return 0 for ok or <0 if there is an error*/
        int loadOSC(const char *filename,
                    rtosc::savefile_dispatcher_t* dispatcher);

        /**Regenerate PADsynth and other non-RT parameters
         * It is NOT SAFE to call this from a RT context*/
        void applyparameters(void) NONREALTIME;

        //This must be called prior-to/at-the-time-of RT insertion
        void initialize_rt(void) REALTIME;

        void getfromXML(XMLwrapper& xml);

        /**get all data to a newly allocated array (used for plugin)
         * @return the datasize*/
        int getalldata(char **data) NONREALTIME;
        /**put all data from the *data array to zynaddsubfx parameters (used for plugin)*/
        void putalldata(const char *data);

        //Midi IN
        void noteOn(char chan, note_t note, char velocity) {
            noteOn(chan, note, velocity, note / 12.0f);
        };
        void noteOn(char chan, note_t note, char velocity, float note_log2_freq);
        void noteOff(char chan, note_t note);
        void polyphonicAftertouch(char chan, note_t note, char velocity);
        void setController(char chan, int type, int par);
        void setController(char chan, int type, note_t note, float value);
        //void NRPN...


        void ShutUp();
        int shutup;

        void vuUpdate(const float *outl, const float *outr);

        //Process a set of OSC events in the bToU buffer
        //This may be called by MiddleWare if we are offline
        //(in this case, the param offline is true)
        bool runOSC(float *outl, float *outr, bool offline=false,
                    Master* master_from_mw = nullptr);

        /**Audio Output*/
        bool AudioOut(float *outl, float *outr) REALTIME;
        /**Audio Output (for callback mode).
         * This allows the program to be controlled by an external program*/
        void GetAudioOutSamples(size_t nsamples,
                                unsigned samplerate,
                                float *outl,
                                float *outr) REALTIME;


        void partonoff(int npart, int what);

        //Set callback to run when master changes
        void setMasterChangedCallback(void(*cb)(void*,Master*),void *ptr);
        //Copy callback to other master
        void copyMasterCbTo(Master* dest);
        bool hasMasterCb() const;

        /**parts \todo see if this can be made to be dynamic*/
        class Part * part[NUM_MIDI_PARTS];

        //parameters
        unsigned char Pkeyshift;
        unsigned char Psysefxvol[NUM_SYS_EFX][NUM_MIDI_PARTS];
        unsigned char Psysefxsend[NUM_SYS_EFX][NUM_SYS_EFX];

        //parameters control
        static float volume127ToFloat(unsigned char volume_);
        void setPkeyshift(char Pkeyshift_);
        void setPsysefxvol(int Ppart, int Pefx, char Pvol);
        void setPsysefxsend(int Pefxfrom, int Pefxto, char Pvol);

        //effects
        class EffectMgr * sysefx[NUM_SYS_EFX]; //system
        class EffectMgr * insefx[NUM_INS_EFX]; //insertion
//      void swapcopyeffects(int what,int type,int neff1,int neff2);

        //HDD recorder
        Recorder HDDRecorder;

        //part that's apply the insertion effect; -1 to disable
        short int Pinsparts[NUM_INS_EFX];


        //peaks for VU-meter
        void vuresetpeaks();

        //peaks for part VU-meters
        float vuoutpeakpartl[NUM_MIDI_PARTS];
        float vuoutpeakpartr[NUM_MIDI_PARTS];
        unsigned char fakepeakpart[NUM_MIDI_PARTS]; //this is used to compute the "peak" when the part is disabled

        AbsTime time;
        Controller ctl;
        bool       swaplr; //if L and R are swapped

        //other objects
        Microtonal microtonal;

        //Strictly Non-RT instrument bank object
        Bank bank;

        class FFTwrapper * fft;

        static const rtosc::Ports &ports;
        float  Volume;

        //Statistics on output levels
        vuData vu;

        //Display info on midi notes
        bool activeNotes[128];

        //Other watchers
        WatchManager watcher;

        //Midi Learn
        rtosc::AutomationMgr automate;
        rtosc::MidiMapperRT midi;

        bool   frozenState;//read-only parameters for threadsafe actions
        Allocator *memory;
        rtosc::ThreadLink *bToU;
        rtosc::ThreadLink *uToB;
        bool pendingMemory;
        const SYNTH_T &synth;
        const int& gzip_compression; //!< value from config

        //Heartbeat for identifying plugin offline modes
        //in units of 10 ms (done s.t. overflow is in 497 days)
        uint32_t last_beat = 0;
        uint32_t last_ack = 0;

        //Buffer to contain the OSC path to the last GUI element
        //on which a drag and drop operation ended
        constexpr static std::size_t dnd_buffer_size = 1024;
        char dnd_buffer[dnd_buffer_size] = {0};

        //Return XML data as string. Must be freed.
        char* getXMLData();
        //Load OSC from OSC savefile
        //Returns 0 if OK, <0 in case of failure
        int loadOSCFromStr(const char *file_content,
                           rtosc::savefile_dispatcher_t* dispatcher);

    private:
        std::atomic<bool> run_osc_in_use = { false };

        float  sysefxvol[NUM_SYS_EFX][NUM_MIDI_PARTS];
        float  sysefxsend[NUM_SYS_EFX][NUM_SYS_EFX];
        int    keyshift;

        //information relevant to generating plugin audio samples
        float *bufl;
        float *bufr;
        off_t  off;
        size_t smps;

        //Callback When Master changes
        void(*mastercb)(void*,Master*);
        void* mastercb_ptr;

        //! apply an OSC event with a DataObj parameter
        //! @note This may be called by MiddleWare if we are offline
        //!   (in this case, the param offline is true)
        //! @return false iff master has been changed
        bool applyOscEvent(const char *event, float *outl, float *outr,
                           bool offline, bool nio,
                           class DataObj& d, int msg_id = -1,
                           Master* master_from_mw = nullptr);

        Value_Smoothing_Filter smoothing;

        Value_Smoothing_Filter smoothing_part_l[NUM_MIDI_PARTS];
        Value_Smoothing_Filter smoothing_part_r[NUM_MIDI_PARTS];
};

class master_dispatcher_t : public rtosc::savefile_dispatcher_t
{
    virtual void vUpdateMaster(Master* m) = 0;
public:
    void updateMaster(Master* m) { vUpdateMaster(m); }
};

}

#endif
