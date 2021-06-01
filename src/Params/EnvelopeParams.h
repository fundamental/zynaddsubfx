/*
  ZynAddSubFX - a software synthesizer

  EnvelopeParams.h - Parameters for Envelope
  Copyright (C) 2002-2005 Nasca Octavian Paul
  Author: Nasca Octavian Paul

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.
*/

#ifndef ENVELOPE_PARAMS_H
#define ENVELOPE_PARAMS_H

#include "../globals.h"
#include "../Misc/XMLwrapper.h"
#include "Presets.h"

namespace zyn {

   enum envmode_enum {
        ADSR_lin = 1,
        ADSR_dB = 2,
        ASR_freqlfo = 3,
        ADSR_filter = 4,
        ASR_bw = 5
    };

class EnvelopeParams:public Presets
{
    public:
        EnvelopeParams(unsigned char Penvstretch_=64,
                       unsigned char Pforcedrelease_=0,
                       const AbsTime *time_ = nullptr);
        ~EnvelopeParams() override;
        void paste(const EnvelopeParams &ep);

        void init(consumer_location_t loc);
        void converttofree();
        void updatefree();
        void updatenonfree();

        void add2XML(XMLwrapper& xml) override;
        void defaults();
        void getfromXML(XMLwrapper& xml);

        float getdt(char i) const;


        int loc; //!< consumer location

        /* MIDI Parameters */
        unsigned char Pfreemode; //1 for free mode, 0 otherwise
        unsigned char Penvpoints;
        unsigned char Penvsustain; //127 for disabled
        float         envdt[MAX_ENVELOPE_POINTS];
        float         envcp[MAX_ENVELOPE_POINTS-1][4];  //control points {rx, ry, lx, ly}
        unsigned char Penvval[MAX_ENVELOPE_POINTS];
        unsigned char Penvstretch; //64=normal stretch (piano-like), 0=no stretch
        unsigned char Pforcedrelease; //0 - OFF, 1 - ON
        unsigned char Plinearenvelope; //1 for linear AMP ENV, 0 otherwise
        unsigned char Prepeating; //0 - OFF, 1 - ON

        float A_dt, D_dt, R_dt;
        unsigned char PA_val, PD_val, PS_val, PR_val;




        envmode_enum Envmode; // 1 for ADSR parameters (linear amplitude)
                     // 2 for ADSR_dB parameters (dB amplitude)
                     // 3 for ASR parameters (frequency LFO)
                     // 4 for ADSR_filter parameters (filter parameters)
                     // 5 for ASR_bw parameters (bandwidth parameters)



        const AbsTime *time;
        int64_t last_update_timestamp;

        static const rtosc::Ports &ports;

        static float env_rap2dB(float rap);
        static float env_dB2rap(float db);

    private:

        void ADSRinit(float A_dt, float D_dt, char S_val, float R_dt);
        void ADSRinit_dB(float A_dt, float D_dt, char S_val, float R_dt);
        void ASRinit(char A_val, float A_dt, char R_val, float R_dt);
        void ADSRinit_filter(char A_val,
                             float A_dt,
                             char D_val,
                             float D_dt,
                             float R_dt,
                             char R_val);
        void ASRinit_bw(char A_val, float A_dt, char R_val, float R_dt);

        void store2defaults();

        /* Default parameters */
        unsigned char Denvstretch;
        unsigned char Dforcedrelease;
        unsigned char Dlinearenvelope;
        unsigned char Drepeating;
        float DA_dt, DD_dt, DR_dt;
        unsigned char DA_val, DD_val, DS_val, DR_val;
};

}

#endif
