/*
  ZynAddSubFX - a software synthesizer

  SEQParams.h - Parameters for SEQ
  Copyright (C) 2021 Michael Kirchner
  Author: Michael Kirchner

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.
*/

#ifndef SEQ_PARAMS_H
#define SEQ_PARAMS_H

#include <Misc/Time.h>
#include <rtosc/ports.h>
#include "Presets.h"

#define MAX_CUTOFF 20.0f

namespace zyn {

class XMLwrapper;

class SEQParams:public Presets
{
    public:
        SEQParams(const AbsTime* time_ = nullptr);
        SEQParams(consumer_location_t loc,
                  const AbsTime* time_ = nullptr);
        SEQParams(float freq_,
                  float cutoff_,
                  float intensity_,
                  int steps_,
                  float delay_,
                  bool continous,
                  float speedratio,
                  consumer_location_t loc,
                  const AbsTime* time_ = nullptr);
        ~SEQParams() override;

        void add2XML(XMLwrapper& xml) override;
        void defaults();
        /**Loads the SEQ from the xml*/
        void getfromXML(XMLwrapper& xml);
        void paste(SEQParams &);


        float freq;
        float cutoff;
        float intensity;
        float delay;
        unsigned char continous; 
        unsigned char ratiofixed; 
        int steps; 
        float         speedratio;
        float*        sequence;

        //! what kind is the SEQ (0 - frequency, 1 - amplitude, 2 - filter)
        consumer_location_type_t fel;
        int loc; //!< consumer location

        const AbsTime *time;
        int64_t last_update_timestamp;

        static const rtosc::Ports &ports;
    private:
        //! common functionality of ctors
        void setup();

        /* Default parameters */
        float         Dfreq;
        float         Dcutoff;
        float         Dintensity;
        float         Ddelay;
        unsigned char Dcontinous;
        float         Dspeedratio;
        int           Dsteps;
        float*        Dsequence;
};

}

#endif
