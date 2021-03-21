/*
  ZynAddSubFX - a software synthesizer

  SEQParams.cpp - Parameters for LFO
  Copyright (C) 2021 Michael Kirchner
  Author: Michael Kirchner

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.
*/

#include <cmath>
#include <cstdio>
#include "../globals.h"
#include "../Misc/Util.h"
#include "../Misc/XMLwrapper.h"
#include "../Misc/Time.h"
#include "SEQParams.h"
#include <rtosc/port-sugar.h>
#include <rtosc/ports.h>
using namespace rtosc;

namespace zyn {


#define rObject SEQParams
#undef rChangeCb
#define rChangeCb if (obj->time) { obj->last_update_timestamp = obj->time->time(); }
#define rBegin [](const char *msg, rtosc::RtData &d) {
#define rEnd }
  
    
#define rCrossOption(name, crossname, setcode, ...) \
  {STRINGIFY(name) "::i:c:S",rProp(parameter) rProp(enumerated) DOC(__VA_ARGS__), NULL, rCrossOptionCb(name, crossname, setcode)}
  
#define rCrossOptionCb(name, cross, setcode) rBOIL_BEGIN \
    { \
        if(!strcmp("i", args) || !strcmp("I", args)) { \
            int var = rtosc_argument(msg, 0).i; \
            rLIMIT(var, atoi) \
            setcode; \
            /* cross broadcast to speedratio port */ \
            char part_loc[128]; \
            strncpy(part_loc, data.loc, sizeof(part_loc)); \
            part_loc[sizeof(part_loc) - 1] = '\0'; \
            char *end = strrchr(part_loc, '/'); \
            if(end) { \
                strcpy(&end[1], STRINGIFY(cross)); \
                data.broadcast(part_loc, "f", obj->cross); \
            } \
            rChangeCb; \
        } \
    } \
    rBOIL_END

static const rtosc::Ports _ports = {
    rSelf(SEQParams),
    rPaste,
    rOption(loc, rProp(internal),
            rOptions(ad_global_amp, ad_global_freq, ad_global_filter,
                     ad_voice_amp, ad_voice_freq, ad_voice_filter, unspecified),
            "location of the filter"),
    rParamF(freq, rShort("freq"), rDefault(2.0f), rUnit(Hz), rLinear(0.1f,10.0f),
            "step frequency"),
    rParamF(cutoff, rShort("lp"), rDefault(0.0f), rUnit(Hz), rLinear(0.0f,MAX_CUTOFF),
            "cutoff of lp-filter for output\n 0.0=off"),
    rParamF(delay, rShort("delay"), rUnit(S),
              rLinear(0.0, 4.0), rDefault(0.0),
              "Delay before SEQ start\n0..4 second delay"),
    rToggle(continous, rShort("c"), rDefault(false),
            "Enable for global operation"),
    rCrossOption(ratiofixed, speedratio, obj->speedratio = speedratios[var], 
        rShort("rat"), 
        rOptions(off, 8/1, 7/1, 6/1, 5/1, 4/1, 3/1, 2/1, 7/4, 3/2, 5/4, 1/1, 7/8, 3/4, 1/2, 1/3, 1/4, 1/5, 1/6, 1/7, 1/8), 
        rLinear(0,20), rDefault(off), "select fixed ratio for BPM sync"),
    rParamF(speedratio, rShort("r"), rLinear(0.0f,8.0f), rDefault(0.0f),
            "ratio for BPM sync"),
    rParamZyn(steps, rShort("steps"), rLinear(1,128), rDefault(8),
            "number of steps"),
    rArrayF(sequence, NUM_SEQ_STEPS, rLinear(-1.0f,1.0f), rDefaultDepends(loc),
              rPreset(ad_global_freq, [0.0 0.0 0.0 ...]),   
              rPreset(ad_global_filter, [0.0 0.0 0.0 ...]),
              rPreset(ad_global_amp, [1.0 1.0 1.0 ...]),
              rPreset(ad_voice_freq, [0.0 0.0 0.0 ...]),   
              rPreset(ad_voice_filter, [0.0 0.0 0.0 ...]),
              rPreset(ad_voice_amp, [1.0 1.0 1.0 ...]),
           "sequence values"),
};
#undef rBegin
#undef rEnd
#undef rChangeCb

const rtosc::Ports &SEQParams::ports = _ports;

void SEQParams::setup()
{
    switch(loc) {
        case loc_unspecified:
            fel = consumer_location_type_t::unspecified;
            break;
        case ad_global_freq:
        case ad_voice_freq:
            fel = consumer_location_type_t::freq;
            setpresettype("Pseqfrequency");
            break;
        case ad_global_amp:
        case ad_voice_amp:
            fel = consumer_location_type_t::amp;
            setpresettype("Pseqamplitude");
            break;
        case ad_global_filter:
        case ad_voice_filter:
            fel = consumer_location_type_t::filter;
            setpresettype("Pseqfilter");
            break;
        default:
            throw std::logic_error("Invalid envelope consumer location");
    }

    defaults();
}

// TODO: reuse
SEQParams::SEQParams(const AbsTime *time_) :
    SEQParams(2.0f, 20.0f, 8, 0.0f, false, 0.0f, loc_unspecified, time_)
{
}

SEQParams::SEQParams(float freq_,
                     float cutoff_,
                     unsigned char steps_,
                     float delay_,
                     bool continous_,
                     float speedratio_,
                     consumer_location_t loc,
                     const AbsTime *time_) : loc(loc),
                                             time(time_),
                                             last_update_timestamp(0) {
    Dfreq       = freq_;
    Dcutoff     = cutoff_;
    Ddelay      = delay_;
    Dcontinous  = continous_;
    Dspeedratio = speedratio_;
    Dsteps      = steps_;

    setup();
}

SEQParams::SEQParams(consumer_location_t loc,
                     const AbsTime *time_) : loc(loc),
                                             time(time_),
                                             last_update_timestamp(0) {

    auto init =
        [&](float freq_, char cutoff_, unsigned char steps_, float val_)
    {
        Dfreq       = freq_;
        Dcutoff     = cutoff_;
        Dsteps      = steps_;
        Ddelay      = 0.0f;
        Dcontinous  = false;
        Dspeedratio = 1.0f;
        if (val_ == 0.0f)
            memset(Dsequence, 0, NUM_SEQ_STEPS*sizeof(*sequence));
        else
            for(int i = 0; i<NUM_SEQ_STEPS; ++i)
                Dsequence[i] = val_;

    };

    switch(loc)
    {                    // (float freq_, char cutoff_, unsigned char steps_, float val_)
        case ad_global_amp:    init(2.0f, MAX_CUTOFF, 8, 1.0f); break;
        case ad_global_freq:   init(2.0f, MAX_CUTOFF, 8, 0.0f); break;
        case ad_global_filter: init(2.0f, MAX_CUTOFF, 8, 0.0f); break;
        case ad_voice_amp:     init(2.0f, MAX_CUTOFF, 8, 1.0f); break;
        case ad_voice_freq:    init(2.0f, MAX_CUTOFF, 8, 0.0f); break;
        case ad_voice_filter:  init(2.0f, MAX_CUTOFF, 8, 0.0f); break;
        default: throw std::logic_error("Invalid SEQ consumer location");
    }
    setup();
}

SEQParams::~SEQParams()
{}

void SEQParams::defaults()
{
    freq       = Dfreq;
    cutoff     = Dcutoff;
    delay      = Ddelay;
    continous  = Dcontinous;
    speedratio = Dspeedratio;
    memcpy(sequence, Dsequence, NUM_SEQ_STEPS*sizeof(*sequence));
}


void SEQParams::add2XML(XMLwrapper& xml)
{
    xml.addparreal("freq", freq);
    xml.addparreal("cutoff", cutoff);
    xml.addpar("steps", steps);
    xml.addparreal("delay", delay);
    xml.addparbool("continous", continous);
    xml.addparreal("speedratio", speedratio);
    xml.beginbranch("SEQUENCE");
    for(int n = 0; n < NUM_SEQ_STEPS; ++n) {
        if((sequence[n] == 0.0f))
            continue;
        xml.beginbranch("STEP", n + 1);
        xml.addparreal("value", sequence[n]);
        xml.endbranch();
    }
    xml.endbranch(); // SEQUENCE
}

void SEQParams::getfromXML(XMLwrapper& xml)
{

    freq       = xml.getparreal("freq", freq);
    cutoff     = xml.getparreal("cutoff", cutoff);
    steps      = xml.getpar("steps", steps,1,127);
    delay      = xml.getparreal("delay", delay);
    continous  = xml.getparbool("continous", continous);
    speedratio = xml.getparreal("speedratio", speedratio);
    
    if(xml.enterbranch("SEQUENCE")) {
        for(int n = 0; n < NUM_SEQ_STEPS; ++n) {
            sequence[n] = 0.0f;
            if(xml.enterbranch("STEP", n + 1) == 0)
                continue;
            sequence[n] = xml.getparreal("value", 0.0f);
            xml.exitbranch();
        }
        xml.exitbranch(); // SEQUENCE
    }

}

#define COPY(y) this->y=x.y
void SEQParams::paste(SEQParams &x)
{
    COPY(freq);
    COPY(cutoff);
    COPY(steps);
    COPY(delay);
    COPY(continous);
    COPY(speedratio);
    memcpy(sequence, x.sequence, NUM_SEQ_STEPS*sizeof(*sequence));

    if ( time ) {
        last_update_timestamp = time->time();
    }
}
#undef COPY

}
