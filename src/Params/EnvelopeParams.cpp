/*
  ZynAddSubFX - a software synthesizer

  EnvelopeParams.cpp - Parameters for Envelope
  Copyright (C) 2002-2005 Nasca Octavian Paul
  Author: Nasca Octavian Paul

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.
*/

#include <cmath>
#include <cstdlib>
#include <cassert>
#include <rtosc/ports.h>
#include <rtosc/port-sugar.h>

#include "zyn-version.h"
#include "EnvelopeParams.h"
#include "../Misc/Util.h"
#include "../Misc/Time.h"

using namespace rtosc;

namespace zyn {

#define rObject EnvelopeParams
#define rBegin [](const char *msg, RtData &d) { \
    (void)* msg; EnvelopeParams *env = (rObject*) d.obj
#define rEnd }

#define dT127(var) limit( (int)roundf(log2f(var*100.0f + 1.0f) * 127.0f/12.0f), 0, 127 )
#define dTREAL(var) (powf(2.0f, var / 127.0f * 12.0f) - 1.0f) / 100.0f

#define rParamDT(name, ...) \
  {"P" STRINGIFY(name) "::i", rProp(alias) rProp(parameter) DOC(__VA_ARGS__), NULL, rParamDTCb(name)}

#define rParamDTCb(name) rBOIL_BEGIN \
        if(!strcmp("", args)) {\
            data.reply(loc, "i", dT127(obj->name)); \
        } else { \
            unsigned char var = rtosc_argument(msg, 0).i; \
            rLIMIT(var, atoi) \
            rCAPPLY(obj->name, "f", obj->name = dTREAL(var)) \
            data.broadcast(loc, "i", dT127(obj->name));\
            rChangeCb \
        } rBOIL_END

#define rParamsDT(name, length, ...) \
rArrayDT(name, length, __VA_ARGS__), \
{"P" STRINGIFY(name) ":", rProp(alias) rDoc("get all data from aliased array"), NULL, rParamsDTCb(name, length)}

#define rParamsDTCb(name, length) rBOIL_BEGIN \
    char varS[length]; \
    for (int i = 0; i < length; i++) {varS[i] = dT127(obj->name[i]);} \
    data.reply(loc, "b", length, varS); \
    rChangeCb \
    rBOIL_END

#define rArrayDT(name, length, ...) \
{"P" STRINGIFY(name) "#" STRINGIFY(length) "::i", rProp(parameter) DOC(__VA_ARGS__), NULL, rArrayDTICb(name)}

#define rArrayDTICb(name) rBOILS_BEGIN \
        if(!strcmp("", args)) {\
            data.reply(loc, "i", dT127(obj->name[idx])); \
        } else { \
            char varI = rtosc_argument(msg, 0).i; \
            float var = dTREAL(varI); \
            rLIMIT(var, atoi) \
            rAPPLY(name[idx], f) \
            data.broadcast(loc, "i", dT127(obj->name[idx]));\
            rChangeCb \
        } rBOILS_END

static const rtosc::Ports localPorts = {
    rSelf(EnvelopeParams),
    rPaste,
#undef  rChangeCb
#define rChangeCb if(!obj->Pfreemode) obj->converttofree(); if (obj->time) { \
            obj->last_update_timestamp = obj->time->time(); }
    rOption(loc, rProp(internal),
            rOptions(ad_global_amp, ad_global_freq, ad_global_filter,
                     ad_voice_amp, ad_voice_freq, ad_voice_filter,
                     ad_voice_fm_amp, ad_voice_fm_freq,
                     sub_freq, sub_bandwidth),
            "location of the envelope"),
    rToggle(Pfreemode, rDefault(false), "Complex Envelope Definitions"),
#undef  rChangeCb
#define rChangeCb if(!obj->Pfreemode) obj->converttofree(); \
                  else obj->updatenonfree(); \
                  if(obj->time) { obj->last_update_timestamp = obj->time->time(); }
    rParamZyn(Penvpoints, rProp(internal), rDefaultDepends(loc),
            rPresetAtMulti(3, ad_global_freq, ad_voice_freq,
                              ad_voice_fm_freq,
                              sub_freq, sub_bandwidth),
            rDefault(4),
            "Number of points in complex definition"),
    rParamZyn(Penvsustain, rDefaultDepends(loc),
            rPresetAtMulti(1, ad_global_freq, ad_voice_freq,
                              ad_voice_fm_freq,
                              sub_freq, sub_bandwidth),
            rDefault(2),
            "Location of the sustain point"),
    rParamsDT(envdt,  MAX_ENVELOPE_POINTS, "Envelope Delay Times"),
    rParams(Penvval, MAX_ENVELOPE_POINTS, "Envelope Values"),
    rParamZyn(Penvstretch,  rShort("stretch"), rDefaultDepends(loc),
            rPresetAtMulti(0, ad_global_freq, ad_global_filter,
                              ad_voice_freq, ad_voice_filter, ad_voice_fm_freq),
            rDefault(64),
            "Stretch with respect to frequency"),
    rToggle(Pforcedrelease, rShort("frcr"), rDefaultDepends(loc),
            rPresetAtMulti(true, ad_global_amp, ad_global_filter, ad_voice_amp,
                                 ad_voice_fm_amp),
            rDefault(false),
            "Force Envelope to fully evaluate"),
    rToggle(Plinearenvelope, rShort("lin/log"), rDefault(false),
            "Linear or Logarithmic Envelopes"),
    rToggle(Prepeating, rShort("repeat"), rDefault(false),
            "Repeat the Envelope"),
#undef  rChangeCb
#define rChangeCb if(!obj->Pfreemode) obj->converttofree(); \
                  else obj->updatefree(); \
                  if(obj->time) { obj->last_update_timestamp = obj->time->time(); }\
                  char part_loc[128]; \
                  strncpy(part_loc, loc, sizeof(part_loc)); \
                  part_loc[sizeof(part_loc) - 1] = '\0'; \
                  char *end = strrchr(part_loc, '/'); \
                  if(end) { \
                      end[1] = '\0'; \
                      data.broadcast("/damage", "s", part_loc);}
    rParamDT(A_dt ,  rShort("a.dt"), rLinear(0,127), "Attack Time"),
    rParamF(A_dt,  rShort("a.dt"), rLog(0.0f,41.0f), rDefaultDepends(loc),
              rPreset(ad_global_freq, 0.254f),   rPreset(ad_global_filter, 0.127f),
              rPreset(ad_voice_freq, 0.127f),    rPreset(ad_voice_filter, 0.970f),
              rPreset(ad_voice_fm_freq, 3.620f), rPreset(ad_voice_fm_amp, 1.876f),
              rPreset(sub_freq, 0.254f),         rPreset(sub_bandwidth, 0.970f),
              rDefault(0.0f),
              "Attack Time"),
    rParamZyn(PA_val, rShort("a.val"), rDefaultDepends(loc),
              rPreset(ad_voice_freq, 30),    rPreset(ad_voice_filter, 90),
              rPreset(ad_voice_fm_freq, 20),
              rPreset(sub_freq, 30),         rPreset(sub_bandwidth, 100),
              rDefault(64),
              "Attack Value"),
    rParamDT(D_dt, rShort("d.dt"), rLinear(0,127), "Decay Time"),
    rParamF(D_dt,  rShort("d.dt"), rLog(0.0f,41.0f),  rDefaultDepends(loc),
              rPreset(ad_global_amp, 0.127f),    rPreset(ad_global_filter, 0.970f),
              rPreset(ad_voice_amp, 6.978f),    rPreset(ad_voice_filter, 0.970f),
              rPreset(ad_voice_fm_amp, 3.620f),
              rDefault(0.009f),
              "Decay Time"),
    rParamZyn(PD_val, rShort("d.val"), rDefaultDepends(loc),
              rDefault(64), rPreset(ad_voice_filter, 40),
              "Decay Value"),
    rParamZyn(PS_val, rShort("s.val"), rDefaultDepends(loc),
              rDefault(64),
              rPresetAtMulti(127, ad_global_amp, ad_voice_amp, ad_voice_fm_amp),
              "Sustain Value"),
    rParamDT(R_dt, rShort("r.dt"), rLinear(0,127), "Release Time"),
    rParamF(R_dt,  rShort("r.dt"), rLog(0.009f,41.0f),  rDefaultDepends(loc),
              rPreset(ad_global_amp, 0.041f),
              rPreset(ad_voice_amp, 6.978f),    rPreset(ad_voice_filter, 0.009f),
              rPreset(ad_voice_fm_freq, 1.876f), rPreset(ad_voice_fm_amp, 6.978f),
              rDefault(0.499f),
              "Release Time"),
    rParamZyn(PR_val, rShort("r.val"), rDefaultDepends(loc),
              rPresetAtMulti(40, ad_voice_filter, ad_voice_fm_freq),
              rDefault(64),
              "Release Value"),
    {"Envmode:", rDoc("Envelope variant type"), NULL,
        rBegin;
        d.reply(d.loc, "i", env->Envmode);
        rEnd},

    {"envdt", rDoc("Envelope Delay Times (ms)"), NULL,
        rBegin;
        const int N = MAX_ENVELOPE_POINTS;
        const int M = rtosc_narguments(msg);
        if(M == 0) {
            rtosc_arg_t args[N];
            char arg_types[N+1] = {0};
            for(int i=0; i<N; ++i) {
                args[i].f    = env->getdt(i)*1000; //answer milliseconds to old gui
                arg_types[i] = 'f';
            }
            d.replyArray(d.loc, arg_types, args);
        } else {
            for(int i=0; i<N && i<M; ++i) {
                env->envdt[i] = (rtosc_argument(msg, i).f)/1000; //store as seconds in member variable
            }
            env->updatenonfree();
              char part_loc[128];
              strncpy(part_loc, d.loc, sizeof(part_loc));
              part_loc[sizeof(part_loc) - 1] = '\0';
              char *end = strrchr(part_loc, '/');
              if(end) {
                  end[1] = '\0';
                  d.broadcast("/damage", "s", part_loc);
                }
        }
        rEnd},
    {"dt", rDoc("Envelope Delay Times (sec)"), NULL,
        rBegin;

        const int N = MAX_ENVELOPE_POINTS;
        const int M = rtosc_narguments(msg);
        if(M == 0) {
            rtosc_arg_t args[N];
            char arg_types[N+1] = {};
            for(int i=0; i<N; ++i) {
                args[i].f    = env->getdt(i);
                arg_types[i] = 'f';
            }
            d.replyArray(d.loc, arg_types, args);
        } else {
            for(int i=0; i<N && i<M; ++i) {
                env->envdt[i] = (rtosc_argument(msg, i).f);
            }
            env->updatenonfree();
              char part_loc[128];
              strncpy(part_loc, d.loc, sizeof(part_loc));
              part_loc[sizeof(part_loc) - 1] = '\0';
              char *end = strrchr(part_loc, '/');
              if(end) {
                  end[1] = '\0';
                  d.broadcast("/damage", "s", part_loc);
              }
        }
        rEnd},
    {"envcp", rShort("bezier") rDoc("Envelope Control Points"), NULL,
        rBegin;
        const int N = MAX_ENVELOPE_POINTS-1;
        const int M = rtosc_narguments(msg);
        if(M == 0) {
            rtosc_arg_t args[N];
            char arg_types[4*(N+1)] = {};
            for(int i=0; i<N; ++i) {
                for (int j=0; i<4; ++j) {
                    args[4*i+j].f    = env->envcp[i][j];
                    arg_types[4*i+j] = 'f';
                }
            }
            d.replyArray(d.loc, arg_types, args);
        } else {
            for(int i=0; i<N && i<M; ++i) {
                for(int j=0; j<4; ++j) {
                    env->envcp[i][j] = limit(roundf(rtosc_argument(msg,4*i+j).f), 0.0f, 1.0f);
                }
            }
            env->updatenonfree();
            char part_loc[128];
            strncpy(part_loc, d.loc, sizeof(part_loc));
            part_loc[sizeof(part_loc) - 1] = '\0';
            char *end = strrchr(part_loc, '/');
            if(end) {
              end[1] = '\0';
              d.broadcast("/damage", "s", part_loc);
            }
        }
        rEnd},
        {"envval", rDoc("Envelope Values"), NULL,
        rBegin;
        const int N = MAX_ENVELOPE_POINTS;
        const int M = rtosc_narguments(msg);
        if(M == 0) {
            rtosc_arg_t args[N];
            char arg_types[N+1] = {};
            for(int i=0; i<N; ++i) {
                args[i].f    = env->Penvval[i]/127.0f;
                arg_types[i] = 'f';
            }
            d.replyArray(d.loc, arg_types, args);
        } else {
            for(int i=0; i<N && i<M; ++i) {
                env->Penvval[i] = limit(roundf(rtosc_argument(msg,i).f*127.0f), 0.0f, 127.0f);
            }
            env->updatenonfree();
              char part_loc[128];
              strncpy(part_loc, d.loc, sizeof(part_loc));
              part_loc[sizeof(part_loc) - 1] = '\0';
              char *end = strrchr(part_loc, '/');
              if(end) {
                  end[1] = '\0';
                  d.broadcast("/damage", "s", part_loc);
                }
        }
        rEnd},
    rArrayFF(envcp, MAX_ENVELOPE_POINTS-1),
    {"addPoint:i", rProp(internal) rDoc("Add point to envelope"), NULL,
        rBegin;
        const int curpoint = rtosc_argument(msg, 0).i;
        //int curpoint=freeedit->lastpoint;
        if (curpoint<0 || curpoint>env->Penvpoints || env->Penvpoints>=MAX_ENVELOPE_POINTS)
            return;
        int i;
        for (i=env->Penvpoints; i>=curpoint+1; i--) {
            env->envdt[i]=env->envdt[i-1];
            env->Penvval[i]=env->Penvval[i-1];
        }
        env->Penvval[i]=env->Penvval[i-1] + 0.2f * (env->Penvval[i+1] - env->Penvval[i-1]);
        
        if (curpoint==0)
            env->envdt[1]=dTREAL(64);

        env->Penvpoints++;
        if (curpoint<=env->Penvsustain)
            env->Penvsustain++;
        rEnd},
    {"delPoint:i", rProp(internal) rDoc("Delete Envelope Point"), NULL,
        rBegin;
        const int curpoint=rtosc_argument(msg, 0).i;
        if(curpoint<1 || curpoint>=env->Penvpoints-1 || env->Penvpoints<=3)
            return;

        for (int i=curpoint+1;i<env->Penvpoints;i++){
            env->envdt[i-1]=env->envdt[i];
            env->Penvval[i-1]=env->Penvval[i];
        };

        env->Penvpoints--;

        if (curpoint<=env->Penvsustain)
            env->Penvsustain--;

        rEnd},
};
#undef  rChangeCb

const rtosc::Ports &EnvelopeParams::ports = localPorts;


EnvelopeParams::EnvelopeParams(unsigned char Penvstretch_,
                               unsigned char Pforcedrelease_,
                               const AbsTime *time_):
        time(time_), last_update_timestamp(0)
{
    A_dt  = 0.009;
    D_dt  = 0.009;
    R_dt  = 0.009;
    PA_val = 64;
    PD_val = 64;
    PS_val = 64;
    PR_val = 64;

    for(int i = 0; i < MAX_ENVELOPE_POINTS; ++i) {
        envdt[i]  = dTREAL(32);
        Penvval[i] = 64;
    }
    envdt[0]        = 0.0f; //not used
    Penvsustain     = 1;
    Penvpoints      = 1;
    Envmode         = ADSR_lin;
    Penvstretch     = Penvstretch_;
    Pforcedrelease  = Pforcedrelease_;
    Pfreemode       = 1;
    Plinearenvelope = 0;
    Prepeating      = 0;

    store2defaults();
}

EnvelopeParams::~EnvelopeParams()
{}

#define COPY(y) this->y = ep.y
void EnvelopeParams::paste(const EnvelopeParams &ep)
{

    COPY(Pfreemode);
    COPY(Penvpoints);
    COPY(Penvsustain);
    for(int i=0; i<MAX_ENVELOPE_POINTS; ++i) {
        this->envdt[i]  = ep.envdt[i];
        this->Penvval[i] = ep.Penvval[i];
    }
    COPY(Penvstretch);
    COPY(Pforcedrelease);
    COPY(Plinearenvelope);

    COPY(A_dt);
    COPY(D_dt);
    COPY(R_dt);
    COPY(PA_val);
    COPY(PD_val);
    COPY(PS_val);
    COPY(PR_val);

    if ( time ) {
        last_update_timestamp = time->time();
    }
}
#undef COPY

void EnvelopeParams::init(zyn::consumer_location_t _loc)
{
    switch(loc = _loc)
    {
        case ad_global_amp:    ADSRinit_dB(0.0f, 0.127f, 127, 0.041f); break;
        case ad_global_freq:   ASRinit(64, 0.254f, 64, 0.499f); break;
        case ad_global_filter:
        case sub_filter:       ADSRinit_filter(64, 0.127f, 64, 0.970f, 0.499f, 64); break;
        case ad_voice_amp:     ADSRinit_dB(0.0f, 6.978f, 127, 6.978f); break;
        case ad_voice_freq:    ASRinit(30, 0.127f, 64, 0.499f); break;
        case ad_voice_filter:  ADSRinit_filter(90, 0.970f, 40, 0.970f, 0.009f, 40); break;
        case ad_voice_fm_freq: ASRinit(20, 3.620f, 40, 1.876f); break;
        case ad_voice_fm_amp:  ADSRinit(1.876f, 3.620f, 127, 6.978f); break;
        case sub_freq:         ASRinit(30, 0.254f, 64, 0.499f); break;
        case sub_bandwidth:    ASRinit_bw(100, 0.970f, 64, 0.499f); break;
        default: throw std::logic_error("Invalid envelope consumer location");
    };
}

float EnvelopeParams::getdt(char i) const
{
    return envdt[(int)i]; //seconds
}

/*
 * ADSR/ASR... initialisations
 */
void EnvelopeParams::ADSRinit(float a_dt, float d_dt, char s_val, float r_dt)
{
    setpresettype("Penvamplitude");
    Envmode   = ADSR_lin;
    A_dt      = a_dt;
    D_dt      = d_dt;
    PS_val    = s_val;
    R_dt      = r_dt;
    Pfreemode = 0;
    converttofree();

    store2defaults();
}

void EnvelopeParams::ADSRinit_dB(float a_dt, float d_dt, char s_val, float r_dt)
{
    setpresettype("Penvamplitude");
    Envmode   = ADSR_dB;
    A_dt      = a_dt;
    D_dt      = d_dt;
    PS_val    = s_val;
    R_dt      = r_dt;
    Pfreemode = 0;
    converttofree();

    store2defaults();
}

void EnvelopeParams::ASRinit(char a_val, float a_dt, char r_val, float r_dt)
{
    setpresettype("Penvfrequency");
    Envmode   = ASR_freqlfo;
    PA_val    = a_val;
    A_dt      = a_dt;
    PR_val    = r_val;
    R_dt      = r_dt;
    Pfreemode = 0;
    converttofree();

    store2defaults();
}

void EnvelopeParams::ADSRinit_filter(char a_val,
                                     float a_dt,
                                     char d_val,
                                     float d_dt,
                                     float r_dt,
                                     char r_val)
{
    setpresettype("Penvfilter");
    Envmode   = ADSR_filter;
    PA_val    = a_val;
    A_dt      = a_dt;
    PD_val    = d_val;
    D_dt      = d_dt;
    R_dt      = r_dt;
    PR_val    = r_val;
    Pfreemode = 0;
    converttofree();
    store2defaults();
}

void EnvelopeParams::ASRinit_bw(char a_val, float a_dt, char r_val, float r_dt)
{
    setpresettype("Penvbandwidth");
    Envmode   = ASR_bw;
    PA_val    = a_val;
    A_dt      = a_dt;
    PR_val    = r_val;
    R_dt      = r_dt;
    Pfreemode = 0;
    converttofree();
    store2defaults();
}

/*
 * Convert the Envelope to freemode
 */
void EnvelopeParams::converttofree()
{
    switch(Envmode) {
        case ADSR_lin:
        case ADSR_dB:
            Penvpoints  = 4;
            Penvsustain = 2;
            Penvval[0]  = 0;
            envdt[1]   = A_dt;
            Penvval[1]  = 127;
            envdt[2]   = D_dt;
            Penvval[2]  = PS_val;
            envdt[3]   = R_dt;
            Penvval[3]  = 0;
            break;
        case ASR_freqlfo:
        case ASR_bw:
            Penvpoints  = 3;
            Penvsustain = 1;
            Penvval[0]  = PA_val;
            envdt[1]   = A_dt;
            Penvval[1]  = 64;
            envdt[2]   = R_dt;
            Penvval[2]  = PR_val;
            break;
        case ADSR_filter:
            Penvpoints  = 4;
            Penvsustain = 2;
            Penvval[0]  = PA_val;
            envdt[1]   = A_dt;
            Penvval[1]  = PD_val;
            envdt[2]   = D_dt;
            Penvval[2]  = 64;
            envdt[3]   = R_dt;
            Penvval[3]  = PR_val;
            break;
    }
}

void EnvelopeParams::updatefree()
{
    Penvval[0]  = PA_val;
    envdt[1]    = A_dt;
    envdt[Penvsustain]   = D_dt;
    Penvval[Penvsustain-1]  = PD_val;
    Penvval[Penvsustain]  = PS_val;
    envdt[Penvpoints-1]   = R_dt;
    Penvval[Penvpoints-1]  = PR_val;
}

void EnvelopeParams::updatenonfree()
{
    PA_val = Penvval[0];
    A_dt = envdt[1];
    D_dt = envdt[Penvsustain];
    PD_val = Penvval[Penvsustain-1];
    PS_val = Penvval[Penvsustain];
    R_dt = envdt[Penvpoints-1];
    PR_val = Penvval[Penvpoints-1];
}


void EnvelopeParams::add2XML(XMLwrapper& xml)
{
    xml.addparbool("free_mode", Pfreemode);
    xml.addpar("env_points", Penvpoints);
    xml.addpar("env_sustain", Penvsustain);
    xml.addpar("env_stretch", Penvstretch);
    xml.addparbool("forced_release", Pforcedrelease);
    xml.addparbool("linear_envelope", Plinearenvelope);
    xml.addparbool("repeating_envelope", Prepeating);
    xml.addparreal("A_dt", A_dt);
    xml.addparreal("D_dt", D_dt);
    xml.addparreal("R_dt", R_dt);
    xml.addpar("A_val", PA_val);
    xml.addpar("D_val", PD_val);
    xml.addpar("S_val", PS_val);
    xml.addpar("R_val", PR_val);

    if((Pfreemode != 0) || (!xml.minimal))
        for(int i = 0; i < Penvpoints; ++i) {
            xml.beginbranch("POINT", i);
            if(i != 0)
                xml.addparreal("dt", envdt[i]);
            xml.addpar("val", Penvval[i]);
            xml.endbranch();
        }
}

float EnvelopeParams::env_dB2rap(float db) {
    return (powf(10.0f, db / 20.0f) - 0.01)/.99f;
}

float EnvelopeParams::env_rap2dB(float rap) {
    return 20.0f * log10f(rap * 0.99f + 0.01);
}

/**
    since commit 5334d94283a513ae42e472aa020db571a3589fb9, i.e. between
    versions 2.4.3 and 2.4.4, the amplitude envelope has been converted
    differently from dB to rap for AmplitudeEnvelope (mode 2)
    this converts the values read from an XML file once
*/
struct version_fixer_t
{
    const bool mismatch;
public:
    int operator()(int input) const
    {
        return (mismatch)
            // The errors occurred when calling env_dB2rap. Let f be the
            // conversion function for mode 2 (see Envelope.cpp), then we
            // load values with (let "o" be the function composition symbol):
            //   f^{-1} o (env_dB2rap^{-1}) o dB2rap o f
            // from the xml file. This results in the following formula:
            ? roundf(127.0f * (0.5f *
                   log10f( 0.01f + 0.99f *
                                       powf(100, input/127.0f - 1))
                               + 1))
            : input;
    }
    version_fixer_t(const version_type& fileversion, int env_mode) :
        mismatch(fileversion < version_type(2,4,4) &&
                 (env_mode == 2)) {}
};

void EnvelopeParams::getfromXML(XMLwrapper& xml)
{
    Pfreemode       = xml.getparbool("free_mode", Pfreemode);
    Penvpoints      = xml.getpar127("env_points", Penvpoints);
    Penvsustain     = xml.getpar127("env_sustain", Penvsustain);
    Penvstretch     = xml.getpar127("env_stretch", Penvstretch);
    Pforcedrelease  = xml.getparbool("forced_release", Pforcedrelease);
    Plinearenvelope = xml.getparbool("linear_envelope", Plinearenvelope);
    Prepeating      = xml.getparbool("repeating_envelope", Prepeating);

    version_fixer_t version_fix(xml.fileversion(), Envmode);

    if(!xml.hasparreal("A_dt")) {
        A_dt = dTREAL(xml.getpar127("A_dt", 0));
        D_dt = dTREAL(xml.getpar127("D_dt", 0));
        R_dt = dTREAL(xml.getpar127("R_dt", 0));
    } else {
        A_dt  = xml.getparreal("A_dt", A_dt);
        D_dt  = xml.getparreal("D_dt", D_dt);
        R_dt  = xml.getparreal("R_dt", R_dt);
    }

    PA_val = version_fix(xml.getpar127("A_val", PA_val));
    PD_val = version_fix(xml.getpar127("D_val", PD_val));
    PS_val = version_fix(xml.getpar127("S_val", PS_val));
    PR_val = version_fix(xml.getpar127("R_val", PR_val));

    for(int i = 0; i < Penvpoints; ++i) {
        if(xml.enterbranch("POINT", i) == 0)
            continue;
        if(i != 0) {
            if(!xml.hasparreal("dt")) {
                int dt = xml.getpar127("dt", dT127(envdt[i]));
                envdt[i] = dTREAL(dt);
            }
            else {
                envdt[i] = xml.getparreal("dt", envdt[i]);
            }
        }
        Penvval[i] = version_fix(xml.getpar127("val", Penvval[i]));
        xml.exitbranch();
    }

    if(!Pfreemode)
        converttofree();
}


void EnvelopeParams::defaults()
{
    Penvstretch     = Denvstretch;
    Pforcedrelease  = Dforcedrelease;
    Plinearenvelope = Dlinearenvelope;
    Prepeating      = Drepeating;
    A_dt     = DA_dt;
    D_dt     = DD_dt;
    R_dt     = DR_dt;
    PA_val    = DA_val;
    PD_val    = DD_val;
    PS_val    = DS_val;
    PR_val    = DR_val;
    Pfreemode = 0;
    converttofree();
}

void EnvelopeParams::store2defaults()
{
    Denvstretch     = Penvstretch;
    Dforcedrelease  = Pforcedrelease;
    Dlinearenvelope = Plinearenvelope;
    Drepeating      = Prepeating;
    DA_dt  = A_dt;
    DD_dt  = D_dt;
    DR_dt  = R_dt;
    DA_val = PA_val;
    DD_val = PD_val;
    DS_val = PS_val;
    DR_val = PR_val;
}

}
