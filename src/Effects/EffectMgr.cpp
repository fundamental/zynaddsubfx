/*
  ZynAddSubFX - a software synthesizer

  EffectMgr.cpp - Effect manager, an interface between the program and effects
  Copyright (C) 2002-2005 Nasca Octavian Paul
  Author: Nasca Octavian Paul

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.
*/

#include <rtosc/ports.h>
#include <rtosc/port-sugar.h>
#include <iostream>
#include <cassert>

#include "EffectMgr.h"
#include "Effect.h"
#include "Alienwah.h"
#include "Reverb.h"
#include "Echo.h"
#include "Chorus.h"
#include "Distorsion.h"
#include "EQ.h"
#include "DynamicFilter.h"
#include "Phaser.h"
#include "../Misc/XMLwrapper.h"
#include "../Misc/Util.h"
#include "../Params/FilterParams.h"
#include "../Misc/Allocator.h"

namespace zyn {

#define rObject EffectMgr
#define rSubtype(name) \
    {STRINGIFY(name)"/", NULL, &name::ports,\
        [](const char *msg, rtosc::RtData &data){\
            rObject &o = *(rObject*)data.obj; \
            data.obj = dynamic_cast<name*>(o.efx); \
            if(!data.obj) \
                return; \
            SNIP \
            name::ports.dispatch(msg, data); \
        }}
static const rtosc::Ports local_ports = {
    rSelf(EffectMgr, rEnabledByCondition(self-enabled)),
    rPasteAndReply,
    rEnabledCondition(self-enabled, obj->geteffect()),
    rRecurp(filterpars, "Filter Parameter for Dynamic Filter"),
    {"Pvolume::i", rProp(parameter) rLinear(0,127) rShort("amt") rDoc("amount of effect"),
        0,
        [](const char *msg, rtosc::RtData &d)
        {
            EffectMgr *eff = (EffectMgr*)d.obj;
            if(!rtosc_narguments(msg))
                d.reply(d.loc, "i", eff->geteffectparrt(0));
            else if(rtosc_type(msg, 0) == 'i'){
                eff->seteffectparrt(0, rtosc_argument(msg, 0).i);
                d.broadcast(d.loc, "i", eff->geteffectparrt(0));
            }
        }},
    {"Ppanning::i", rProp(parameter) rLinear(0,127) rShort("pan") rDoc("panning"),
        0,
        [](const char *msg, rtosc::RtData &d)
        {
            EffectMgr *eff = (EffectMgr*)d.obj;
            if(!rtosc_narguments(msg))
                d.reply(d.loc, "i", eff->geteffectparrt(1));
            else if(rtosc_type(msg, 0) == 'i'){
                eff->seteffectparrt(1, rtosc_argument(msg, 0).i);
                d.broadcast(d.loc, "i", eff->geteffectparrt(1));
            }
        }},
    {"parameter#128::i:T:F", rProp(parameter) rProp(alias) rLinear(0,127) rDoc("Parameter Accessor"),
        NULL,
        [](const char *msg, rtosc::RtData &d)
        {
            EffectMgr *eff = (EffectMgr*)d.obj;
            const char *mm = msg;
            while(!isdigit(*mm))++mm;

            if(!rtosc_narguments(msg))
                d.reply(d.loc, "i", eff->geteffectparrt(atoi(mm)));
            else if(rtosc_type(msg, 0) == 'i'){
                eff->seteffectparrt(atoi(mm), rtosc_argument(msg, 0).i);
                d.broadcast(d.loc, "i", eff->geteffectparrt(atoi(mm)));
            } else if(rtosc_type(msg, 0) == 'T'){
                eff->seteffectparrt(atoi(mm), 127);
                d.broadcast(d.loc, "i", eff->geteffectparrt(atoi(mm)));
            } else if(rtosc_type(msg, 0) == 'F'){
                eff->seteffectparrt(atoi(mm), 0);
                d.broadcast(d.loc, "i", eff->geteffectparrt(atoi(mm)));
            }
        }},
    {"preset::i", rProp(parameter) rProp(alias) rDoc("Effect Preset Selector")
        rDefault(0), NULL,
        [](const char *msg, rtosc::RtData &d)
        {
            char loc[1024];
            EffectMgr *eff = (EffectMgr*)d.obj;
            if(!rtosc_narguments(msg))
                d.reply(d.loc, "i", eff->getpreset());
            else {
                eff->changepresetrt(rtosc_argument(msg, 0).i);
                d.broadcast(d.loc, "i", eff->getpreset());

                //update parameters as well
                fast_strcpy(loc, d.loc, sizeof(loc));
                char *tail = strrchr(loc, '/');
                if(!tail)
                    return;
                for(int i=0;i<128;++i) {
                    sprintf(tail+1, "parameter%d", i);
                    d.broadcast(loc, "i", eff->geteffectparrt(i));
                }
            }
        }},
    {"eq-coeffs:", rProp(internal) rDoc("Get equalizer Coefficients"), NULL,
        [](const char *, rtosc::RtData &d)
        {
            EffectMgr *eff = (EffectMgr*)d.obj;
            if(eff->nefx != 7)
                return;
            EQ *eq = (EQ*)eff->efx;
            float a[MAX_EQ_BANDS*MAX_FILTER_STAGES*3];
            float b[MAX_EQ_BANDS*MAX_FILTER_STAGES*3];
            memset(a, 0, sizeof(a));
            memset(b, 0, sizeof(b));
            eq->getFilter(a,b);
            d.reply(d.loc, "bb", sizeof(a), a, sizeof(b), b);
        }},
    {"efftype::i:c:S", rOptions(Disabled, Reverb, Echo, Chorus,
     Phaser, Alienwah, Distortion, EQ, DynFilter) rDefault(Disabled)
     rProp(parameter) rDoc("Get Effect Type"), NULL,
     rCOptionCb(obj->nefx, obj->changeeffectrt(var))},
    {"efftype:b", rProp(internal) rDoc("Pointer swap EffectMgr"), NULL,
        [](const char *msg, rtosc::RtData &d)
        {
            printf("OBSOLETE METHOD CALLED\n");
            EffectMgr *eff  = (EffectMgr*)d.obj;
            EffectMgr *eff_ = *(EffectMgr**)rtosc_argument(msg,0).b.data;

            //Lets trade data
            std::swap(eff->nefx,eff_->nefx);
            std::swap(eff->efx,eff_->efx);
            std::swap(eff->filterpars,eff_->filterpars);
            std::swap(eff->efxoutl, eff_->efxoutl);
            std::swap(eff->efxoutr, eff_->efxoutr);

            //Return the old data for destruction
            d.reply("/free", "sb", "EffectMgr", sizeof(EffectMgr*), &eff_);
        }},
    rSubtype(Alienwah),
    rSubtype(Chorus),
    rSubtype(Distorsion),
    rSubtype(DynamicFilter),
    rSubtype(Echo),
    rSubtype(EQ),
    rSubtype(Phaser),
    rSubtype(Reverb),
};

const rtosc::Ports &EffectMgr::ports = local_ports;

EffectMgr::EffectMgr(Allocator &alloc, const SYNTH_T &synth_,
                     const bool insertion_, const AbsTime *time_)
    :insertion(insertion_),
      efxoutl(new float[synth_.buffersize]),
      efxoutr(new float[synth_.buffersize]),
      filterpars(new FilterParams(in_effect, time_)),
      nefx(0),
      efx(NULL),
      time(time_),
      dryonly(false),
      memory(alloc),
      synth(synth_)
{
    setpresettype("Peffect");
    memset(efxoutl, 0, synth.bufferbytes);
    memset(efxoutr, 0, synth.bufferbytes);
    memset(settings, 255, sizeof(settings));
    defaults();
}


EffectMgr::~EffectMgr()
{
    memory.dealloc(efx);
    delete filterpars;
    delete [] efxoutl;
    delete [] efxoutr;
}

void EffectMgr::defaults(void)
{
    changeeffect(0);
    setdryonly(false);
}

//Change the effect
void EffectMgr::changeeffectrt(int _nefx, bool avoidSmash)
{
    cleanup();
    if(nefx == _nefx && efx != NULL)
        return;
    nefx = _nefx;
    preset = 0;
    memset(efxoutl, 0, synth.bufferbytes);
    memset(efxoutr, 0, synth.bufferbytes);
    memory.dealloc(efx);
    EffectParams pars(memory, insertion, efxoutl, efxoutr, 0,
            synth.samplerate, synth.buffersize, filterpars, avoidSmash);
    try {
        switch (nefx) {
            case 1:
                efx = memory.alloc<Reverb>(pars);
                break;
            case 2:
                efx = memory.alloc<Echo>(pars);
                break;
            case 3:
                efx = memory.alloc<Chorus>(pars);
                break;
            case 4:
                efx = memory.alloc<Phaser>(pars);
                break;
            case 5:
                efx = memory.alloc<Alienwah>(pars);
                break;
            case 6:
                efx = memory.alloc<Distorsion>(pars);
                break;
            case 7:
                efx = memory.alloc<EQ>(pars);
                break;
            case 8:
                efx = memory.alloc<DynamicFilter>(pars);
                break;
            //put more effect here
            default:
                efx = NULL;
                break; //no effect (thru)
        }
    } catch (std::bad_alloc &ba) {
        std::cerr << "failed to change effect " << _nefx << ": " << ba.what() << std::endl;
        return;
    }

    if(!avoidSmash)
        for(int i = 0; i != 128; i++)
            settings[i] = geteffectparrt(i);
}

void EffectMgr::changeeffect(int _nefx)
{
    nefx = _nefx;
    //preset    = 0;
}

//Obtain the effect number
int EffectMgr::geteffect(void)
{
    return nefx;
}

void EffectMgr::changesettingsrt(const short int *p_value)
{
    for(int i = 0; i != 128; i++) {
        short int value = p_value[i];
        /* check if setting is missing */
        if(value == -1) {
            if(efx)
                value = efx->getpresetpar(preset, i);
            else
                value = 0;
        }
        /* update settings */
        seteffectparrt(i, value);
    }
}

// Initialize An Effect in RT context
void EffectMgr::init(void)
{
    kill();
    changeeffectrt(nefx, true);
    changepresetrt(preset, true);
    changesettingsrt(settings);
}

//Strip effect manager of it's realtime memory
void EffectMgr::kill(void)
{
    //printf("Killing Effect(%d)\n", nefx);
    memory.dealloc(efx);
}

// Cleanup the current effect
void EffectMgr::cleanup(void)
{
    if(efx)
        efx->cleanup();
}


// Get the preset of the current effect
unsigned char EffectMgr::getpreset(void)
{
    if(efx)
        return efx->Ppreset;
    else
        return 0;
}

// Change the preset of the current effect
void EffectMgr::changepreset(unsigned char npreset)
{
    preset = npreset;
}

// Change the preset of the current effect
void EffectMgr::changepresetrt(unsigned char npreset, bool avoidSmash)
{
    preset = npreset;
    if(avoidSmash && dynamic_cast<DynamicFilter*>(efx)) {
        efx->Ppreset = npreset;
        return;
    }
    if(efx)
        efx->setpreset(npreset);
    if(!avoidSmash)
        for(int i = 0; i != 128; i++)
            settings[i] = geteffectparrt(i);
}

//Change a parameter of the current effect
void EffectMgr::seteffectparrt(int npar, unsigned char value)
{
    if(npar < 0 || npar >= 128)
        return;
    settings[npar] = value;

    if(!efx)
        return;
    try {
        efx->changepar(npar, value);
    } catch (std::bad_alloc &ba) {
        std::cerr << "failed to change effect parameter " << npar << " to " << value << ": " << ba.what() << std::endl;
    }
}

unsigned char EffectMgr::geteffectparrt(int npar)
{
    if(!efx)
        return 0;
    return efx->getpar(npar);
}

// Apply the effect
void EffectMgr::out(float *smpsl, float *smpsr)
{
    if(!efx) {
        if(!insertion)
            for(int i = 0; i < synth.buffersize; ++i) {
                smpsl[i]   = 0.0f;
                smpsr[i]   = 0.0f;
                efxoutl[i] = 0.0f;
                efxoutr[i] = 0.0f;
            }
        return;
    }
    for(int i = 0; i < synth.buffersize; ++i) {
        smpsl[i]  += synth.denormalkillbuf[i];
        smpsr[i]  += synth.denormalkillbuf[i];
        efxoutl[i] = 0.0f;
        efxoutr[i] = 0.0f;
    }
    efx->out(smpsl, smpsr);

    float volume = efx->volume;

    if(nefx == 7) { //this is need only for the EQ effect
        memcpy(smpsl, efxoutl, synth.bufferbytes);
        memcpy(smpsr, efxoutr, synth.bufferbytes);
        return;
    }

    //Insertion effect
    if(insertion != 0) {
        float v1, v2;
        if(volume < 0.5f) {
            v1 = 1.0f;
            v2 = volume * 2.0f;
        }
        else {
            v1 = (1.0f - volume) * 2.0f;
            v2 = 1.0f;
        }
        if((nefx == 1) || (nefx == 2))
            v2 *= v2;  //for Reverb and Echo, the wet function is not liniar

        if(dryonly)   //this is used for instrument effect only
            for(int i = 0; i < synth.buffersize; ++i) {
                smpsl[i]   *= v1;
                smpsr[i]   *= v1;
                efxoutl[i] *= v2;
                efxoutr[i] *= v2;
            }
        else // normal instrument/insertion effect
            for(int i = 0; i < synth.buffersize; ++i) {
                smpsl[i] = smpsl[i] * v1 + efxoutl[i] * v2;
                smpsr[i] = smpsr[i] * v1 + efxoutr[i] * v2;
            }
    }
    else // System effect
        for(int i = 0; i < synth.buffersize; ++i) {
            efxoutl[i] *= 2.0f * volume;
            efxoutr[i] *= 2.0f * volume;
            smpsl[i]    = efxoutl[i];
            smpsr[i]    = efxoutr[i];
        }
}


// Get the effect volume for the system effect
float EffectMgr::sysefxgetvolume(void)
{
    return efx ? efx->outvolume : 1.0f;
}


// Get the EQ response
float EffectMgr::getEQfreqresponse(float freq)
{
    return (nefx == 7) ? efx->getfreqresponse(freq) : 0.0f;
}


void EffectMgr::setdryonly(bool value)
{
    dryonly = value;
}

void EffectMgr::paste(EffectMgr &e)
{
    changeeffectrt(e.nefx, true);
    changepresetrt(e.preset, true);
    changesettingsrt(e.settings);
    if(dynamic_cast<DynamicFilter*>(efx)) {
        std::swap(filterpars, e.filterpars);
        efx->filterpars = filterpars;
    }
    cleanup(); // cleanup the effect and recompute its parameters
}

void EffectMgr::add2XML(XMLwrapper& xml)
{
    xml.addpar("type", geteffect());

    if(!geteffect())
        return;
    xml.addpar("preset", preset);

    xml.beginbranch("EFFECT_PARAMETERS");
    for(int n = 0; n != 128; n++) {
        int par;
        int def;
        if(efx) {
            par = efx->getpar(n);
            def = efx->getpresetpar(preset, n);
        } else {
            par = settings[n];
            def = -1;
        }
        /* don't store default values */
        if(par == def)
            continue;
        xml.beginbranch("par_no", n);
        xml.addpar("par", par);
        xml.endbranch();
    }
    assert(filterpars);
    if(nefx == 8) {
        xml.beginbranch("FILTER");
        filterpars->add2XML(xml);
        xml.endbranch();
    }
    xml.endbranch();
}

void EffectMgr::getfromXML(XMLwrapper& xml)
{
    changeeffect(xml.getpar127("type", geteffect()));

    if(!geteffect())
        return;

    preset = xml.getpar127("preset", preset);

    if(xml.enterbranch("EFFECT_PARAMETERS")) {
        for(int n = 0; n != 128; n++) {
            if(xml.enterbranch("par_no", n) == 0) {
                /*
                 * XXX workaround for old presets:
                 *
                 * All effect parameters have a default value.
                 * Default values are skipped when storing parameters,
                 * and must appear as the default value when loading.
                 *
                 * Up until recently it was assumed that the default
                 * value of all parameters is zero. This is no longer
                 * true, but when loading old presets we need to
                 * preserve this behaviour! Else sounds may change.
                 */
                if (xml.fileversion() < version_type(3,0,6) &&
                    /* XXX old presets don't have DC offset */
                    (geteffect() != 6 || n < 11)) {
                        settings[n] = 0;
                } else {
                        settings[n] = -1; /* use parameter default */
                }
            } else {
                settings[n] = xml.getpar127("par", 0);
                xml.exitbranch();
            }
        }
        assert(filterpars);
        if(xml.enterbranch("FILTER")) {
            filterpars->getfromXML(xml);
            xml.exitbranch();
        }
        xml.exitbranch();
    }
    cleanup();
}

}
