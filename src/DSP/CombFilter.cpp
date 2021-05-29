#include <cassert>
#include <cstdio>
#include <cmath>
#include <stdio.h>

#include "../Misc/Allocator.h"
#include "../Misc/Util.h"
#include "CombFilter.h"

// theory from `Introduction to Digital Filters with Audio Applications'', by Julius O. Smith III, (September 2007 Edition).
// https://www.dsprelated.com/freebooks/filters/Analysis_Digital_Comb_Filter.html

namespace zyn{

CombFilter::CombFilter(Allocator *alloc, unsigned char Ftype, float Ffreq, float Fq,
    unsigned int srate, int bufsize)
    :Filter(srate, bufsize), memory(*alloc), sr(srate), gain(1.0f), type(Ftype)
{
    mem_size = sr/30;
    input = (float*)memory.alloc_mem(mem_size);
    output = (float*)memory.alloc_mem(mem_size);
    memset(input, 0, sizeof(input));
    memset(output, 0, sizeof(output));
    
    delayfwd_smoothing.sample_rate(srate);
    delayfwd_smoothing.reset(sr/1000.0f);
    delaybwd_smoothing.sample_rate(srate);
    delaybwd_smoothing.reset(sr/1000.0f);
}

CombFilter::~CombFilter(void)
{
    memory.dealloc(input);
    memory.dealloc(output);
}


inline float CombFilter::tanhX(const float x)
{
    // Pade approximation of tanh(x) bound to [-1 .. +1]
    // https://mathr.co.uk/blog/2017-09-06_approximating_hyperbolic_tangent.html
    float x2 = x*x;
    return (x*(105.0f+10.0f*x2)/(105.0f+(45.0f+x2)*x2)); //
}

inline float CombFilter::sampleLerp(float *smp, float pos) {
    int poshi = (int)pos; // integer part
    float poslo = pos - (float) poshi; // decimal part
    // linear interpolation between samples
    return smp[poshi] + poslo * (smp[poshi+1]-smp[poshi]); 
}

inline float interp_cubic(float x0, float x1, float x2, float x3, float mu) {
   // http://paulbourke.net/miscellaneous/interpolation/
   float mu2 = mu*mu;
   float a0 = -0.5*x0 + 1.5*x1 - 1.5*x2 + 0.5*x3;
   float a1 = x0 - 2.5*x1 + 2*x2 - 0.5*x3;
   float a2 = -0.5*x0 + 0.5*x2;
   float a3 = x1;
   return(a0*mu*mu2+a1*mu2+a2*mu+a3);
}

inline float CombFilter::sampleHermite(float *smp, float pos) {
    int poshi = (int)pos; // integer part
    float poslo = pos - (float) poshi; // decimal part
    // linear interpolation between samples
    return smp[poshi] + poslo * (smp[poshi+1]-smp[poshi]); 
    interp_cubic(smp[poshi-1], smp[poshi], smp[poshi+1], smp[poshi+2], poslo); // smp[poshi+2] could be beyond borders at highes freq
}

void CombFilter::filterout(float *smp)
{
    float delayfwdbuf[buffersize];
    float delaybwdbuf[buffersize];
    delayfwd_smoothing.apply(delayfwdbuf, buffersize, delayfwd);
    delaybwd_smoothing.apply(delaybwdbuf, buffersize, delaybwd);

    memmove(&input[0], &input[buffersize-1], mem_size-buffersize);
    memmove(&input[mem_size-1-buffersize], smp, buffersize);
    for (int i = 0; i < buffersize; i ++)
    {
        smp[i] = smp[i]*gain + 
            gainfwd * sampleHermite( input, float(mem_size-buffersize+i)-delayfwdbuf[i]) + 
            gainbwd * sampleHermite(output, float(mem_size-buffersize+i)-delayfwdbuf[i]); 
        smp[i] *= outgain;
    }
    memmove(&output[0], &output[buffersize-1], sizeof(output)-buffersize);
    memmove(&output[sizeof(output)-1-buffersize], smp, buffersize);
    
}

void CombFilter::setfreq_and_q(float frequency, float q)
{
    setfreq(frequency/sr);
    setq(q);
}

void CombFilter::setfreq(float ff)
{
    delayfwd = sr/ff;
    delaybwd = sr/ff;
}

void CombFilter::setq(float q_)
{
    q = 10.0f*q_;
}

void CombFilter::setgain(float dBgain)
{
    gain = dB2rap(dBgain);
}

void CombFilter::settype(unsigned char type_)
{
    type = type_;
    switch (type)
    {
        case 0:
        default:
            gainfwd = 0.0f;
            gainbwd = q;
            break;
        case 1:
            gainfwd = q;
            gainbwd = 0.0f;
            break;
        case 2:
            gainfwd = q;
            gainbwd = q;
            break;
    }
}

};
