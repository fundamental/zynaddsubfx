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
    //worst case: looking back from smps[0] at 25Hz using higher order interpolation
    mem_size = sr/25 + buffersize+2; 
    input = (float*)memory.alloc_mem(mem_size);
    output = (float*)memory.alloc_mem(mem_size);
    memset(input, 0, mem_size);
    memset(output, 0, mem_size);
    
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

inline float sampleLagrange3o4p(float *smp, float pos) {
    // 4-point, 3rd-order Lagrange (x-form)
    float c0 = smp[0];
    float c1 = smp[1] - 0.33333333333f*smp[-1] - 0.5f*smp[0] - 0.16666666666f*smp[2];
    float c2 = 0.5f*(smp[-1]+smp[1]) - smp[0];
    float c3 = 0.16666666666f*(smp[2]-smp[-1]) + 0.5f*(smp[0]-smp[1]);
    return 16.0f*(((c3*pos+c2)*pos+c1)*pos+c0);
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
            gainfwd * sampleLagrange3o4p( input, float(mem_size-buffersize+i)-delayfwdbuf[i]) + 
            gainbwd * sampleLagrange3o4p(output, float(mem_size-buffersize+i)-delayfwdbuf[i]); 
        smp[i] *= outgain;
    }
    memmove(&output[0], &output[buffersize-1], mem_size-buffersize);
    memmove(&output[mem_size-1-buffersize], smp, buffersize);
    
}

void CombFilter::setfreq_and_q(float frequency, float q)
{
    setfreq(frequency/sr);
    setq(q);
}

void CombFilter::setfreq(float freq)
{
    float ff = limit(freq, 25.0f, 40000.0f);
    delayfwd = sr/ff;
    delaybwd = sr/ff;
}

void CombFilter::setq(float q_)
{
    q = 100.0f*q_;
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
