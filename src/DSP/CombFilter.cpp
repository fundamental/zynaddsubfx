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
    int poshi = (int)floorf(pos);
    float poslo = pos - (float) poshi;
    return (1.0f - poslo) * smp[poshi] + poslo * smp[poshi+1];
}

void CombFilter::filterout(float *smp)
{
    memmove(&input[0], &input[buffersize-1], sizeof(input)-buffersize);
    memmove(&input[sizeof(input)-1-buffersize], smp, buffersize);
    for (int i = 0; i < buffersize; i ++)
    {
        smp[i] = smp[i]*gain + 
            gainfwd * sampleLerp(input, i-delayfwd) + 
            gainbwd * sampleLerp(output, i-delaybwd); 
        smp[i] *= outgain;
    }
    memmove(&output[0], &output[buffersize-1], sizeof(output)-buffersize);
    memmove(&output[sizeof(output)-1-buffersize], smp, buffersize);
    
}

void CombFilter::setfreq_and_q(float frequency, float q_)
{
    setfreq(frequency/sr);
    setq(q_);
}

void CombFilter::setfreq(float ff)
{
    delayfwd = sr/ff;
    delaybwd = sr/ff;
}

void CombFilter::setq(float q_)
{
    q = q_;
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
