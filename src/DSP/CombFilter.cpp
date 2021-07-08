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
    :Filter(srate, bufsize), memory(*alloc), gain(1.0f), type(Ftype)
{
    //worst case: looking back from smps[0] at 25Hz using higher order interpolation
    mem_size = (int)ceilf((float)samplerate/25.0) + buffersize + 2; // 2178 at 48000Hz and 256Samples
    input = (float*)memory.alloc_mem(mem_size*sizeof(float));
    output = (float*)memory.alloc_mem(mem_size*sizeof(float));
    memset(input, 0, mem_size*sizeof(float));
    memset(output, 0, mem_size*sizeof(float));
    
    //~ delay_smoothing.sample_rate(srate);
    //~ delay_smoothing.reset(samplerate/1000.0f);
    setfreq_and_q(Ffreq, Fq);
    settype(type);
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
    int poshi = (int)floorf(pos); // integer part
    float poslo = pos - (float) poshi; // decimal part
    // linear interpolation between samples
    return smp[poshi] + poslo * (smp[poshi+1]-smp[poshi]); 
}

//~ inline float sampleLagrange3o4p(float *smp, float pos) {
    //~ // 4-point, 3rd-order Lagrange (x-form)
    //~ float c0 = smp[0];
    //~ float c1 = smp[1] - 0.33333333333f*smp[-1] - 0.5f*smp[0] - 0.16666666666f*smp[2];
    //~ float c2 = 0.5f*(smp[-1]+smp[1]) - smp[0];
    //~ float c3 = 0.16666666666f*(smp[2]-smp[-1]) + 0.5f*(smp[0]-smp[1]);
    //~ return 16.0f*(((c3*pos+c2)*pos+c1)*pos+c0);
//~ }

void CombFilter::filterout(float *smp)
{
    //~ float delaybuf[buffersize];

    //~ if(delayfwd_smoothing.apply(delaybuf, buffersize, delay))
        //~ for(int i=0; i<buffersize; ++i)
            //~ delayfwdbuf[i] = delayfwd;

    // shift the buffer content to the left
    memmove(&input[0], &input[buffersize], (mem_size-buffersize)*sizeof(float));
    // copy new input samples to the right end of the buffer
    memcpy(&input[mem_size-buffersize], smp, buffersize*sizeof(float));
    for (int i = 0; i < buffersize; i ++)
    {
        // calculate the feedback sample positions in the buffer
        float pos = float(mem_size-buffersize+i)-delay;
        // add the fwd and bwd feedback samples to current sample
        smp[i] = smp[i]*gain + tanhX(
            gainfwd * sampleLerp( input, pos) - 
            gainbwd * sampleLerp(output, pos)); 
        // apply output gain
        smp[i] *= outgain;
    }
    // shift the buffer content one buffersize to the left
    memmove(&output[0], &output[buffersize], (mem_size-buffersize)*sizeof(float));
    // copy new output samples to the right end of the buffer
    memcpy(&output[mem_size-buffersize], smp, buffersize*sizeof(float));
}

void CombFilter::setfreq_and_q(float freq, float q)
{
    setfreq(freq);
    setq(q);
}

void CombFilter::setfreq(float freq)
{
    float ff = limit(freq, 25.0f, 40000.0f);
    delay = ((float)samplerate)/ff;
}

void CombFilter::setq(float q_)
{
    q = cbrtf(0.0015f*q_);
    settype(type);
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
