/*
  ZynAddSubFX - a software synthesizer

  Moog Filter.h - moog style multimode filter (lowpass, highpass...)
  Copyright (C) 2020-2020 Michael Kirchner
  Author: Michael Kirchner

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.
*/

#pragma once
#include "Filter.h"

namespace zyn {

class MoogFilter:public Filter
{
    public:
        //! @param Fq resonance, range [0.1,1000], logscale
        MoogFilter(unsigned char Ftype, float Ffreq, float Fq,
                unsigned int srate, int bufsize);
        ~MoogFilter() override;
        void filterout(float *smp) override;
        void setfreq(float /*frequency*/) override;
        void setfreq_and_q(float frequency, float q_) override;
        void setq(float /*q_*/) override;
        void setgain(float dBgain) override;
        void settype(unsigned char ftype);

    private:
        unsigned sr;
        float gain;

        float step(float x);

        float tanhXdivX(const float x) const;
        float tanhX(const float x) const;
        float tan_2(const float x) const;

        float feedbackGain;
        // aN: multimode coefficients for LP,BP,HP configurations
        float a0, a1, a2, a3, a4;
        float state[4] = {0.0f,0.0f,0.0f,0.0f};
        float passbandCompensation;
        // c*: cutoff frequency c and precalced products (times t) and powers(p) of it
        float c, ct2, cp2, cp3, cp4;
};

}
