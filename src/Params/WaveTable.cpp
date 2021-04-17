/*
  ZynAddSubFX - a software synthesizer

  WaveTable.cpp - WaveTable implementation
  Copyright (C) 2020-2020 Johannes Lorenz

  This program is free software; you can redistribute it and/or modify
  it under the terms of version 2 of the GNU General Public License
  as published by the Free Software Foundation.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License (version 2 or later) for more details.

  You should have received a copy of the GNU General Public License (version 2)
  along with this program; if not, write to the Free Software Foundation,
  Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA

*/

#include <cmath>
#include "WaveTable.h"

namespace zyn {

const Tensor1<WaveTable::float32>& WaveTable::get(float32 freq)
{
    tensor_size_t num_freqs = freqs.size(), i;
    float bestFreqDist = 9999999;
    tensor_size_t bestI = 0;
    for(i = 0; i < num_freqs; ++i) // TODO: std::lower_bound?
    {
        float curFreqDist = fabsf(freqs[i] - freq);
        if(bestFreqDist > curFreqDist)
        {
            bestFreqDist = curFreqDist;
            bestI = i;
        }
        else
        {
            assert(i > 0);
            bestI = i - 1; // it got worse, so take the previous freq
            break;
        }
    }
    assert(bestI < num_freqs);

    AbstractRingbuffer& rb = data.ringbuffers[bestI];
    assert(data[rb.read_pos()].size() == freqs.size());


    // TODO: hide this logic in WaveTable struct
    const Tensor1<WaveTable::float32>& res =
        (rb.read_space() == 0)
            // previous element
            ? data[(rb.read_pos() + data.size() - 1) % data.size()][bestI]
            // normal case
            : data[rb.read_pos()][bestI];
    if(mode() == WtMode::freqseed_smps && rb.read_space()) {
#ifdef DBG_WAVETABLES
        rb.dump();
#endif
        rb.inc_read_pos();
    }
    else {
        // there was no read space (we just read the previous
        // element),
        // or we just do not consume any seeds in general (peak at read pos)
    }
    return res;
}

WaveTable::WaveTable(tensor_size_t buffersize) :
    semantics(Shape1{1}, Shape1{0}),
    freqs(Shape1{1}, Shape1{0}),
    data(Shape3{1, 1, buffersize}, Shape3{0, 0, 0})
{
    setMode(WtMode::freq_smps);
}

WaveTable::WaveTable(tensor_size_t nsemantics, tensor_size_t nfreqs) :
    semantics(Shape1{0}, Shape1{0}),
    freqs(Shape1{0}, Shape1{0}),
    data(Shape3{nsemantics, nfreqs, 0}, Shape3{0, 0, 0})
{
    setMode(WtMode::freq_smps);
}

}

