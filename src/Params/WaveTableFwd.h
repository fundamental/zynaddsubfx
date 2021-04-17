/*
  ZynAddSubFX - a software synthesizer

  WaveTable.h - WaveTable forward declarations
  Copyright (C) 2020-2020 Johannes Lorenz
  Author: Johannes Lorenz

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.
*/

#ifndef WAVETABLEFWD_H
#define WAVETABLEFWD_H

#include <cstddef>
#include <cstdint>

namespace zyn
{

typedef uint32_t prng_t;

namespace wavetable_types
{
    using float32 = float;
    static_assert (sizeof (float32) == 4, "float does not have size 4");

    union IntOrFloat
    {
        prng_t intVal;
        float32 floatVal;
    };

    enum class WtMode
    {
        freq_smps, // (freq)->samples
        freqseed_smps, // (seed, freq)->samples
        freqwave_smps // (wave param, freq)->samples
    };
}

using tensor_size_t = std::size_t;

//! Tensor class for all dimensions != 1
template <std::size_t N, class T>
class Tensor;

//! Tensor class for dimension 1
template <class T>
class Tensor<1, T>;

template<class T> using Tensor1 = Tensor<1, T>;
template<class T> using Tensor2 = Tensor<2, T>;
template<class T> using Tensor3 = Tensor<3, T>;

class WaveTable;

}

#endif // WAVETABLEFWD_H
