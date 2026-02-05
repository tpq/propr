#pragma once

#include <cute/tensor.hpp>
#include <cute/atom/mma_traits.hpp>

#include <propr/kernels/cuda/detail/cutlass/omega/atoms/harmonic_atom_sm80.cuh>
#include <propr/utils/constants.h>
#include <propr/utils/preprocessor.cuh>


namespace cute {

    template <>
    struct MMA_Traits<propr::cutlass_atoms::SM80_16x8x16_F32F16F16F32_TN_HARMONIC_1>
    : cute::MMA_Traits<cute::SM80_16x8x16_F32F16F16F32_TN>
    {
        using ValTypeD = float;
        using ValTypeA = cute::half_t;
        using ValTypeB = cute::half_t;
        using ValTypeC = float;
    };

    template <>
    struct MMA_Traits<propr::cutlass_atoms::SM80_16x8x16_F32F16F16F32_TN_HARMONIC_2>
    : cute::MMA_Traits<cute::SM80_16x8x16_F32F16F16F32_TN>
    {
        using ValTypeD = float;
        using ValTypeA = cute::half_t;
        using ValTypeB = cute::half_t;
        using ValTypeC = float;
    };
}