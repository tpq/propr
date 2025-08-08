#pragma once
#include <cute/tensor.hpp>
#include <cutlass/numeric_conversion.h>

#include <propr/kernels/cuda/detail/cutlass/omega/atoms/harmonic_atom_sm80.cuh>
#include <propr/kernels/cuda/detail/cutlass/omega/traits/harmonic_traits_sm80.cuh>


namespace propr {
    namespace kernels{
        namespace cutlass_impl {
            struct OmegaConfig {
                public:
                    static constexpr auto BLK_M = cute::Int<128>{};
                    static constexpr auto BLK_K = cute::Int<32>{};
                private:
                    using BlockShapeA    = cute::Shape<cute::Int<BLK_M>, cute::Int<BLK_K>>;
                    using BlockShapeB    = cute::Shape<cute::Int<BLK_M>, cute::Int<BLK_K>>;
                    using SmemLayoutAtom = decltype(cute::composition(cute::Swizzle<3, 3, 3>{},
                                                    cute::make_layout(cute::make_shape(cute::Int<8>{}, cute::Int<BLK_K>{}),
                                                    cute::make_stride(cute::Int<BLK_K>{}, cute::Int<1>{}))));

                public:
                    using SmemLayoutA = decltype(cute::tile_to_shape(SmemLayoutAtom{}, BlockShapeA{}));
                    using SmemLayoutB = decltype(cute::tile_to_shape(SmemLayoutAtom{}, BlockShapeB{}));
                    struct SharedStorage {
                        cute::ArrayEngine<cute::half_t, cute::cosize_v<SmemLayoutA>> A;
                        cute::ArrayEngine<cute::half_t, cute::cosize_v<SmemLayoutB>> B;
                    };

                private:
                    using GmemCopyTraits = cute::Copy_Traits<cute::SM80_CP_ASYNC_CACHEGLOBAL<cute::uint128_t>>;
                    using GmemCopyAtom   = cute::Copy_Atom<GmemCopyTraits, cute::half_t>;
                    using GmemCopyThreadLayout = decltype(cute::make_layout(
                                                             cute::make_shape (cute::Int<32>{}, cute::Int<4>{}),    
                                                             cute::make_stride(cute::Int<4>{} , cute::Int<1>{})));
                    using GmemCopyValLayout = decltype(cute::make_layout(cute::make_shape(cute::Int<1>{} , cute::Int<8>{})));
                public:
                    using GmemCopyA = decltype(cute::make_tiled_copy(GmemCopyAtom{}, GmemCopyThreadLayout{}, GmemCopyValLayout{}));
                    using GmemCopyB = GmemCopyA;

                private:
                    using SmemCopyAtom  = cute::Copy_Atom<cute::Copy_Traits<cute::SM75_U32x4_LDSM_N>, cute::half_t>;
                 // using SmemCopyAtom  = cute::Copy_Atom<cute::Copy_Traits<cute::SM80_CP_ASYNC_CACHEGLOBAL<cute::uint128_t>>, cute::half_t>;

                    using MmaTraits = cute::MMA_Traits<propr::cutlass_atoms::SM80_16x8x16_F32F16F16F32_TN_HARMONIC_1>;
                    using MmaAtom   = cute::MMA_Atom<MmaTraits>;

                    using MmaTraits2 = cute::MMA_Traits<propr::cutlass_atoms::SM80_16x8x16_F32F16F16F32_TN_HARMONIC_2>;
                    using MmaAtom2   = cute::MMA_Atom<MmaTraits2>;

                    static constexpr int kMmaEURepeatM = 2;
                    static constexpr int kMmaEURepeatN = 2;
                    static constexpr int kMmaEURepeatK = 1;
                    using MmaAtomShape = MmaTraits::Shape_MNK;

                    static constexpr int kMmaPM = 2 * kMmaEURepeatM * cute::get<0>(MmaAtomShape{});
                    static constexpr int kMmaPN = 2 * kMmaEURepeatN * cute::get<1>(MmaAtomShape{});
                    static constexpr int kMmaPK = 1 * kMmaEURepeatK * cute::get<2>(MmaAtomShape{});

                    using MMA_EU_RepeatT = decltype(cute::make_layout(
                                                        cute::make_shape(cute::Int<kMmaEURepeatM>{}, 
                                                                         cute::Int<kMmaEURepeatN>{}, 
                                                                         cute::Int<kMmaEURepeatK>{}
                                                                    )));
                    using MMA_P_T = cute::Tile<cute::Int<kMmaPM>, cute::Int<kMmaPN>, cute::Int<kMmaPK>>;
                public:
                    using TiledMMA   = decltype(cute::make_tiled_mma(MmaAtom{}, MMA_EU_RepeatT{}, MMA_P_T{}));
                    using TiledMMA2  = decltype(cute::make_tiled_mma(MmaAtom2{}, MMA_EU_RepeatT{}, MMA_P_T{}));

                    using SmemCopyA = decltype(cute::make_tiled_copy_A(SmemCopyAtom{}, TiledMMA{}));
                    using SmemCopyB = decltype(cute::make_tiled_copy_B(SmemCopyAtom{}, TiledMMA{}));
                
                public:
                    using NumericConverter     = decltype(cutlass::NumericConverter<cute::half_t, float, cutlass::FloatRoundStyle::round_to_nearest>());
                    using BackNumericConverter = decltype(cutlass::NumericConverter<float,cute::half_t , cutlass::FloatRoundStyle::round_to_nearest>());
                };
        }
    }
}