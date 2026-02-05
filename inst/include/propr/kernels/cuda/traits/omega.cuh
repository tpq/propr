#pragma once
#include <cub/cub.cuh>
#include <cuda_runtime.h>

#include <propr/utils/preprocessor.cuh>


namespace propr {
    namespace cuda {
            namespace traits {

                struct omega_global_config {
                    const static int NUM_ACC = 2;

                    constexpr static inline int BLK_M = 128;
                    constexpr static inline int BLK_K = 8;
                    constexpr static inline int TH_Y  = 8;
                    constexpr static inline int TH_X  = 8;

                    constexpr static inline cub::CacheLoadModifier  LoadModifer  = cub::LOAD_CG;
                    constexpr static inline cub::CacheStoreModifier StoreModifer = cub::STORE_CG;

                    PROPR_DEVICE 
                    static 
                    PROPR_INLINE 
                    void update(const int tidx,const int tidy, 
                                float acc[NUM_ACC][TH_Y][TH_X], 
                                const float& a,const float& b) {
                        float n = 2.0f * (a / (a + b + FLT_EPSILON)) * b;
                        acc[0][tidy][tidx] += n; 
                        acc[1][tidy][tidx] += n * n;
                    }

                    PROPR_DEVICE 
                    static 
                    PROPR_INLINE 
                    float finalize(const int tidx,const int tidy,
                                    const float acc[NUM_ACC][TH_Y][TH_X]) {
                        float n = acc[0][tidy][tidx];
                        float s = acc[1][tidy][tidx];
                        return n - s / (n + FLT_EPSILON);
                    }
                };


                struct omega_population_config {
                    const static int NUM_ACC = 1; // number of acummilators

                    const static int BLK_M = 128;
                    const static int BLK_K = 8;
                    const static int TH_Y  = 8;
                    const static int TH_X  = 8;


                    static const cub::CacheLoadModifier  LoadModifer  = cub::LOAD_CG;
                    static const cub::CacheStoreModifier StoreModifer = cub::STORE_CG;

                    PROPR_DEVICE 
                    static 
                    PROPR_INLINE 
                    void  update(const int tidx, const int tidy,
                                float acc[NUM_ACC][TH_Y][TH_X], 
                                const float& a, const float& b) {
                        float n = 2.0f * (a / (a + b + FLT_EPSILON)) * b;
                        acc[0][tidy][tidx] += n;
                    }

                    PROPR_DEVICE 
                    static 
                    PROPR_INLINE 
                    float finalize(const int tidx, const int tidy,
                                    const float acc[NUM_ACC][TH_Y][TH_X]) {
                        return acc[0][tidy][tidx];
                    }
                };


            }
    }
}