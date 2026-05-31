#pragma once

#include <Rcpp.h>
#include <cuda_runtime.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <type_traits>
#include <vector>

#include <propr/context.h>
#include <propr/data/math.cuh>
#include <propr/utils/common/preprocessor.cuh>
#include <propr/utils/cuda/cuda_checks.h>

namespace propr {
    namespace detail {
        namespace cuda {

            enum class threshold_direction {less, greater,};

            template <class Config>
            struct threshold_cutoff {
                using cutoff_t = typename Config::cutoff_t;

                cutoff_t value;
                R_xlen_t index;
            };

            template <class Config>
            struct uniform_cutoffs {
                using cutoff_t = typename Config::cutoff_t;

                bool valid = false;
                cutoff_t first = static_cast<cutoff_t>(0);
                cutoff_t step  = static_cast<cutoff_t>(1);
                cutoff_t inv_step = static_cast<cutoff_t>(1);
            };

            template <class Config>
            struct cuda_threshold_group {
                using cutoff_t = typename Config::cutoff_t;
                using accumulator_t = typename Config::accumulator_t;

                std::vector<threshold_cutoff<Config>> cutoffs;
                std::vector<cutoff_t> values;
                uniform_cutoffs<Config> uniform;
                cutoff_t* d_cutoffs = nullptr;
                accumulator_t* d_buckets = nullptr;
            };

            template <class Config>
            struct cuda_threshold_counter {
                R_xlen_t ncutoffs = 0;
                cuda_threshold_group<Config> less;
                cuda_threshold_group<Config> greater;
            };

            template <class Config>
            inline void validate_threshold_config() {
                static_assert(std::is_floating_point_v<typename Config::value_t>,  "Config::value_t must be a floating point type.");
                static_assert(std::is_floating_point_v<typename Config::cutoff_t>, "Config::cutoff_t must be a floating point type.");
                static_assert(std::is_unsigned_v<typename Config::block_count_t>,  "Config::block_count_t must be an unsigned integer type.");
                static_assert(std::is_unsigned_v<typename Config::accumulator_t>,  "Config::accumulator_t must be an unsigned integer type.");

                static_assert(Config::BLK_X > 0,                                   "Config::BLK_X must be positive.");
                static_assert(Config::BLOCKS_PER_SM > 0,                           "Config::BLOCKS_PER_SM must be positive.");
                static_assert(Config::SHARED_CUTOFF_PAD_INTERVAL > 0,              "Config::SHARED_CUTOFF_PAD_INTERVAL must be positive.");
            }

            template <class Config>
            inline uniform_cutoffs<Config> 
            describe_uniform_cutoffs(const std::vector<typename Config::cutoff_t>& cutoffs) {
                using cutoff_t = typename Config::cutoff_t;
                uniform_cutoffs<Config> uniform{};

                if (cutoffs.empty()) return uniform;

                uniform.first = cutoffs[0];

                if (cutoffs.size() == 1) {
                    uniform.valid = true;
                    return uniform;
                }

                uniform.step = cutoffs[1] - cutoffs[0];
                if (!(std::isfinite(static_cast<double>(uniform.step)) && uniform.step > static_cast<cutoff_t>(0))) {
                    return uniform;
                }

                uniform.inv_step = static_cast<cutoff_t>(1) / uniform.step;

                for (size_t i = 2; i < cutoffs.size(); ++i) {
                    const cutoff_t expected = uniform.first + uniform.step * static_cast<cutoff_t>(i);
                    if (!propr::math::nearly_equal( static_cast<double>(cutoffs[i]), static_cast<double>(expected))) {
                        return uniform;
                    }
                }

                uniform.valid = true;
                return uniform;
            }

            inline size_t align_up_bytes(size_t bytes, size_t alignment) {
                return ((bytes + alignment - 1) / alignment) * alignment;
            }

            template <class Config>
            inline size_t bucket_shared_bytes(int ncutoffs) {
                return static_cast<size_t>(ncutoffs + 1) * sizeof(typename Config::accumulator_t);
            }

            template <class Config>
            inline size_t bucket_and_cutoff_shared_bytes(int ncutoffs) {
                return bucket_shared_bytes<Config>(ncutoffs) + static_cast<size_t>(ncutoffs) * sizeof(typename Config::cutoff_t);
            }

            template <class Config>
            inline size_t bucket_and_cutoff_block_shared_bytes(int ncutoffs) {
                using block_count_t = typename Config::block_count_t;
                using cutoff_t = typename Config::cutoff_t;
                const size_t bucket_bytes = static_cast<size_t>(ncutoffs + 1) * sizeof(block_count_t);
                const size_t cutoff_offset = align_up_bytes(bucket_bytes, alignof(cutoff_t));
                return cutoff_offset + static_cast<size_t>(ncutoffs) * sizeof(cutoff_t);
            }

            template <class Config>
            inline size_t bucket_and_cutoff_block_padded_shared_bytes(int ncutoffs) {
                using block_count_t = typename Config::block_count_t;
                using cutoff_t = typename Config::cutoff_t;
                const size_t bucket_bytes = static_cast<size_t>(ncutoffs + 1) * sizeof(block_count_t);
                const size_t cutoff_offset = align_up_bytes(bucket_bytes, alignof(cutoff_t));
                const size_t padded_cutoffs = ncutoffs == 0 ? 0  : static_cast<size_t>(ncutoffs) +
                                                                   static_cast<size_t>(ncutoffs - 1) / Config::SHARED_CUTOFF_PAD_INTERVAL;
                return cutoff_offset + padded_cutoffs * sizeof(cutoff_t);
            }

            // begin: move outside this file
            inline size_t dynamic_shared_capacity(const cudaDeviceProp& prop) {
                return std::max(
                    static_cast<size_t>(prop.sharedMemPerBlock),
                    static_cast<size_t>(prop.sharedMemPerBlockOptin));
            }

            template <typename Kernel>
            inline void prepare_dynamic_shared_memory(
                Kernel kernel,
                size_t shared_bytes,
                const cudaDeviceProp& prop) {
                if (shared_bytes <= static_cast<size_t>(prop.sharedMemPerBlock)) return;

                PROPR_CUDA_CHECK(cudaFuncSetAttribute(kernel, cudaFuncAttributeMaxDynamicSharedMemorySize,  static_cast<int>(shared_bytes)));
                PROPR_CUDA_CHECK(cudaFuncSetAttribute(kernel, cudaFuncAttributePreferredSharedMemoryCarveout, cudaSharedmemCarveoutMaxShared));
            }
            // end: move outside this file

            template <class Config>
            inline void 
            prepare_group(cuda_threshold_group<Config>& group, cudaStream_t stream) {
                validate_threshold_config<Config>();
                using cutoff_t = typename Config::cutoff_t;
                using accumulator_t = typename Config::accumulator_t;

                std::stable_sort(group.cutoffs.begin(), group.cutoffs.end(),
                                 [](const threshold_cutoff<Config>& a, const threshold_cutoff<Config>& b) {
                                     return a.value < b.value;
                                 });

                group.values.resize(group.cutoffs.size());
                for (size_t i = 0; i < group.cutoffs.size(); ++i) {
                    group.values[i] = group.cutoffs[i].value;
                }

                group.uniform = describe_uniform_cutoffs<Config>(group.values);

                if (group.cutoffs.empty()) return;

                PROPR_CUDA_CHECK(cudaMalloc(&group.d_cutoffs, group.values.size() * sizeof(cutoff_t)));
                PROPR_CUDA_CHECK(cudaMemcpyAsync(group.d_cutoffs, group.values.data(),
                                                 group.values.size() * sizeof(cutoff_t),
                                                 cudaMemcpyHostToDevice,stream));

                PROPR_CUDA_CHECK(cudaMalloc( &group.d_buckets, (group.values.size() + 1) * sizeof(accumulator_t)));
                PROPR_CUDA_CHECK(cudaMemsetAsync(group.d_buckets, 0, (group.values.size() + 1) * sizeof(accumulator_t), stream));
            }

             // begin: move outside this file
            template <class Config>
            inline void 
            destroy_group(cuda_threshold_group<Config>& group) {
                if (group.d_cutoffs != nullptr) {
                    PROPR_CUDA_CHECK(cudaFree(group.d_cutoffs));
                    group.d_cutoffs = nullptr;
                }

                if (group.d_buckets != nullptr) {
                    PROPR_CUDA_CHECK(cudaFree(group.d_buckets));
                    group.d_buckets = nullptr;
                }
            }

            template <class Config>
            PROPR_DEVICE 
            PROPR_FORCE_INLINE
            int shared_cutoff_padded_index(int i) {
                return i + i / Config::SHARED_CUTOFF_PAD_INTERVAL;
            }

            template <typename cutoff_t, typename value_t>
            PROPR_DEVICE 
            PROPR_FORCE_INLINE
            int upper_bound_cutoff(
                const cutoff_t* __restrict__ cutoffs,
                int ncutoffs,
                value_t value) {
                int lo = 0;
                int hi = ncutoffs;

                while (lo < hi) {
                    const int mid = lo + ((hi - lo) >> 1);
                    if (value < cutoffs[mid]) {
                        hi = mid;
                    } else {
                        lo = mid + 1;
                    }
                }

                return lo;
            }

            template <typename cutoff_t, typename value_t>
            PROPR_DEVICE 
            PROPR_FORCE_INLINE
            int lower_bound_cutoff(
                const cutoff_t* __restrict__ cutoffs,
                int ncutoffs,
                value_t value) {
                int lo = 0;
                int hi = ncutoffs;

                while (lo < hi) {
                    const int mid = lo + ((hi - lo) >> 1);
                    if (cutoffs[mid] < value) {
                        lo = mid + 1;
                    } else {
                        hi = mid;
                    }
                }

                return lo;
            }

            template <class Config, typename cutoff_t, typename value_t>
            PROPR_DEVICE 
            PROPR_FORCE_INLINE
            int lower_bound_padded_cutoff(const cutoff_t* __restrict__ cutoffs,
                                          int ncutoffs, value_t value) {
                int lo = 0, hi = ncutoffs;
                while (lo < hi) {
                    const int mid = lo + ((hi - lo) >> 1);
                    if (cutoffs[shared_cutoff_padded_index<Config>(mid)] < value) {
                        lo = mid + 1;
                    } else {
                        hi = mid;
                    }
                }
                return lo;
            }

            template <class Config, typename cutoff_t, typename value_t>
            PROPR_DEVICE 
            PROPR_FORCE_INLINE
            int upper_bound_padded_cutoff(const cutoff_t* __restrict__ cutoffs, int ncutoffs, value_t value) {
                int lo = 0, hi = ncutoffs;
                while (lo < hi) {
                    const int mid = lo + ((hi - lo) >> 1);
                    if (value < cutoffs[shared_cutoff_padded_index<Config>(mid)]) {
                        hi = mid;
                    } else {
                        lo = mid + 1;
                    }
                }
                return lo;
            }
            // end: move outside this file

            template <class Config, threshold_direction direction>
            PROPR_DEVICE 
            PROPR_FORCE_INLINE
            int bucket_for_direction(
                typename Config::value_t value,
                const typename Config::cutoff_t* __restrict__ cutoffs,
                int ncutoffs) {
                if (propr::math::is_nan(value)) return direction == threshold_direction::less ? ncutoffs : 0;
                if constexpr (direction == threshold_direction::greater) {
                    return lower_bound_cutoff(cutoffs, ncutoffs, value);
                } else {
                    return upper_bound_cutoff(cutoffs, ncutoffs, value);
                }
            }

            template <class Config, threshold_direction direction>
            PROPR_DEVICE 
            PROPR_FORCE_INLINE
            int bucket_for_padded_cutoffs_direction(typename Config::value_t value,
                                                    const typename Config::cutoff_t* __restrict__ cutoffs, int ncutoffs) {
                if (propr::math::is_nan(value)) return direction == threshold_direction::less ? ncutoffs : 0;

                if constexpr (direction == threshold_direction::greater) {
                    return lower_bound_padded_cutoff<Config>(cutoffs, ncutoffs, value);
                } else {
                    return upper_bound_padded_cutoff<Config>(cutoffs, ncutoffs, value);
                }
            }

            template <class Config>
            PROPR_DEVICE 
            PROPR_FORCE_INLINE
            int uniform_bucket_less(typename Config::value_t value,
                                    typename Config::cutoff_t first,
                                    typename Config::cutoff_t step,
                                    typename Config::cutoff_t inv_step,
                                    int ncutoffs) {
                if (propr::math::is_nan(value)) return ncutoffs;
                if (ncutoffs == 1) return value < first ? 0 : 1;

                const auto scaled = (value - first) * inv_step;
                int bucket = static_cast<int>(floor(static_cast<double>(scaled))) + 1;

                if (bucket < 0) return 0;
                if (bucket > ncutoffs) return ncutoffs;

                while (bucket < ncutoffs && value >= first + step * static_cast<typename Config::cutoff_t>(bucket)) {
                    ++bucket;
                }

                while (bucket > 0 && value < first + step * static_cast<typename Config::cutoff_t>(bucket - 1)) {
                    --bucket;
                }

                return bucket;
            }

            template <class Config>
            PROPR_DEVICE 
            PROPR_FORCE_INLINE
            int uniform_bucket_greater(
                typename Config::value_t value,
                typename Config::cutoff_t first,
                typename Config::cutoff_t step,
                typename Config::cutoff_t inv_step,
                int ncutoffs) {
                if (propr::math::is_nan(value)) return 0;
                if (ncutoffs == 1) return value <= first ? 0 : 1;

                const auto scaled = (value - first) * inv_step;
                int bucket = static_cast<int>(ceil(static_cast<double>(scaled)));

                if (bucket < 0) return 0;
                if (bucket > ncutoffs) return ncutoffs;

                while (bucket < ncutoffs && value > first + step * static_cast<typename Config::cutoff_t>(bucket)) {
                    ++bucket;
                }

                while (bucket > 0 && value <= first + step * static_cast<typename Config::cutoff_t>(bucket - 1)) {
                    --bucket;
                }

                return bucket;
            }

            template <class Config, threshold_direction direction>
            PROPR_DEVICE 
            PROPR_FORCE_INLINE
            int uniform_bucket_for_direction(typename Config::value_t value,
                                             typename Config::cutoff_t first,
                                             typename Config::cutoff_t step,
                                             typename Config::cutoff_t inv_step,
                                             int ncutoffs) {
                if constexpr (direction == threshold_direction::greater) {
                    return uniform_bucket_greater<Config>(value, first, step, inv_step, ncutoffs);
                } else {
                    return uniform_bucket_less<Config>(value, first, step, inv_step, ncutoffs);
                }
            }

            template <class Config, threshold_direction direction>
            __global__ void bucket_counts_shared_cutoffs_shared_buckets_kernel(
                const typename Config::value_t* __restrict__ values,   size_t nvalues,
                const typename Config::cutoff_t* __restrict__ cutoffs, int    ncutoffs,
                typename Config::accumulator_t* __restrict__ buckets) {
                using accumulator_t = typename Config::accumulator_t;
                using cutoff_t = typename Config::cutoff_t;

                extern __shared__ accumulator_t shared_storage[];

                accumulator_t* shared_buckets = shared_storage;
                cutoff_t* shared_cutoffs =
                    reinterpret_cast<cutoff_t*>(shared_buckets + ncutoffs + 1);

                for (int i = threadIdx.x; i < ncutoffs + 1; i += blockDim.x) {
                    shared_buckets[i] = 0;
                }
                for (int i = threadIdx.x; i < ncutoffs; i += blockDim.x) {
                    shared_cutoffs[i] = cutoffs[i];
                }
                __syncthreads();

                const size_t stride = static_cast<size_t>(blockDim.x) * gridDim.x;
                const size_t gid = static_cast<size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
                for (size_t i = gid; i < nvalues; i += stride) {
                    const int bucket = bucket_for_direction<Config, direction>(values[i], shared_cutoffs, ncutoffs);
                    atomicAdd(&shared_buckets[bucket], static_cast<accumulator_t>(1));
                }

                __syncthreads();

                for (int i = threadIdx.x; i < ncutoffs + 1; i += blockDim.x) {
                    const accumulator_t count = shared_buckets[i];
                    if (count != 0) atomicAdd(&buckets[i], count);
                }
            }

            template <class Config, threshold_direction direction>
            __global__ void bucket_counts_shared_cutoffs_block_buckets_kernel(
                const typename Config::value_t* __restrict__ values,
                size_t nvalues,
                const typename Config::cutoff_t* __restrict__ cutoffs,
                int ncutoffs,
                typename Config::accumulator_t* __restrict__ buckets) {
                using accumulator_t = typename Config::accumulator_t;
                using block_count_t = typename Config::block_count_t;
                using cutoff_t = typename Config::cutoff_t;

                extern __shared__ accumulator_t shared_storage[];

                auto* shared_buckets       = reinterpret_cast<block_count_t*>(shared_storage);
                const size_t bucket_bytes  = static_cast<size_t>(ncutoffs + 1) * sizeof(block_count_t);
                const size_t cutoff_offset = ((bucket_bytes + alignof(cutoff_t) - 1) / alignof(cutoff_t)) * alignof(cutoff_t);
                cutoff_t* shared_cutoffs   = reinterpret_cast<cutoff_t*>(reinterpret_cast<char*>(shared_storage) + cutoff_offset);

                for (int i = threadIdx.x; i < ncutoffs + 1; i += blockDim.x) {
                    shared_buckets[i] = 0;
                }
                for (int i = threadIdx.x; i < ncutoffs; i += blockDim.x) {
                    shared_cutoffs[i] = cutoffs[i];
                }
                __syncthreads();

                const size_t stride = static_cast<size_t>(blockDim.x) * gridDim.x;
                const size_t gid = static_cast<size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
                for (size_t i = gid; i < nvalues; i += stride) {
                    const int bucket = bucket_for_direction<Config, direction>(values[i], shared_cutoffs, ncutoffs);
                    atomicAdd(&shared_buckets[bucket], static_cast<block_count_t>(1));
                }

                __syncthreads();

                for (int i = threadIdx.x; i < ncutoffs + 1; i += blockDim.x) {
                    const block_count_t count = shared_buckets[i];
                    if (count != 0) atomicAdd(&buckets[i], static_cast<accumulator_t>(count));
                }
            }

            template <class Config, threshold_direction direction>
            __global__ void bucket_counts_shared_cutoffs_block_buckets_padded_cutoffs_kernel(
                const typename Config::value_t* __restrict__ values, size_t nvalues,
                const typename Config::cutoff_t* __restrict__ cutoffs,  int ncutoffs,
                typename Config::accumulator_t* __restrict__ buckets) {
                using accumulator_t = typename Config::accumulator_t;
                using block_count_t = typename Config::block_count_t;
                using cutoff_t = typename Config::cutoff_t;

                extern __shared__ accumulator_t shared_storage[];

                auto* shared_buckets = reinterpret_cast<block_count_t*>(shared_storage);
                const size_t bucket_bytes =
                    static_cast<size_t>(ncutoffs + 1) * sizeof(block_count_t);
                const size_t cutoff_offset =
                    ((bucket_bytes + alignof(cutoff_t) - 1) / alignof(cutoff_t)) * alignof(cutoff_t);
                cutoff_t* shared_cutoffs =
                    reinterpret_cast<cutoff_t*>(
                        reinterpret_cast<char*>(shared_storage) + cutoff_offset);

                for (int i = threadIdx.x; i < ncutoffs + 1; i += blockDim.x) {
                    shared_buckets[i] = 0;
                }
                for (int i = threadIdx.x; i < ncutoffs; i += blockDim.x) {
                    shared_cutoffs[shared_cutoff_padded_index<Config>(i)] = cutoffs[i];
                }
                __syncthreads();

                const size_t stride = static_cast<size_t>(blockDim.x) * gridDim.x;
                const size_t gid = static_cast<size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
                for (size_t i = gid; i < nvalues; i += stride) {
                    const int bucket =
                        bucket_for_padded_cutoffs_direction<Config, direction>(values[i], shared_cutoffs, ncutoffs);
                    atomicAdd(&shared_buckets[bucket], static_cast<block_count_t>(1));
                }

                __syncthreads();

                for (int i = threadIdx.x; i < ncutoffs + 1; i += blockDim.x) {
                    const block_count_t count = shared_buckets[i];
                    if (count != 0) atomicAdd(&buckets[i], static_cast<accumulator_t>(count));
                }
            }

            template <class Config, threshold_direction direction>
            __global__ void bucket_counts_global_cutoffs_shared_buckets_kernel(
                const typename Config::value_t* __restrict__ values, size_t nvalues,
                const typename Config::cutoff_t* __restrict__ cutoffs, int ncutoffs,
                typename Config::accumulator_t* __restrict__ buckets) {
                using accumulator_t = typename Config::accumulator_t;

                extern __shared__ accumulator_t shared_buckets[];

                for (int i = threadIdx.x; i < ncutoffs + 1; i += blockDim.x) {
                    shared_buckets[i] = 0;
                }

                __syncthreads();

                const size_t stride = static_cast<size_t>(blockDim.x) * gridDim.x;
                const size_t gid = static_cast<size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
                for (size_t i = gid; i < nvalues; i += stride) {
                    const int bucket = bucket_for_direction<Config, direction>(values[i], cutoffs, ncutoffs);
                    atomicAdd(&shared_buckets[bucket], static_cast<accumulator_t>(1));
                }

                __syncthreads();

                for (int i = threadIdx.x; i < ncutoffs + 1; i += blockDim.x) {
                    const accumulator_t count = shared_buckets[i];
                    if (count != 0) atomicAdd(&buckets[i], count);
                }
            }

            template <class Config, threshold_direction direction>
            __global__ void bucket_counts_uniform_cutoffs_shared_buckets_kernel(
                const typename Config::value_t* __restrict__ values,
                size_t nvalues,
                typename Config::cutoff_t first,
                typename Config::cutoff_t step,
                typename Config::cutoff_t inv_step,
                int ncutoffs,
                typename Config::accumulator_t* __restrict__ buckets) {
                using accumulator_t = typename Config::accumulator_t;

                extern __shared__ accumulator_t shared_buckets[];

                for (int i = threadIdx.x; i < ncutoffs + 1; i += blockDim.x) {
                    shared_buckets[i] = 0;
                }
                __syncthreads();

                const size_t stride = static_cast<size_t>(blockDim.x) * gridDim.x;
                const size_t gid = static_cast<size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
                for (size_t i = gid; i < nvalues; i += stride) {
                    const int bucket =
                        uniform_bucket_for_direction<Config, direction>(values[i], first, step, inv_step, ncutoffs);
                    atomicAdd(&shared_buckets[bucket], static_cast<accumulator_t>(1));
                }

                __syncthreads();

                for (int i = threadIdx.x; i < ncutoffs + 1; i += blockDim.x) {
                    const accumulator_t count = shared_buckets[i];
                    if (count != 0) atomicAdd(&buckets[i], count);
                }
            }

            template <class Config, threshold_direction direction>
            __global__ void bucket_counts_global_cutoffs_global_buckets_kernel(
                const typename Config::value_t* __restrict__ values,
                size_t nvalues,
                const typename Config::cutoff_t* __restrict__ cutoffs,
                int ncutoffs,
                typename Config::accumulator_t* __restrict__ buckets) {
                using accumulator_t = typename Config::accumulator_t;

                const size_t stride = static_cast<size_t>(blockDim.x) * gridDim.x;
                const size_t gid = static_cast<size_t>(blockIdx.x) * blockDim.x + threadIdx.x;

                for (size_t i = gid; i < nvalues; i += stride) {
                    const int bucket = bucket_for_direction<Config, direction>(values[i], cutoffs, ncutoffs);
                    atomicAdd(&buckets[bucket], static_cast<accumulator_t>(1));
                }
            }

            template <class Config, threshold_direction direction>
            __global__ void bucket_counts_uniform_cutoffs_global_buckets_kernel(
                const typename Config::value_t* __restrict__ values,
                size_t nvalues,
                typename Config::cutoff_t first,
                typename Config::cutoff_t step,
                typename Config::cutoff_t inv_step,
                int ncutoffs,
                typename Config::accumulator_t* __restrict__ buckets) {
                using accumulator_t = typename Config::accumulator_t;

                const size_t stride = static_cast<size_t>(blockDim.x) * gridDim.x;
                const size_t gid = static_cast<size_t>(blockIdx.x) * blockDim.x + threadIdx.x;

                for (size_t i = gid; i < nvalues; i += stride) {
                    const int bucket =
                        uniform_bucket_for_direction<Config, direction>(values[i], first, step, inv_step, ncutoffs);
                    atomicAdd(&buckets[bucket], static_cast<accumulator_t>(1));
                }
            }

        }
    }
}
