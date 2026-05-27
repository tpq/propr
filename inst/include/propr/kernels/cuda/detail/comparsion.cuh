#pragma once

#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <vector>

#include <propr/context.h>
#include <propr/data/math.cuh>

namespace propr {
    namespace detail {
        namespace cuda {

            using threshold_count_t = unsigned long long;

            enum class threshold_direction {less, greater,};

            struct threshold_cutoff {
                double value; R_xlen_t index;
            };

            struct uniform_cutoffs {
                bool valid = false;
                double first = 0.0;
                double step = 1.0;
                double inv_step = 1.0;
            };

            struct cuda_threshold_group {
                std::vector<threshold_cutoff> cutoffs;
                std::vector<double> values;
                uniform_cutoffs uniform;
                double* d_cutoffs = nullptr;
                threshold_count_t* d_buckets = nullptr;
            };

            struct cuda_threshold_counter {
                R_xlen_t ncutoffs = 0;
                cuda_threshold_group less;
                cuda_threshold_group greater;
            };

            inline uniform_cutoffs 
            describe_uniform_cutoffs(const std::vector<double>& cutoffs) {
                uniform_cutoffs uniform{};

                if (cutoffs.empty()) return uniform;
                

                uniform.first = cutoffs[0];

                if (cutoffs.size() == 1) {
                    uniform.valid = true;
                    return uniform;
                }

                uniform.step = cutoffs[1] - cutoffs[0];
                if (!(std::isfinite(uniform.step) && uniform.step > 0.0)) {
                    return uniform;
                }

                uniform.inv_step = 1.0 / uniform.step;

                for (size_t i = 2; i < cutoffs.size(); ++i) {
                    const double expected = uniform.first + uniform.step * static_cast<double>(i);
                    if (!propr::math::nearly_equal(cutoffs[i], expected)) {
                        return uniform;
                    }
                }
                uniform.valid = true;
                return uniform;
            }

            inline void 
            prepare_group(cuda_threshold_group& group, cudaStream_t stream) {
                std::stable_sort( group.cutoffs.begin(), group.cutoffs.end(),
                    [](const threshold_cutoff& a, const threshold_cutoff& b) {
                        return a.value < b.value;
                    });

                group.values.resize(group.cutoffs.size());
                for (size_t i = 0; i < group.cutoffs.size(); ++i) {
                    group.values[i] = group.cutoffs[i].value;
                }

                group.uniform = describe_uniform_cutoffs(group.values);

                if (group.cutoffs.empty()) return;

                PROPR_CUDA_CHECK(cudaMalloc(&group.d_cutoffs,group.values.size() * sizeof(double)));
                PROPR_CUDA_CHECK(cudaMemcpyAsync(
                    group.d_cutoffs,
                    group.values.data(), group.values.size() * sizeof(double),
                    cudaMemcpyHostToDevice, stream));

                PROPR_CUDA_CHECK(cudaMalloc( &group.d_buckets, (group.values.size() + 1) * sizeof(threshold_count_t)));
                PROPR_CUDA_CHECK(cudaMemsetAsync(group.d_buckets, 0, (group.values.size() + 1) * sizeof(threshold_count_t), stream));
            }

            inline void 
            destroy_group(cuda_threshold_group& group) {
                if (group.d_cutoffs != nullptr) {
                    PROPR_CUDA_CHECK(cudaFree(group.d_cutoffs));
                    group.d_cutoffs = nullptr;
                }

                if (group.d_buckets != nullptr) {
                    PROPR_CUDA_CHECK(cudaFree(group.d_buckets));
                    group.d_buckets = nullptr;
                }
            }

            // these two calls donot semantically belong here 
            PROPR_DEVICE 
            PROPR_FORCE_INLINE
            int upper_bound_cutoff(double value, const double* __restrict__ cutoffs, int ncutoffs) {
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

            PROPR_DEVICE 
            PROPR_FORCE_INLINE
            int lower_bound_cutoff(double value, const double* __restrict__ cutoffs, int ncutoffs) {
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

            template <threshold_direction direction>
            PROPR_DEVICE 
            PROPR_FORCE_INLINE
            int bucket_for_direction(double value, const double* __restrict__ cutoffs, int ncutoffs) {
                if (propr::math::is_nan(value)) return direction == threshold_direction::less ? ncutoffs : 0;
                if constexpr (direction == threshold_direction::greater) {
                    return lower_bound_cutoff(value, cutoffs, ncutoffs);
                } else {
                    return upper_bound_cutoff(value, cutoffs, ncutoffs);
                }
            }

            PROPR_DEVICE 
            PROPR_FORCE_INLINE
            int uniform_bucket_less(double value, double first, double step, double inv_step, int ncutoffs) {
                if (propr::math::is_nan(value)) return ncutoffs;
                if (ncutoffs == 1) return value < first ? 0 : 1;

                const double scaled = (value - first) * inv_step;
                int bucket = static_cast<int>(floor(scaled)) + 1;

                if (bucket < 0) return 0;
                if (bucket > ncutoffs) return ncutoffs;

                while (bucket < ncutoffs && value >= first + step * static_cast<double>(bucket)) {
                    ++bucket;
                }

                while (bucket > 0 && value < first + step * static_cast<double>(bucket - 1)) {
                    --bucket;
                }

                return bucket;
            }

            PROPR_DEVICE 
            PROPR_FORCE_INLINE
            int uniform_bucket_greater(double value, double first, double step, double inv_step, int ncutoffs) {
                if (propr::math::is_nan(value)) return 0;
                if (ncutoffs == 1) return value <= first ? 0 : 1;

                const double scaled = (value - first) * inv_step;
                int bucket = static_cast<int>(ceil(scaled));

                if (bucket < 0) return 0;
                if (bucket > ncutoffs) return ncutoffs;

                while (bucket < ncutoffs && value > first + step * static_cast<double>(bucket)) {
                    ++bucket;
                }
                while (bucket > 0 && value <= first + step * static_cast<double>(bucket - 1)) {
                    --bucket;
                }
                return bucket;
            }

            template <threshold_direction direction>
            PROPR_DEVICE 
            PROPR_FORCE_INLINE
            int uniform_bucket_for_direction(double value, double first, double step, double inv_step, int ncutoffs) {
                if constexpr (direction == threshold_direction::greater) {
                    return uniform_bucket_greater(value, first, step, inv_step, ncutoffs);
                } else {
                    return uniform_bucket_less(value, first, step, inv_step, ncutoffs);
                }
            }

            template <threshold_direction direction>
            __global__ void bucket_counts_shared_cutoffs_shared_buckets_kernel(
                const double* __restrict__ values,  size_t nvalues,
                const double* __restrict__ cutoffs, int   ncutoffs, // ncutoffs is not expected to be really big
                threshold_count_t* __restrict__ buckets) {
                extern __shared__ threshold_count_t shared_storage[];

                threshold_count_t* shared_buckets = shared_storage;
                double* shared_cutoffs = reinterpret_cast<double*>(shared_buckets + ncutoffs + 1);

                for (int i = threadIdx.x; i < ncutoffs; i += blockDim.x) {
                    shared_buckets[i] = 0;
                    shared_cutoffs[i] = cutoffs[i];
                }
                if (threadIdx.x == 0) {
                    shared_buckets[ncutoffs] = 0;
                }
                __syncthreads();

                const size_t stride = static_cast<size_t>(blockDim.x) * gridDim.x;
                const size_t gid = static_cast<size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
                for (size_t i = gid; i < nvalues; i += stride) {
                    const int bucket = bucket_for_direction<direction>(values[i], shared_cutoffs, ncutoffs);
                    atomicAdd(&shared_buckets[bucket], static_cast<threshold_count_t>(1));
                }

                __syncthreads();

                for (int i = threadIdx.x; i < ncutoffs + 1; i += blockDim.x) {
                    const threshold_count_t count = shared_buckets[i];
                    if (count != 0)  atomicAdd(&buckets[i], count);
                }
            }

            template <threshold_direction direction>
            __global__ void bucket_counts_global_cutoffs_shared_buckets_kernel(
                const double* __restrict__ values,  size_t nvalues,
                const double* __restrict__ cutoffs, int ncutoffs,
                threshold_count_t* __restrict__ buckets) {
                extern __shared__ threshold_count_t shared_buckets[];

                for (int i = threadIdx.x; i < ncutoffs + 1; i += blockDim.x) {
                    shared_buckets[i] = 0;
                }

                __syncthreads();
                size_t gid = static_cast<size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
                const size_t stride = static_cast<size_t>(blockDim.x) * gridDim.x;
                for (size_t i = gid; i < nvalues; i += stride) {
                    const int bucket = bucket_for_direction<direction>(values[i], cutoffs, ncutoffs);
                    atomicAdd(&shared_buckets[bucket], static_cast<threshold_count_t>(1));
                }

                __syncthreads();

                for (int i = threadIdx.x; i < ncutoffs + 1; i += blockDim.x) {
                    const threshold_count_t count = shared_buckets[i];
                    if (count != 0) atomicAdd(&buckets[i], count);
                }
            }

            template <threshold_direction direction>
            __global__ void bucket_counts_uniform_cutoffs_shared_buckets_kernel(
                const double* __restrict__ values, size_t nvalues,
                double first, double step, double inv_step, int ncutoffs,
                threshold_count_t* __restrict__ buckets) {
                extern __shared__ threshold_count_t shared_buckets[];

                for (int i = threadIdx.x; i < ncutoffs + 1; i += blockDim.x) {
                    shared_buckets[i] = 0;
                }
                __syncthreads();

                size_t gtid =  static_cast<size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
                const size_t stride = static_cast<size_t>(blockDim.x) * gridDim.x;
                for (size_t i = gtid; i < nvalues; i += stride) {
                    const int bucket = uniform_bucket_for_direction<direction>(values[i], first, step, inv_step, ncutoffs);
                    atomicAdd(&shared_buckets[bucket], static_cast<threshold_count_t>(1));
                }

                __syncthreads();

                for (int i = threadIdx.x; i < ncutoffs + 1; i += blockDim.x) {
                    const threshold_count_t count = shared_buckets[i];
                    if (count != 0)  atomicAdd(&buckets[i], count);
                }
            }

            template <threshold_direction direction>
            __global__ void bucket_counts_global_cutoffs_global_buckets_kernel(
                const double* __restrict__ values,  size_t nvalues,
                const double* __restrict__ cutoffs, int ncutoffs,
                threshold_count_t* __restrict__ buckets) {
                const size_t stride = static_cast<size_t>(blockDim.x) * gridDim.x;
                const size_t gid = static_cast<size_t>(blockIdx.x) * blockDim.x + threadIdx.x;

                for (size_t i = gid; i < nvalues; i += stride) {
                    const int bucket = bucket_for_direction<direction>(values[i], cutoffs, ncutoffs);
                    atomicAdd(&buckets[bucket], static_cast<threshold_count_t>(1));
                }
            }

            template <threshold_direction direction>
            __global__ void bucket_counts_uniform_cutoffs_global_buckets_kernel(
                const double* __restrict__ values,
                size_t nvalues, double first,
                double step, double inv_step,
                int ncutoffs,
                threshold_count_t* __restrict__ buckets) {
                const size_t stride = static_cast<size_t>(blockDim.x) * gridDim.x;
                const size_t gid = static_cast<size_t>(blockIdx.x) * blockDim.x + threadIdx.x;

                for (size_t i = gid; i < nvalues; i += stride) {
                    const int bucket = uniform_bucket_for_direction<direction>(values[i], first, step, inv_step, ncutoffs);
                    atomicAdd(&buckets[bucket], static_cast<threshold_count_t>(1));
                }
            }

        }
    }
}
