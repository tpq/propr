#pragma once

#include <cuda_runtime.h>
#include <cub/cub.cuh>
#include <limits>
#include <climits> 

#include <propr/utils/common/constants.h>
#include <propr/utils/common/preprocessor.cuh>
#include <propr/utils/common/cuda_helpers.cuh>


#include <propr/internal/device/cuda/warp/reduction.cuh>
#include <propr/data/radix.cuh>

namespace propr {
    namespace dispatch {
        namespace cuda {

            template <int THREADS_PER_BLOCK, int PAIRS_PER_THREAD>
            __global__ void genewise_connectivity_stats(
                const int*   __restrict__ partner,
                const int*   __restrict__ pair,
                const float* __restrict__ theta,
                const float* __restrict__ fdr,
                int          num_edges,
                float        fdr_thresh,
                int          sort_end_bit,
                int*         __restrict__ per_gene_count,
                int*         __restrict__ per_gene_conn,
                float*       __restrict__ per_gene_wconn)
            {
                static_assert(THREADS_PER_BLOCK % 32 == 0, "THREADS_PER_BLOCK must be a multiple of warpSize");
                constexpr int GENES_PER_THREAD = 2 * PAIRS_PER_THREAD;
                constexpr int INVALID_KEY      = INT_MAX;

                struct __align__(16) Stats {
                    int   count = 0;   // total occurrences
                    int   conn  = 0;   // significant occurrences
                    float wconn = 0;   // sum(1-theta) for significant occurrences
                    int   _pad  = 0;
                };

                struct SegItem {
                    Stats stat;
                    int   head;   // 1 if the position is the start of a run
                };

                struct SegOp {
                    __device__ __forceinline__ SegItem operator()(SegItem const& a, SegItem const& b) const {
                        if (b.head) return b;
                        SegItem out;
                        out.stat.count = a.stat.count + b.stat.count;
                        out.stat.conn  = a.stat.conn  + b.stat.conn;
                        out.stat.wconn = a.stat.wconn + b.stat.wconn;
                        out.head    = a.head;
                        return out;
                    }
                };

                using BlockSortT = cub::BlockRadixSort<int, THREADS_PER_BLOCK, GENES_PER_THREAD, Stats>;
                using BlockDiscT = cub::BlockDiscontinuity<int, THREADS_PER_BLOCK, GENES_PER_THREAD>;
                using BlockScanT = cub::BlockScan<SegItem, THREADS_PER_BLOCK, cub::BLOCK_SCAN_WARP_SCANS, GENES_PER_THREAD>;

                __shared__ union {
                    typename BlockSortT::TempStorage sort;
                    typename BlockDiscT::TempStorage disc;
                    typename BlockScanT::TempStorage scan;
                } temp;

                const int tid = (int)threadIdx.x;
                const int blk_pair_offset = (int)blockIdx.x * THREADS_PER_BLOCK * PAIRS_PER_THREAD;

                int   keys[GENES_PER_THREAD];
                Stats vals[GENES_PER_THREAD];

                PROPR_UNROLL
                for (int j = 0; j < PAIRS_PER_THREAD; ++j) {
                    const int e = blk_pair_offset + j * THREADS_PER_BLOCK + tid;

                    int   a   = INVALID_KEY, b = INVALID_KEY;
                    int   sig = 0;
                    float w   = 0.0f;

                    if (e < num_edges) {
                        a = partner[e];
                        b = pair[e];
                        if (fdr[e] < fdr_thresh) {
                            sig = 1; w   = 1.0f - theta[e];  // mark the pair as sig
                        }
                    }

                    keys[2*j + 0]       = a;
                    vals[2*j + 0].count = (a != INVALID_KEY) ?   1 : 0;
                    vals[2*j + 0].conn  = (a != INVALID_KEY) ? sig : 0;
                    vals[2*j + 0].wconn = (a != INVALID_KEY && sig) ? w : 0.0f;

                    keys[2*j + 1]       = b;
                    vals[2*j + 1].count = (b != INVALID_KEY) ?   1 : 0;
                    vals[2*j + 1].conn  = (b != INVALID_KEY) ? sig : 0;
                    vals[2*j + 1].wconn = (b != INVALID_KEY && sig) ? w : 0.0f;
                }

                BlockSortT(temp.sort).Sort(keys, vals, 0, sort_end_bit);
                __syncthreads();

                int head_flags[GENES_PER_THREAD];
                BlockDiscT(temp.disc).FlagHeads(head_flags, keys, cub::Inequality());
                __syncthreads();

                SegItem items[GENES_PER_THREAD];
                PROPR_UNROLL
                for (int i = 0; i < GENES_PER_THREAD; ++i) {
                    items[i].stat = vals[i];
                    items[i].head = head_flags[i];
                }

                BlockScanT(temp.scan).InclusiveScan(items, items, SegOp{});
                __syncthreads();

                int tail_flags[GENES_PER_THREAD];
                BlockDiscT(temp.disc).FlagTails(tail_flags, keys, cub::Inequality());
                __syncthreads();

                PROPR_UNROLL
                for (int i = 0; i < GENES_PER_THREAD; ++i) {
                    const int k = keys[i];
                    if (k == INVALID_KEY) continue;
                    if (tail_flags[i]) {
                        atomicAdd(&per_gene_count[k], items[i].stat.count);
                        if (items[i].stat.conn)          atomicAdd(&per_gene_conn[k],  items[i].stat.conn);
                        if (items[i].stat.wconn != 0.0f) atomicAdd(&per_gene_wconn[k], items[i].stat.wconn);
                    }
                }
            }

            constexpr int RADIX_BITS = 4;
            constexpr int RADIX_SIZE = 1 << RADIX_BITS; // 16
            constexpr int RADIX_MASK = RADIX_SIZE - 1;  // 15

            // Packed order index for edge (I<J), 0-based genes, with partner-major storage.
            // idx = J*(J-1)/2 + I  (because partner=J+1 has J entries: pairs 0..J-1)
            static __host__ __device__ __forceinline__ long long packed_idx_partner_major(int I, int J) {
                return (long long)J * (J - 1LL) / 2LL + (long long)I;
            }

            // this is the part that allows us to hack
            // For gene g, virtual slice of length N-1: neighbors != g
            // pos -> neighbor = pos if pos < g else pos+1
            template <typename T>
            static __device__ __forceinline__ T load_theta_incident(
                int gene_id, int pos, int num_genes, const T* __restrict__ theta_edges) {
                int nbr = (pos < gene_id) ? pos     : (pos + 1);
                int I   = (nbr < gene_id) ? nbr     : gene_id;
                int J   = (nbr < gene_id) ? gene_id : nbr;   // J=max
                long long idx = packed_idx_partner_major(I, J);
                return theta_edges[idx];
            }

            template <typename T>
            __device__ void count_radix_using_mask(
                int counts[RADIX_SIZE],
                int *smem_counts, // shared [RADIX_SIZE]
                unsigned_type_t<T> desired,
                unsigned_type_t<T> desired_mask,
                int digit_pos,
                int num_edges,
                int gene_id,
                int num_genes,
                const T *__restrict__ theta_edges,
                T *smem_sum, // shared scalar (valid only if compute_sum=true)
                bool compute_sum)
            {
                using RadixT = unsigned_type_t<T>;

                PROPR_UNROLL
                for (int d = 0; d < RADIX_SIZE; ++d) counts[d] = 0;

                if (threadIdx.x < RADIX_SIZE) smem_counts[threadIdx.x] = 0;
                if (compute_sum && threadIdx.x == 0) *smem_sum = T(0);
                __syncthreads();

                T local_sum = T(0);
                int iters = propr::round_up(num_edges, int(blockDim.x));
                
                for (int pos = threadIdx.x; pos < iters; pos += blockDim.x) {
                    bool in_range = (pos < num_edges);
                    T fv = in_range ? load_theta_incident<T>(gene_id, pos, num_genes, theta_edges) : T(0);

                    if (compute_sum && in_range) local_sum += fv;

                    RadixT rv  = propr::radix::radix_convert(fv);
                    bool has_val = in_range && ((rv & desired_mask) == desired);
                    if (!has_val) continue;
                    RadixT digit = get_bitfield<RadixT>(rv, digit_pos, RADIX_BITS);
                    counts[static_cast<int>(digit)] += 1;
                }

                int lane = threadIdx.x % PROPR_WARP_SIZE;

                if (compute_sum) {
                    local_sum = propr::cuda::internal::warp::warp_reduce(local_sum, propr::ReduceSum<T>{});
                    if (lane == 0) atomicAdd(smem_sum, local_sum);
                }

                PROPR_UNROLL
                for (int d = 0; d < RADIX_SIZE; ++d) {
                    int v = propr::cuda::internal::warp::warp_reduce(counts[d], propr::ReduceSum<int>{});
                    if (lane == 0) atomicAdd(&smem_counts[d], v);
                }

                __syncthreads();

                PROPR_UNROLL
                for (int d = 0; d < RADIX_SIZE; ++d) counts[d] = smem_counts[d];

                __syncthreads();
            }

            template <typename T>
            __device__ T find_pattern(
                T *smem_flag_val, // shared [2] : [0]=flag, [1]=value
                int num_edges,
                int gene_id,
                int num_genes,
                const T *__restrict__ theta_edges,
                unsigned_type_t<T> desired,
                unsigned_type_t<T> desired_mask) {
                if (threadIdx.x < 2) smem_flag_val[threadIdx.x] = T(0);
                __syncthreads();

                int iters = propr::round_up(num_edges, int(blockDim.x));
                for (int pos = threadIdx.x; pos < iters; pos += blockDim.x) {
                    bool in_range = (pos < num_edges);
                    T v = in_range ? load_theta_incident<T>(gene_id, pos, num_genes, theta_edges) : T(0);
                    if (in_range && ((radix::radix_convert(v) & desired_mask) == desired)) {
                            smem_flag_val[0] = T(1);
                            smem_flag_val[1] = v;
                    }
                    __syncthreads();
                    T found = smem_flag_val[0];
                    T val   = smem_flag_val[1];
                    __syncthreads();
                    if (found != T(0)) return val;
                }

                return std::numeric_limits<T>::quiet_NaN();
            }

            template<typename T, bool EARLY_EXIT=false>
            __device__ void radix_select(
                int k_1based, // 1-based
                int num_edges,
                int gene_id,
                int num_genes,
                int   *smem_counts,   // shared [RADIX_SIZE]
                T *smem_flag_val, // shared [2]
                T *smem_sum,      // shared scalar
                const T *__restrict__ theta_edges,
                T *median_val,
                T *theta_sum) {
                
                using RadixT = unsigned_type_t<T>;
                constexpr int RADIX_TOTAL_BITS = static_cast<int>(sizeof(RadixT) * 8);

                int counts[RADIX_SIZE];
                RadixT desired = RadixT(0);
                RadixT desired_mask = RadixT(0);
                int k_to_find = k_1based;

                bool first = true;
               
                for (int digit_pos = RADIX_TOTAL_BITS - RADIX_BITS; digit_pos >= 0; digit_pos -= RADIX_BITS) {
                    count_radix_using_mask<T>(
                        counts, smem_counts, desired, desired_mask,
                        digit_pos, num_edges, gene_id, num_genes, theta_edges,
                        smem_sum, /*compute_sum=*/first);

                    if (first) {
                        if (threadIdx.x == 0) *theta_sum = *smem_sum; // sum over ALL theta values
                        first = false;
                    }

                    if constexpr (EARLY_EXIT) {
                        auto found_unique = [&](int digit, int count) -> bool {
                            if (count == 1 && k_to_find == 1) {
                                desired      = set_bitfield<RadixT>(desired, static_cast<RadixT>(digit), digit_pos, RADIX_BITS);
                                desired_mask = set_bitfield<RadixT>(desired_mask, static_cast<RadixT>(RADIX_MASK), digit_pos, RADIX_BITS);
                                *median_val = find_pattern<T>(smem_flag_val, num_edges, gene_id, num_genes, theta_edges, desired, desired_mask);
                                return true;
                            }
                            return false;
                        };
                        auto found_non_unique = [&](int digit, int count) -> bool {
                            if (count >= k_to_find) {
                                desired      = set_bitfield<RadixT>(desired, static_cast<RadixT>(digit), digit_pos, RADIX_BITS);
                                desired_mask = set_bitfield<RadixT>(desired_mask, static_cast<RadixT>(RADIX_MASK), digit_pos, RADIX_BITS);
                                return true;
                            }
                            k_to_find -= count;
                            return false;
                        };
                        // k-th smallest
                        PROPR_UNROLL
                        for (int d = 0; d < RADIX_SIZE; ++d) {
                            int c = counts[d];
                            if (found_unique(d, c)) return;
                            if (found_non_unique(d, c)) break;
                        }
                    } else {
                        PROPR_UNROLL
                        for (int d = 0; d < RADIX_SIZE; ++d) {
                            int c = counts[d];
                            if (c >= k_to_find) {
                                desired      = set_bitfield<RadixT>(desired, static_cast<RadixT>(d), digit_pos, RADIX_BITS);
                                desired_mask = set_bitfield<RadixT>(desired_mask, static_cast<RadixT>(RADIX_MASK), digit_pos, RADIX_BITS);
                                break;
                            }
                            k_to_find -= c;
                        }
                    }

                }
                *median_val = radix::radix_deconvert<T>(desired);
            }

            // One block per gene: mean + lower median over N-1 incident edges
            template <typename T>
            __global__ 
            void genewise_theta_stats(
                const T *__restrict__ theta_edges,
                int num_genes,
                T *__restrict__ out_mean,
                T *__restrict__ out_median) {
                static_assert(std::is_same_v<T, float> || std::is_same_v<T, double>, "genewise_theta_stats only supports float and double");

                __shared__   int smem_counts[RADIX_SIZE];
                __shared__ T smem_flag_val[2]; // keep find
                __shared__ T smem_sum; // sum for first pass
                __shared__ T sum_all;

                int gene_id = (int) blockIdx.x;
                if (gene_id >= num_genes) return;

                int num_edges = num_genes - 1;
                if (num_edges <= 0) {
                    if (threadIdx.x == 0) {
                        out_mean[gene_id] = std::numeric_limits<T>::quiet_NaN();
                        out_median[gene_id] = std::numeric_limits<T>::quiet_NaN();
                    }
                    return;
                }

                int k0 = (num_edges - 1) / 2; // lower median
                T med = T(0);

                radix_select<T>(
                    k0 + 1, num_edges, gene_id, num_genes,
                    smem_counts, smem_flag_val, &smem_sum,
                    theta_edges,
                    &med,
                    &sum_all);

                if (threadIdx.x == 0) {
                    out_median[gene_id] = med;
                    out_mean[gene_id] = sum_all / static_cast<T>(num_edges);
                }
            }

        } // namespace cuda 
    } // namespace dispatch 
} // namespace propr
