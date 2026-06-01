#pragma once

#include <cuda_runtime.h>
#include <cub/cub.cuh>
#include <cooperative_groups.h>

#include <limits>
#include <cmath>
#include <type_traits>
#include <climits> 

#include <propr/utils/common/constants.h>
#include <propr/utils/common/preprocessor.cuh>
#include <propr/utils/common/cuda_helpers.cuh>

#include <propr/kernels/cuda/traits/genewise.cuh>

#include <propr/data/radix.cuh>
#include <propr/internal/device/cuda/warp/reduction.cuh>

namespace propr {
    namespace dispatch {
        namespace cuda {
            template <typename Real, int THREADS_PER_BLOCK, int PAIRS_PER_THREAD>
            __global__ void genewise_connectivity_stats(
                const int*   __restrict__ partner,
                const int*   __restrict__ pair,
                const Real* __restrict__ theta,
                const Real* __restrict__ fdr,
                int          num_edges,
                Real         fdr_thresh,
                int          sort_end_bit,
                int*         __restrict__ per_gene_count,
                int*         __restrict__ per_gene_conn,
                Real*        __restrict__ per_gene_wconn,
                Real*        __restrict__ per_mean_fdr
            ){
                static_assert(THREADS_PER_BLOCK % 32 == 0, "THREADS_PER_BLOCK must be a multiple of warpSize");
                constexpr int GENES_PER_THREAD = 2 * PAIRS_PER_THREAD;
                constexpr int INVALID_KEY      = INT_MAX;

                struct __align__(16) Stats {
                    int   count = 0;   // total occurrences
                    int   conn  = 0;   // significant occurrences
                    Real  wconn = 0;   // sum(1 - theta) for significant occurrences
                    Real  fdr   = 0;
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
                        out.stat.fdr   = a.stat.fdr   + b.stat.fdr;
                        out.head       = a.head;
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

                    int   a    = INVALID_KEY, b = INVALID_KEY;
                    int   sig  = 0;
                    Real  w    = Real(0);
                    Real  fdr_ = Real(0);
                    int   has_fdr = 0;

                    if (e < num_edges) {
                        a    = partner[e];
                        b    = pair[e];
                        fdr_ = fdr[e];
                        has_fdr = !isnan(fdr_);
                        if (has_fdr && fdr_ > Real(0) && fdr_ < fdr_thresh) {
                            sig = 1; 
                            const Real theta_ = theta[e];
                            // I am not sure to be honest if it is divide or - lets wait for them
                            if (!isnan(theta_)) {
                                w = Real(1) - theta_;
                            }
                        }
                    }

                    keys[2*j + 0]       = a;
                    vals[2*j + 0].count = (a != INVALID_KEY && has_fdr)  ? 1 : 0;
                    vals[2*j + 0].conn  = (a != INVALID_KEY)  ? sig     : 0;
                    vals[2*j + 0].fdr   = (a != INVALID_KEY && has_fdr) ? fdr_ : Real(0);
                    vals[2*j + 0].wconn = (a != INVALID_KEY && sig) ? w : Real(0);

                    keys[2*j + 1]       = b;
                    vals[2*j + 1].count = (b != INVALID_KEY && has_fdr) ?   1      : 0;
                    vals[2*j + 1].conn  = (b != INVALID_KEY) ? sig      : 0;
                    vals[2*j + 1].fdr   = (b != INVALID_KEY && has_fdr) ? fdr_ : Real(0);
                    vals[2*j + 1].wconn = (b != INVALID_KEY && sig) ? w : Real(0);
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
                        atomicAdd(&per_mean_fdr[k],   items[i].stat.fdr);
                        if (items[i].stat.conn)          atomicAdd(&per_gene_conn[k],  items[i].stat.conn);
                        if (items[i].stat.wconn != Real(0)) atomicAdd(&per_gene_wconn[k], items[i].stat.wconn);
                    }
                }
            }

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
                int J   = (nbr < gene_id) ? gene_id : nbr;   // J = max
                long long idx = packed_idx_partner_major(I, J);
                return theta_edges[idx];
            }

            template <typename T, typename Config = typename propr::cuda::traits::genewise_theta_stats_config_for<T>>
            __device__ void count_radix_using_mask(
                int counts[Config::RADIX_SIZE],
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
                for (int d = 0; d < Config::RADIX_SIZE; ++d) counts[d] = 0;

                if (threadIdx.x < Config::RADIX_SIZE) smem_counts[threadIdx.x] = 0;
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
                    RadixT digit = get_bitfield<RadixT>(rv, digit_pos, Config::RADIX_BITS);
                    counts[static_cast<int>(digit)] += 1;
                }

                int lane = threadIdx.x % PROPR_WARP_SIZE;
                if (compute_sum) {
                    local_sum = propr::cuda::internal::warp::warp_reduce(local_sum, propr::ReduceSum<T>{});
                    if (lane == 0) atomicAdd(smem_sum, local_sum);
                }

                PROPR_UNROLL
                for (int d = 0; d < Config::RADIX_SIZE; ++d) {
                    int v = propr::cuda::internal::warp::warp_reduce(counts[d], propr::ReduceSum<int>{});
                    if (lane == 0) atomicAdd(&smem_counts[d], v);
                }

                __syncthreads();

                PROPR_UNROLL
                for (int d = 0; d < Config::RADIX_SIZE; ++d) counts[d] = smem_counts[d];

                __syncthreads();
            }

            template <typename T, typename Config = typename propr::cuda::traits::genewise_theta_stats_config_for<T>>
            __device__ void count_radix_contiguous(
                int counts[Config::RADIX_SIZE],
                int *smem_counts, // shared [RADIX_SIZE]
                unsigned_type_t<T> desired,
                unsigned_type_t<T> desired_mask,
                int digit_pos,
                int num_values,
                const T *__restrict__ values)
            {
                using RadixT = unsigned_type_t<T>;

                PROPR_UNROLL
                for (int d = 0; d < Config::RADIX_SIZE; ++d) counts[d] = 0;

                if (threadIdx.x < Config::RADIX_SIZE) smem_counts[threadIdx.x] = 0;
                __syncthreads();

                int iters = propr::round_up(num_values, int(blockDim.x));
                for (int pos = threadIdx.x; pos < iters; pos += blockDim.x) {
                    bool in_range = (pos < num_values);
                    T v = in_range ? values[pos] : T(0);
                    RadixT rv = propr::radix::radix_convert(v);
                    bool has_val = in_range && ((rv & desired_mask) == desired);
                    if (!has_val) continue;
                    RadixT digit = get_bitfield<RadixT>(rv, digit_pos, Config::RADIX_BITS);
                    counts[static_cast<int>(digit)] += 1;
                }

                int lane = threadIdx.x % PROPR_WARP_SIZE;

                PROPR_UNROLL
                for (int d = 0; d < Config::RADIX_SIZE; ++d) {
                    int v = propr::cuda::internal::warp::warp_reduce(counts[d], propr::ReduceSum<int>{});
                    if (lane == 0) atomicAdd(&smem_counts[d], v);
                }

                __syncthreads();
                PROPR_UNROLL
                for (int d = 0; d < Config::RADIX_SIZE; ++d) counts[d] = smem_counts[d];
                __syncthreads();
            }

            template <typename T>
            __device__ int compact_masked_bucket_to_shared(
                T *cache,
                int cache_cap,
                int *s_out_ptr, // shared scalar
                int num_edges,
                int gene_id,
                int num_genes,
                const T *__restrict__ theta_edges,
                unsigned_type_t<T> desired,
                unsigned_type_t<T> desired_mask)
            {
                using RadixT = unsigned_type_t<T>;

                if (threadIdx.x == 0) *s_out_ptr = 0;
                __syncthreads();

                int iters = propr::round_up(num_edges, int(blockDim.x));
                for (int pos = threadIdx.x; pos < iters; pos += blockDim.x) {
                    bool in_range = (pos < num_edges);
                    T fv       = in_range ? load_theta_incident<T>(gene_id, pos, num_genes, theta_edges) : T(0);
                    RadixT rv  = propr::radix::radix_convert(fv);
                    bool match = in_range && ((rv & desired_mask) == desired);

                    if (match) {
                        auto grp  = cooperative_groups::coalesced_threads();
                        int  n    = grp.size();
                        int  rank = grp.thread_rank();
                        int  base;
                        if (rank == 0) base = atomicAdd(s_out_ptr, n);
                        base = grp.shfl(base, 0);
                        int idx = base + rank;
                        if (idx < cache_cap) cache[idx] = fv;
                    }
                }

                __syncthreads();
                int out = *s_out_ptr;
                __syncthreads();
                return out;
            }

            template <typename T>
            __device__ void load_all_edges_to_cache(
                T *cache,
                int num_edges,
                int gene_id,
                int num_genes,
                const T *__restrict__ theta_edges,
                T *smem_sum) // shared scalar
            {
                if (threadIdx.x == 0) *smem_sum = T(0);
                __syncthreads();

                T local_sum = T(0);
                for (int pos = threadIdx.x; pos < num_edges; pos += blockDim.x) {
                    T fv = load_theta_incident<T>(gene_id, pos, num_genes, theta_edges);
                    local_sum += fv;
                    cache[pos] = fv;
                }

                int lane = threadIdx.x % PROPR_WARP_SIZE;
                local_sum = propr::cuda::internal::warp::warp_reduce(local_sum, propr::ReduceSum<T>{});
                if (lane == 0) atomicAdd(smem_sum, local_sum);
                __syncthreads();
            }

            template <typename T, typename Config = typename propr::cuda::traits::genewise_theta_stats_config_for<T>>
            __device__ void radix_select(
                int k_1based,       // 1-based
                int num_edges,
                int gene_id,
                int num_genes,
                int *smem_counts,   // shared [RADIX_SIZE]
                T *smem_sum,        // shared scalar
                T *cache,           // dynamic shared cache
                int cache_cap,
                int *s_compact_out, // shared scalar
                const T *__restrict__ theta_edges,
                T *median_val,
                T *theta_sum)
            {
                using RadixT = unsigned_type_t<T>;
                constexpr int RADIX_TOTAL_BITS = static_cast<int>(sizeof(RadixT) * 8);

                int counts[Config::RADIX_SIZE];
                RadixT desired = RadixT(0);
                RadixT desired_mask = RadixT(0);
                int k_to_find = k_1based;

                bool use_cache = false;
                int cache_n = 0;

                // If all edges fit in the shared-memory cache, load them once
                // and compute the sum in our single global-memory pass. Every
                // subsequent radix pass then counts from shared memory, reducing
                // total global memory traffic to exactly once ( at the start of the computation)
                if (num_edges > 0 && num_edges <= cache_cap && cache_cap > 0) {
                    load_all_edges_to_cache<T>(cache, num_edges, gene_id, num_genes, theta_edges, smem_sum);
                    cache_n = num_edges;
                    if (threadIdx.x == 0) *theta_sum = *smem_sum;
                    use_cache = true;
                }

                for (int digit_pos = RADIX_TOTAL_BITS - Config::RADIX_BITS; digit_pos >= 0; digit_pos -= Config::RADIX_BITS) {
                    if (!use_cache) {
                        bool compute_sum = (digit_pos == RADIX_TOTAL_BITS - Config::RADIX_BITS);
                        count_radix_using_mask<T>( counts, smem_counts, desired, desired_mask,
                                                   digit_pos,
                                                   num_edges, gene_id, num_genes,
                                                   theta_edges, smem_sum, /*compute_sum=*/compute_sum);
                        if (compute_sum && threadIdx.x == 0) *theta_sum = *smem_sum;
                    } else {
                        count_radix_contiguous<T>( counts, smem_counts, desired, desired_mask, digit_pos, cache_n, cache);
                    }

                    int chosen_digit = 0;
                    PROPR_UNROLL
                    for (int d = 0; d < Config::RADIX_SIZE; ++d) {
                        int count = counts[d];
                        if (count >= k_to_find) {
                            chosen_digit = d;
                            break;
                        }
                        k_to_find -= count;
                    }

                    desired = set_bitfield<RadixT>(desired, static_cast<RadixT>(chosen_digit), digit_pos, Config::RADIX_BITS);
                    desired_mask = set_bitfield<RadixT>(desired_mask, static_cast<RadixT>(Config::RADIX_MASK), digit_pos, Config::RADIX_BITS);

                    // Attempt compaction into shared memory at every global-memory pass 
                    // and not just the first. Each 4-bit digit narrows the bucket by 16x, 
                    // so after the first pass the bucket is typically N/256 which is very 
                    // well within the 2048-value cache for gene counts up to 500k
                    // (Even more for the now larger shared memory cache)
                    // The condition (digit_pos >= 2*RADIX_BITS) ensures at least 2 passes remain 
                    // as the compact itself costs one full global read of N-1 edges, 
                    // so we need the saved cache-passes to outweigh that cost
                    if (!use_cache && digit_pos >= 2 * Config::RADIX_BITS) {
                        int bucket_count = counts[chosen_digit];
                        if (bucket_count > 0 && bucket_count <= cache_cap && cache_cap > 0) {
                            cache_n = compact_masked_bucket_to_shared<T>(cache, cache_cap, s_compact_out,
                                                                         num_edges, gene_id, num_genes, theta_edges,
                                                                         desired, desired_mask);
                            use_cache = (cache_n > 0);
                        }
                    }
                }

                *median_val = radix::radix_deconvert<T>(desired);
            }

            // One block per gene: mean + lower median over N-1 incident edges
            template <typename T, typename Config = typename propr::cuda::traits::genewise_theta_stats_config_for<T>>
            __global__ 
            void genewise_theta_stats(
                const T *__restrict__ theta_edges,
                int num_genes,
                T *__restrict__ out_mean,
                T *__restrict__ out_median,
                int cache_cap_values = 0) {
                static_assert(std::is_same_v<T, float> || std::is_same_v<T, double>, "genewise_theta_stats only supports float and double");

                __shared__   int smem_counts[Config::RADIX_SIZE];
                __shared__ T smem_sum; // sum for first pass
                __shared__ T sum_all;
                __shared__ int s_compact_out;

                extern __shared__ unsigned char cache_raw[];
                T *cache = reinterpret_cast<T *>(cache_raw);

                int gene_id = (int) blockIdx.x;
                if (gene_id >= num_genes) return;

                int num_edges = num_genes - 1;
                if (num_edges <= 0) {
                    if (threadIdx.x == 0) {
                        out_mean[gene_id]   = std::numeric_limits<T>::quiet_NaN(); // i think best replace with R's nans
                        out_median[gene_id] = std::numeric_limits<T>::quiet_NaN(); // i think best replace with R's nans
                    }
                    return;
                }

                int k0 = (num_edges - 1) / 2; // lower median
                T med = T(0);

                radix_select<T>(
                    k0 + 1, num_edges, gene_id, num_genes,
                    smem_counts, &smem_sum,
                    cache, cache_cap_values, &s_compact_out,
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
