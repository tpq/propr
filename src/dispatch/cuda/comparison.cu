#include <Rcpp.h>

#include <thrust/count.h>
#include <thrust/device_ptr.h>
#include <thrust/execution_policy.h>
#include <cub/device/device_scan.cuh>
#include <thrust/iterator/constant_iterator.h>
#include <thrust/iterator/reverse_iterator.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/tuple.h>

#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>
#include <algorithm>
#include <type_traits>


#include <propr/context.h>
#include <propr/data/compare.cuh>
#include <propr/runtime/dispatch.hpp>

#include <propr/kernels/cuda/detail/comparsion.cuh>
#include <propr/kernels/cuda/dispatch/comparison.cuh>
#include <propr/kernels/cuda/traits/comparison.cuh>

#include <propr/utils/rcpp/rcpp_checks.h>
#include <propr/utils/rcpp/rcpp_cuda.cuh>
#include <propr/utils/cuda/cuda_checks.h>
#include <propr/utils/profilers/cuda_profiler.cuh>


using namespace Rcpp;
// using namespace propr;

template <template <typename> class Comparator, typename Real>
struct zipped_comparison {
    template <typename Tuple>
    PROPR_FORCE_INLINE
    PROPR_HOST_DEVICE
    bool operator()(const Tuple& values) const {
        return Comparator<Real>{}(thrust::get<0>(values), thrust::get<1>(values));
    }
};

template <template <typename> class Comparator>
inline int count_comparison(Rcpp::NumericVector& x,
                            double cutoff,
                            propr::propr_context context) {
    return propr::runtime::with_precision([&](auto tag) -> int {
        using Real = typename decltype(tag)::type;

        const int n = x.size();
        Real* d_x = RcppVectorToDevice<Real>(x, n);

        thrust::device_ptr<Real> dev_ptr(d_x);
        auto policy = thrust::cuda::par.on(context.stream);
        const Real typed_cutoff = static_cast<Real>(cutoff);
        auto cutoffs = thrust::make_constant_iterator(typed_cutoff);
        auto first   = thrust::make_zip_iterator(thrust::make_tuple(dev_ptr, cutoffs));
        auto last    = thrust::make_zip_iterator(thrust::make_tuple(dev_ptr + n, cutoffs + n));
        int cnt = 0;
        {
            PROPR_PROFILE_CUDA("kernel", context.stream);
            cnt = thrust::count_if(
                policy,
                first, last,
                zipped_comparison<Comparator, Real>{}
            );
            PROPR_STREAM_SYNCHRONIZE(context);
        }

        PROPR_CUDA_CHECK(cudaFree(d_x));
        return cnt;
    });
}

int 
propr::dispatch::cuda::count_less_than(Rcpp::NumericVector& x, double cutoff, propr::propr_context context) {
    return count_comparison<propr::compare::less_than>(x, cutoff, context);
}

int 
propr::dispatch::cuda::count_greater_than(Rcpp::NumericVector& x, double cutoff, propr::propr_context context) {
    return count_comparison<propr::compare::greater_than>(x, cutoff, context);
}

int propr::dispatch::cuda::count_less_equal_than(Rcpp::NumericVector& x, double cutoff, propr::propr_context context) {
    return count_comparison<propr::compare::less_equal_than>(x, cutoff, context);
}

int propr::dispatch::cuda::count_greater_equal_than(Rcpp::NumericVector& x, double cutoff, propr::propr_context context) {
    return count_comparison<propr::compare::greater_equal_than>(x, cutoff, context);
}

template <typename InputIt, typename OutputIt>
inline void inclusive_scan_counts( InputIt input, OutputIt output, int nitems, cudaStream_t stream) {
    if (nitems <= 0) return;
    void* d_temp = nullptr;
    size_t temp_bytes = 0;
    PROPR_CUDA_CHECK(cub::DeviceScan::InclusiveSum(d_temp, temp_bytes, input, output, nitems, stream));
    PROPR_CUDA_CHECK(cudaMalloc(&d_temp, temp_bytes));
    PROPR_CUDA_CHECK(cub::DeviceScan::InclusiveSum(d_temp, temp_bytes, input, output, nitems, stream));
    PROPR_CUDA_CHECK(cudaFree(d_temp));
}

template <class Config, propr::detail::cuda::threshold_direction direction>
inline void launch_bucket_kernel(
    propr::detail::cuda::cuda_threshold_group<Config>& group,
    const typename Config::value_t* d_values,
    size_t nvalues,
    propr::propr_context context) {
    if (group.cutoffs.empty()) return;

    int device = 0;
    PROPR_CUDA_CHECK(cudaGetDevice(&device));

    cudaDeviceProp prop{};
    PROPR_CUDA_CHECK(cudaGetDeviceProperties(&prop, device));

    constexpr int threads = Config::BLK_X;
    int blocks = static_cast<int>((nvalues + threads - 1) / threads);
    blocks = std::max(1, std::min(blocks, int(prop.multiProcessorCount * Config::BLOCKS_PER_SM)));

    const int ncutoffs                            = static_cast<int>(group.cutoffs.size());
    const size_t shared_capacity                  = propr::detail::cuda::dynamic_shared_capacity(prop);
    const size_t bucket_shared_bytes              = propr::detail::cuda::bucket_shared_bytes<Config>(ncutoffs);
    const size_t cutoff_shared_bytes              = propr::detail::cuda::bucket_and_cutoff_shared_bytes<Config>(ncutoffs);
    const size_t cutoff_block_shared_bytes        = propr::detail::cuda::bucket_and_cutoff_block_shared_bytes<Config>(ncutoffs);
    const size_t cutoff_block_padded_shared_bytes = propr::detail::cuda::bucket_and_cutoff_block_padded_shared_bytes<Config>(ncutoffs);

    const bool bucket_shared_fit = bucket_shared_bytes <= shared_capacity;
    const bool cutoffs_fit = cutoff_shared_bytes <= shared_capacity;
    const bool block_cutoffs_fit = cutoff_block_shared_bytes <= shared_capacity;
    const bool block_padded_cutoffs_fit = cutoff_block_padded_shared_bytes <= shared_capacity;

    const bool block_bucket_counts_safe = nvalues <= static_cast<size_t>(std::numeric_limits<typename Config::block_count_t>::max());
    const bool use_shared_cutoffs = cutoffs_fit || (block_bucket_counts_safe && (block_padded_cutoffs_fit || block_cutoffs_fit));

    if (group.uniform.valid) {
        if (bucket_shared_fit) {
            propr::detail::cuda::prepare_dynamic_shared_memory(
                propr::detail::cuda::bucket_counts_uniform_cutoffs_shared_buckets_kernel<Config, direction>,
                bucket_shared_bytes,
                prop);
            propr::detail::cuda::bucket_counts_uniform_cutoffs_shared_buckets_kernel<Config, direction>
                <<<blocks, threads, bucket_shared_bytes, context.stream>>>(
                    d_values,
                    nvalues,
                    group.uniform.first,
                    group.uniform.step,
                    group.uniform.inv_step,
                    ncutoffs,
                    group.d_buckets);
        } else {
            propr::detail::cuda::bucket_counts_uniform_cutoffs_global_buckets_kernel<Config, direction>
                <<<blocks, threads, 0, context.stream>>>(
                    d_values,
                    nvalues,
                    group.uniform.first,
                    group.uniform.step,
                    group.uniform.inv_step,
                    ncutoffs,
                    group.d_buckets);
        }
    } else if (use_shared_cutoffs) {
        if (block_bucket_counts_safe && block_padded_cutoffs_fit) {
            propr::detail::cuda::prepare_dynamic_shared_memory(
                propr::detail::cuda::bucket_counts_shared_cutoffs_block_buckets_padded_cutoffs_kernel<Config, direction>,
                cutoff_block_padded_shared_bytes,
                prop);
            
            propr::detail::cuda::bucket_counts_shared_cutoffs_block_buckets_padded_cutoffs_kernel<Config, direction>
                <<<blocks, threads, cutoff_block_padded_shared_bytes, context.stream>>>(
                    d_values,
                    nvalues,
                    group.d_cutoffs,
                    ncutoffs,
                    group.d_buckets);
        } else if (block_bucket_counts_safe && block_cutoffs_fit) {
            propr::detail::cuda::prepare_dynamic_shared_memory(
                propr::detail::cuda::bucket_counts_shared_cutoffs_block_buckets_kernel<Config, direction>,
                cutoff_block_shared_bytes,
                prop);
            
            propr::detail::cuda::bucket_counts_shared_cutoffs_block_buckets_kernel<Config, direction>
                <<<blocks, threads, cutoff_block_shared_bytes, context.stream>>>(
                    d_values,
                    nvalues,
                    group.d_cutoffs,
                    ncutoffs,
                    group.d_buckets);
        } else {
            propr::detail::cuda::prepare_dynamic_shared_memory(
                propr::detail::cuda::bucket_counts_shared_cutoffs_shared_buckets_kernel<Config, direction>,
                cutoff_shared_bytes,
                prop);
            
            propr::detail::cuda::bucket_counts_shared_cutoffs_shared_buckets_kernel<Config, direction>
                <<<blocks, threads, cutoff_shared_bytes, context.stream>>>(
                    d_values,
                    nvalues,
                    group.d_cutoffs,
                    ncutoffs,
                    group.d_buckets);
        }
    } else if (bucket_shared_fit) {
        propr::detail::cuda::prepare_dynamic_shared_memory(
            propr::detail::cuda::bucket_counts_global_cutoffs_shared_buckets_kernel<Config, direction>,
            bucket_shared_bytes,
            prop);
        
        propr::detail::cuda::bucket_counts_global_cutoffs_shared_buckets_kernel<Config, direction>
            <<<blocks, threads, bucket_shared_bytes, context.stream>>>(
                d_values,
                nvalues,
                group.d_cutoffs,
                ncutoffs,
                group.d_buckets);
    } else {
        propr::detail::cuda::bucket_counts_global_cutoffs_global_buckets_kernel<Config, direction>
            <<<blocks, threads, 0, context.stream>>>(
                d_values,
                nvalues,
                group.d_cutoffs,
                ncutoffs,
                group.d_buckets);
    }

    PROPR_CUDA_CHECK(cudaGetLastError());
}

template <class Config, propr::detail::cuda::threshold_direction direction>
inline std::vector<typename Config::accumulator_t> scan_group(
    propr::detail::cuda::cuda_threshold_group<Config>& group,
    propr::propr_context context) {
    using accumulator_t = typename Config::accumulator_t;

    std::vector<accumulator_t> counts(group.cutoffs.size(), 0);

    if (group.cutoffs.empty()) return counts;

    accumulator_t* d_counts = nullptr;
    PROPR_CUDA_CHECK(cudaMalloc(&d_counts, group.cutoffs.size() * sizeof(accumulator_t)));
    if constexpr (direction == propr::detail::cuda::threshold_direction::greater) {
        auto input  = thrust::make_reverse_iterator(group.d_buckets + group.cutoffs.size() + 1);
        auto output = thrust::make_reverse_iterator(d_counts + group.cutoffs.size());
        inclusive_scan_counts( input, output, static_cast<int>(group.cutoffs.size()), context.stream);
    } else {
        inclusive_scan_counts( group.d_buckets, d_counts, static_cast<int>(group.cutoffs.size()), context.stream);
    }

    PROPR_CUDA_CHECK(cudaMemcpyAsync(counts.data(), d_counts, group.cutoffs.size() * sizeof(accumulator_t), 
                                     cudaMemcpyDeviceToHost, context.stream));
    PROPR_STREAM_SYNCHRONIZE(context);
    PROPR_CUDA_CHECK(cudaFree(d_counts));
    return counts;
}

struct cuda_threshold_counter_handle {
    propr::runtime::Precision precision;
    void* counter = nullptr;
};

template <typename Real>
using threshold_config_for = propr::cuda::traits::count_values_beyond_thresholds_config_for<Real>;

template <typename Real>
void destroy_threshold_counter_typed(void* raw_counter) {
    using Config = threshold_config_for<Real>;
    auto* counter = static_cast<propr::detail::cuda::cuda_threshold_counter<Config>*>(raw_counter);
    if (counter == nullptr) return;
    propr::detail::cuda::destroy_group<Config>(counter->less);
    propr::detail::cuda::destroy_group<Config>(counter->greater);
    delete counter;
}

template <typename Real>
void* make_threshold_counter_typed(Rcpp::NumericVector& cutoffs,
                                   bool direct,
                                   propr::propr_context context) {
    using Config = threshold_config_for<Real>;
    using cutoff_t = typename Config::cutoff_t;

    auto* counter = new propr::detail::cuda::cuda_threshold_counter<Config>{};

    try {
        counter->ncutoffs = cutoffs.size();
        for (R_xlen_t i = 0; i < cutoffs.size(); ++i) {
            const double cutoff = cutoffs[i];
            if (std::isnan(cutoff)) continue;
            const cutoff_t typed_cutoff = static_cast<cutoff_t>(cutoff);
            if (direct && typed_cutoff >= static_cast<cutoff_t>(0)) {
                counter->greater.cutoffs.push_back({typed_cutoff, i});
            } else {
                counter->less.cutoffs.push_back({typed_cutoff, i});
            }
        }
        propr::detail::cuda::prepare_group<Config>(counter->less, context.stream);
        propr::detail::cuda::prepare_group<Config>(counter->greater, context.stream);
        PROPR_STREAM_SYNCHRONIZE(context);
        return counter;
    } catch (...) {
        destroy_threshold_counter_typed<Real>(counter);
        throw;
    }
}

template <typename Real>
void accumulate_threshold_counter_typed(
    void* raw_counter,
    Rcpp::NumericVector& values,
    propr::propr_context context) {
    using Config = threshold_config_for<Real>;
    using value_t = typename Config::value_t;

    if (values.size() == 0) return;

    auto* counter =static_cast<propr::detail::cuda::cuda_threshold_counter<Config>*>(raw_counter);
    const size_t nvalues = static_cast<size_t>(values.size());
    value_t* d_values = RcppVectorToDevice<value_t>(values, nvalues);

    {
        PROPR_PROFILE_CUDA("threshold_counting", context.stream);
        launch_bucket_kernel<Config, propr::detail::cuda::threshold_direction::less>( counter->less, d_values, nvalues, context);
        launch_bucket_kernel<Config, propr::detail::cuda::threshold_direction::greater>(counter->greater, d_values, nvalues, context);
        PROPR_STREAM_SYNCHRONIZE(context);
    }

    PROPR_CUDA_CHECK(cudaFree(d_values));
}

template <typename Real>
Rcpp::NumericVector finish_threshold_counter_typed(
    void* raw_counter,
    propr::propr_context context) {
    using Config = threshold_config_for<Real>;
    auto* counter = static_cast<propr::detail::cuda::cuda_threshold_counter<Config>*>(raw_counter);
    Rcpp::NumericVector out(counter->ncutoffs, NA_REAL);
    auto less_counts = scan_group<Config, propr::detail::cuda::threshold_direction::less>(counter->less, context);
    for (size_t i = 0; i < counter->less.cutoffs.size(); ++i) {
        out[counter->less.cutoffs[i].index] = static_cast<double>(less_counts[i]);
    }
    auto greater_counts = scan_group<Config, propr::detail::cuda::threshold_direction::greater>( counter->greater, context);
    for (size_t i = 0; i < counter->greater.cutoffs.size(); ++i) {
        out[counter->greater.cutoffs[i].index] = static_cast<double>(greater_counts[i]);
    }
    return out;
}

template <typename F>
decltype(auto) with_counter_precision(propr::runtime::Precision precision, F&& f) {
    switch (precision) {
        case propr::runtime::Precision::Float32: return f(propr::runtime::precision_tag<float>{});
        case propr::runtime::Precision::Float64: return f(propr::runtime::precision_tag<double>{});
    }
    throw std::logic_error("Internal error: unknown propr.precision.");
}

void* propr::dispatch::cuda::count_values_beyond_thresholds_begin(
    Rcpp::NumericVector& cutoffs,
    bool direct,
    propr::propr_context context) {
    const auto precision = propr::runtime::resolve_precision();
    auto* handle = new cuda_threshold_counter_handle{precision, nullptr};

    try {
        // i kinda like this syntax, even though its inconsistent
        handle->counter =
            with_counter_precision(precision, [&](auto tag) -> void* {
                using Real = typename decltype(tag)::type;
                return make_threshold_counter_typed<Real>(cutoffs, direct, context);
            });
        return handle;
    } catch (...) {
        delete handle;
        throw;
    }
}

void propr::dispatch::cuda::count_values_beyond_thresholds_accumulate(
    void* raw_counter,
    Rcpp::NumericVector& values,
    propr::propr_context context) {
    auto* handle = static_cast<cuda_threshold_counter_handle*>(raw_counter);
    with_counter_precision(handle->precision, [&](auto tag) {
        using Real = typename decltype(tag)::type;
        accumulate_threshold_counter_typed<Real>(handle->counter, values, context);
    });
}

Rcpp::NumericVector propr::dispatch::cuda::count_values_beyond_thresholds_end(
    void* raw_counter,
    propr::propr_context context) {
    auto* handle = static_cast<cuda_threshold_counter_handle*>(raw_counter);
    return with_counter_precision(
        handle->precision,
        [&](auto tag) -> Rcpp::NumericVector {
            using Real = typename decltype(tag)::type;
            return finish_threshold_counter_typed<Real>(handle->counter, context);
        });
}

void propr::dispatch::cuda::count_values_beyond_thresholds_destroy(void* raw_counter) {
    auto* handle = static_cast<cuda_threshold_counter_handle*>(raw_counter);
    if (handle == nullptr) return;
    try {
        with_counter_precision(handle->precision, [&](auto tag) {
            using Real = typename decltype(tag)::type;
            destroy_threshold_counter_typed<Real>(handle->counter);
        });
    } catch (...) {
        delete handle;
        throw;
    }
    delete handle;
}
