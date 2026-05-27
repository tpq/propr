#include <Rcpp.h>

#include <cub/device/device_scan.cuh>
#include <thrust/count.h>
#include <thrust/device_ptr.h>
#include <thrust/execution_policy.h>
#include <thrust/iterator/reverse_iterator.h>

#include <algorithm>
#include <cmath>
#include <vector>

#include <propr/context.h>

#include <propr/kernels/cuda/detail/comparsion.cuh>
#include <propr/kernels/cuda/dispatch/comparison.cuh>

#include <propr/utils/rcpp/rcpp_checks.h>
#include <propr/utils/cuda/cuda_checks.h>
#include <propr/utils/profilers/cuda_profiler.cuh>


using namespace Rcpp;
using namespace propr;


int 
dispatch::cuda::count_less_than(Rcpp::NumericVector& x,
                                double cutoff,
                                propr::propr_context context) {
    const int n = x.size();
    double* d_x = nullptr;
    PROPR_CUDA_CHECK(cudaMalloc(&d_x, n * sizeof(double)));
    PROPR_CUDA_CHECK(cudaMemcpy(
        d_x, x.begin(), n * sizeof(double),
        cudaMemcpyHostToDevice
    ));

    thrust::device_ptr<double> dev_ptr(d_x);
    auto policy = thrust::cuda::par.on(context.stream);
    int cnt = 0;
    {
        PROPR_PROFILE_CUDA("kernel", context.stream);
        cnt = thrust::count_if(
            policy,
            dev_ptr, dev_ptr + n,
            [cutoff] __device__ (double v) { return v < cutoff; }
        );
        PROPR_STREAM_SYNCHRONIZE(context);
    }

    PROPR_CUDA_CHECK(cudaFree(d_x));
    return cnt;
}

int 
dispatch::cuda::count_greater_than(Rcpp::NumericVector& x,
                                   double cutoff,
                                   propr::propr_context context) {
    const int n = x.size();
    double* d_x = nullptr;
    PROPR_CUDA_CHECK(cudaMalloc(&d_x, n * sizeof(double)));
    PROPR_CUDA_CHECK(cudaMemcpy(
        d_x, x.begin(), n * sizeof(double),
        cudaMemcpyHostToDevice
    ));

    thrust::device_ptr<double> dev_ptr(d_x);
    auto policy = thrust::cuda::par.on(context.stream);
    int cnt = 0;
    {
        PROPR_PROFILE_CUDA("kernel", context.stream);
        cnt = thrust::count_if(
            policy,
            dev_ptr, dev_ptr + n,
            [cutoff] __device__ (double v) { return v > cutoff; }
        );
        PROPR_STREAM_SYNCHRONIZE(context);
    }
    PROPR_CUDA_CHECK(cudaFree(d_x));
    return cnt;
}

int dispatch::cuda::count_less_equal_than(Rcpp::NumericVector& x,
                                          double cutoff,
                                          propr::propr_context context) {
    const int n = x.size();
    double* d_x = nullptr;
    PROPR_CUDA_CHECK(cudaMalloc(&d_x, n * sizeof(double)));
    PROPR_CUDA_CHECK(cudaMemcpy(
        d_x, x.begin(), n * sizeof(double),
        cudaMemcpyHostToDevice
    ));

    thrust::device_ptr<double> dev_ptr(d_x);
    auto policy = thrust::cuda::par.on(context.stream);
    int cnt =  0;
    {    
        PROPR_PROFILE_CUDA("kernel", context.stream);
        cnt = thrust::count_if(
            policy,
            dev_ptr, dev_ptr + n,
            [cutoff] __device__ (double v) { return v <= cutoff; }
        );
        PROPR_STREAM_SYNCHRONIZE(context);
    }

    PROPR_CUDA_CHECK(cudaFree(d_x));
    return cnt;
}

int dispatch::cuda::count_greater_equal_than(Rcpp::NumericVector& x,
                                             double cutoff,
                                             propr::propr_context context) {
    const int n = x.size();
    double* d_x = nullptr;
    PROPR_CUDA_CHECK(cudaMalloc(&d_x, n * sizeof(double)));
    PROPR_CUDA_CHECK(cudaMemcpy(
        d_x, x.begin(), n * sizeof(double),
        cudaMemcpyHostToDevice
    ));

    thrust::device_ptr<double> dev_ptr(d_x);
    auto policy = thrust::cuda::par.on(context.stream);
    int cnt = 0;
    {
        PROPR_PROFILE_CUDA("kernel", context.stream);
        cnt = thrust::count_if(
            policy,
            dev_ptr, dev_ptr + n,
            [cutoff] __device__ (double v) { return v >= cutoff; }
        );
        PROPR_STREAM_SYNCHRONIZE(context);
    
    }
    PROPR_CUDA_CHECK(cudaFree(d_x));
    return cnt;
}

template <typename InputIt, typename OutputIt>
inline void inclusive_scan_counts(
    InputIt input,
    OutputIt output,
    int nitems,
    cudaStream_t stream) {
    if (nitems <= 0) return;
    void* d_temp = nullptr;
    size_t temp_bytes = 0;
    PROPR_CUDA_CHECK(cub::DeviceScan::InclusiveSum(d_temp, temp_bytes, input, output, nitems, stream));
    PROPR_CUDA_CHECK(cudaMalloc(&d_temp, temp_bytes));
    PROPR_CUDA_CHECK(cub::DeviceScan::InclusiveSum(d_temp, temp_bytes, input, output, nitems, stream));
    PROPR_CUDA_CHECK(cudaFree(d_temp));
}

template <detail::cuda::threshold_direction direction>
inline void launch_bucket_kernel(
    detail::cuda::cuda_threshold_group& group,
    const double* d_values,
    size_t nvalues,
    propr::propr_context context) {
    if (group.cutoffs.empty()) return;
    int device = 0;
    PROPR_CUDA_CHECK(cudaGetDevice(&device));

    cudaDeviceProp prop{};
    PROPR_CUDA_CHECK(cudaGetDeviceProperties(&prop, device));

    constexpr int threads = 256;
    constexpr int blocks_per_sm = 8;
    int blocks = static_cast<int>((nvalues + threads - 1) / threads);
    blocks = std::max(1, std::min(blocks, prop.multiProcessorCount * blocks_per_sm)); // this is a heuristic for now

    const int ncutoffs = static_cast<int>(group.cutoffs.size());
    const size_t bucket_shared_bytes = static_cast<size_t>(ncutoffs + 1) * sizeof(detail::cuda::threshold_count_t);
    const size_t cutoff_shared_bytes = bucket_shared_bytes + static_cast<size_t>(ncutoffs) * sizeof(double);
    const bool buckets_fit = bucket_shared_bytes <= static_cast<size_t>(prop.sharedMemPerBlock);
    const bool cutoffs_fit = cutoff_shared_bytes <= static_cast<size_t>(prop.sharedMemPerBlock);

    if (group.uniform.valid && buckets_fit) {
        detail::cuda::bucket_counts_uniform_cutoffs_shared_buckets_kernel<direction>
            <<<blocks, threads, bucket_shared_bytes, context.stream>>>(
                d_values,
                nvalues,
                group.uniform.first,
                group.uniform.step,
                group.uniform.inv_step,
                ncutoffs,
                group.d_buckets);
    } else if (cutoffs_fit) {
        detail::cuda::bucket_counts_shared_cutoffs_shared_buckets_kernel<direction>
            <<<blocks, threads, cutoff_shared_bytes, context.stream>>>(
                d_values,
                nvalues,
                group.d_cutoffs,
                ncutoffs,
                group.d_buckets);
    } else if (buckets_fit) {
        detail::cuda::bucket_counts_global_cutoffs_shared_buckets_kernel<direction>
            <<<blocks, threads, bucket_shared_bytes, context.stream>>>(
                d_values,
                nvalues,
                group.d_cutoffs,
                ncutoffs,
                group.d_buckets);
    } else if (group.uniform.valid) {
        detail::cuda::bucket_counts_uniform_cutoffs_global_buckets_kernel<direction>
            <<<blocks, threads, 0, context.stream>>>(
                d_values,
                nvalues,
                group.uniform.first,
                group.uniform.step,
                group.uniform.inv_step,
                ncutoffs,
                group.d_buckets);
    } else {
        detail::cuda::bucket_counts_global_cutoffs_global_buckets_kernel<direction>
            <<<blocks, threads, 0, context.stream>>>(
                d_values,
                nvalues,
                group.d_cutoffs,
                ncutoffs,
                group.d_buckets);
    }

    PROPR_CUDA_CHECK(cudaGetLastError());
}

template <detail::cuda::threshold_direction direction>
inline std::vector<detail::cuda::threshold_count_t> scan_group(
    detail::cuda::cuda_threshold_group& group,
    propr::propr_context context) {
    std::vector<detail::cuda::threshold_count_t> counts(group.cutoffs.size(), 0);

    if (group.cutoffs.empty()) return counts;

    detail::cuda::threshold_count_t* d_counts = nullptr;
    PROPR_CUDA_CHECK(cudaMalloc(&d_counts, group.cutoffs.size() * sizeof(detail::cuda::threshold_count_t)));
    if constexpr (direction == detail::cuda::threshold_direction::greater) {
        auto input  = thrust::make_reverse_iterator(group.d_buckets + group.cutoffs.size() + 1);
        auto output = thrust::make_reverse_iterator(d_counts + group.cutoffs.size());
        inclusive_scan_counts(
            input,
            output,
            static_cast<int>(group.cutoffs.size()),
            context.stream);
    } else {
        inclusive_scan_counts(
            group.d_buckets,
            d_counts,
            static_cast<int>(group.cutoffs.size()),
            context.stream);
    }

    PROPR_CUDA_CHECK(cudaMemcpyAsync(
        counts.data(),
        d_counts,
        group.cutoffs.size() * sizeof(detail::cuda::threshold_count_t),
        cudaMemcpyDeviceToHost,
        context.stream));
    PROPR_STREAM_SYNCHRONIZE(context);

    PROPR_CUDA_CHECK(cudaFree(d_counts));
    return counts;
}

void* dispatch::cuda::count_values_beyond_thresholds_begin(
    Rcpp::NumericVector& cutoffs,
    bool direct,
    propr::propr_context context) {
    auto* counter = new detail::cuda::cuda_threshold_counter{};
    counter->ncutoffs = cutoffs.size();

    for (R_xlen_t i = 0; i < cutoffs.size(); ++i) {
        const double cutoff = cutoffs[i];
        if (std::isnan(cutoff)) continue;
        if (direct && cutoff >= 0.0) {
            counter->greater.cutoffs.push_back({cutoff, i});
        } else {
            counter->less.cutoffs.push_back({cutoff, i});
        }
    }

    detail::cuda::prepare_group(counter->less, context.stream);
    detail::cuda::prepare_group(counter->greater, context.stream);
    PROPR_STREAM_SYNCHRONIZE(context);

    return counter;
}

void dispatch::cuda::count_values_beyond_thresholds_accumulate(
    void* raw_counter,
    Rcpp::NumericVector& values,
    propr::propr_context context) {
    auto* counter = static_cast<detail::cuda::cuda_threshold_counter*>(raw_counter);

    if (values.size() == 0) return;

    double* d_values = nullptr;
    const size_t nvalues = static_cast<size_t>(values.size());

    PROPR_CUDA_CHECK(cudaMalloc(&d_values, nvalues * sizeof(double)));
    PROPR_CUDA_CHECK(cudaMemcpyAsync(
        d_values,
        values.begin(),
        nvalues * sizeof(double),
        cudaMemcpyHostToDevice,
        context.stream));

    {
        PROPR_PROFILE_CUDA("threshold_counting", context.stream);
        launch_bucket_kernel<detail::cuda::threshold_direction::less>( counter->less, d_values, nvalues, context);
        launch_bucket_kernel<detail::cuda::threshold_direction::greater>(counter->greater, d_values, nvalues, context);
        PROPR_STREAM_SYNCHRONIZE(context);
    }

    PROPR_CUDA_CHECK(cudaFree(d_values));
}

Rcpp::NumericVector dispatch::cuda::count_values_beyond_thresholds_end( void* raw_counter, propr::propr_context context) {
    auto* counter = static_cast<detail::cuda::cuda_threshold_counter*>(raw_counter);
    Rcpp::NumericVector out(counter->ncutoffs, NA_REAL);

    auto less_counts = scan_group<detail::cuda::threshold_direction::less>(counter->less, context);
    for (size_t i = 0; i < counter->less.cutoffs.size(); ++i) {
        out[counter->less.cutoffs[i].index] = static_cast<double>(less_counts[i]);
    }

    auto greater_counts = scan_group<detail::cuda::threshold_direction::greater>(counter->greater, context);
    for (size_t i = 0; i < counter->greater.cutoffs.size(); ++i) {
        out[counter->greater.cutoffs[i].index] = static_cast<double>(greater_counts[i]);
    }

    return out;
}

void dispatch::cuda::count_values_beyond_thresholds_destroy(void* raw_counter) {
    auto* counter = static_cast<detail::cuda::cuda_threshold_counter*>(raw_counter);
    if (counter == nullptr)  return;
    detail::cuda::destroy_group(counter->less);
    detail::cuda::destroy_group(counter->greater);
    delete counter;
}
