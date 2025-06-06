#include <cfloat>
#include <math.h>
#include <stdio.h>
#include "../../../include/kernels/cuda/detail/lrv.cuh"

#define LRV_BLOCK_SIZE 16

__global__
void
propr::detail::cuda::lrv_basic(float* __restrict__ d_Y,
                               float* __restrict__ d_variances,
                               int nb_samples, int nb_genes) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i < nb_genes && j < i) {
        float2 accum = make_float2(0.0f, 0.0f);
        int k = 0;
        #pragma unroll
        for (; k < (nb_samples/4)*4; k += 4) {
            float4 y_i = *reinterpret_cast<float4*>(&d_Y[k + i * nb_samples]);
            float4 y_j = *reinterpret_cast<float4*>(&d_Y[k + j * nb_samples]);

            #pragma unroll
            for (int m = 0; m < 4; ++m) {
                float ratio = __fdividef((&y_i.x)[m], (&y_j.x)[m]);
                float log_val = __logf(ratio);

                accum.x = __fmaf_rn(1.0f, log_val,    accum.x);
                accum.y = __fmaf_rn(log_val, log_val, accum.y); 
            }
        }

        for (; k < nb_samples; ++k) {
            const float yi = d_Y[k + i * nb_samples];
            const float yj = d_Y[k + j * nb_samples];

            const float ratio   = __fdividef(yi, yj);
            const float log_val = __logf(ratio);

            accum.x = __fmaf_rn(1.0f, log_val,    accum.x);
            accum.y = __fmaf_rn(log_val, log_val, accum.y);
        }

        float inv_n    = __frcp_rn(static_cast<float>(nb_samples));
        float mean     = accum.x * inv_n;
        float variance = (accum.y - __fmul_rn(nb_samples, __fmul_rn(mean, mean))) * __frcp_rn(static_cast<float>(nb_samples - 1));

        int pair_index = (i * (i - 1)) / 2 + j;
        d_variances[pair_index] = variance;
    }
}


__global__
void
propr::detail::cuda::lrv_weighted(float* __restrict__ d_Y,
                        float* __restrict__ d_W,
                        float* __restrict__ d_variances, 
                        int nb_samples, int nb_genes) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i < nb_genes && j < i) {
        // accum.x = w_sum, accum.y = w_sum2, accum.z = mean, accum.w = numerator (S)
        float4 accum = make_float4(0.0f, 0.0f, 0.0f, 0.0f);
        int k = 0;

        #pragma unroll
        for (; k <= nb_samples - 4; k += 4) {
            float4 y_i = *reinterpret_cast<float4*>(&d_Y[k + i * nb_samples]);
            float4 y_j = *reinterpret_cast<float4*>(&d_Y[k + j * nb_samples]);
            float4 w_i = *reinterpret_cast<float4*>(&d_W[k + i * nb_samples]);
            float4 w_j = *reinterpret_cast<float4*>(&d_W[k + j * nb_samples]);

            #pragma unroll
            for (int m = 0; m < 4; ++m) {
                float ratio = __fdividef((&y_i.x)[m], (&y_j.x)[m]);
                float log_val = __logf(ratio);
                float mean_old = accum.z;
                float w = (&w_i.x)[m] * (&w_j.x)[m];

                accum.x = __fmaf_rn(1.0f, w, accum.x);
                accum.y = __fmaf_rn(w, w, accum.y);

                float delta = log_val - mean_old;
                float w_ratio = __fdividef(w, accum.x);
                accum.z = __fmaf_rn(w_ratio, delta, mean_old);

                float delta_new = log_val - accum.z;
                accum.w = __fmaf_rn(w, delta * delta_new, accum.w);
            }
        }

        for (; k < nb_samples; ++k) {
            float y_ik = d_Y[k + i * nb_samples];
            float y_jk = d_Y[k + j * nb_samples];
            float w_ik = d_W[k + i * nb_samples];
            float w_jk = d_W[k + j * nb_samples];
            float w_k = w_ik * w_jk;

            float ratio = __fdividef(y_ik, y_jk);
            float log_val = __logf(ratio);
            float mean_old = accum.z;

            // Update w_sum and w_sum2
            accum.x += w_k;
            accum.y += w_k * w_k;

            // Compute new mean
            float delta = log_val - mean_old;
            float w_ratio = w_k / accum.x;
            accum.z += w_ratio * delta;

            // Update numerator (S)
            float delta_new = log_val - accum.z;
            accum.w += w_k * delta * delta_new;
        }

        float S_total = accum.x;
        float Q_total = accum.y;
        float numerator = accum.w;

        float denominator = S_total - (Q_total / S_total);
        float lrv = numerator / denominator;

        int pair_index = (i * (i - 1)) / 2 + j;
        d_variances[pair_index] = lrv;
    }
}


__global__
void
propr::detail::cuda::lrv_alpha(float* __restrict__ d_Y,
                     float* __restrict__ d_Yfull,
                     float a,
                     float* __restrict__ d_variances,
                     int nb_samples, int nb_genes) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i < nb_genes && j < i) {                    
        float4 mu_m = make_float4(0.0f, 0.0f, 0.0f, 0.0f);
        float2 acc  = make_float2(0.0f, 0.0f);
        float C = 0.0f;
        int n = 0, k = 0;

        #pragma unroll
        for (; k <= nb_samples - 4; k += 4) {
            float4 y_i     = *reinterpret_cast<float4*>(&d_Y[k + i * nb_samples]);
            float4 y_j     = *reinterpret_cast<float4*>(&d_Y[k + j * nb_samples]);
            float4 yfull_i = *reinterpret_cast<float4*>(&d_Yfull[k + i * nb_samples]);
            float4 yfull_j = *reinterpret_cast<float4*>(&d_Yfull[k + j * nb_samples]);

            #pragma unroll
            for (int m = 0; m < 4; ++m) {
                n++;
                float inv_n = __frcp_rn(static_cast<float>(n));
                
                float X_i      = __powf(reinterpret_cast<float*>(&y_i)[m], a);
                float X_j      = __powf(reinterpret_cast<float*>(&y_j)[m], a);
                float X_full_i = __powf(reinterpret_cast<float*>(&yfull_i)[m], a);
                float X_full_j = __powf(reinterpret_cast<float*>(&yfull_j)[m], a);

                float prev_mu_i = mu_m.x;
                float dx_i = __fsub_rn(X_i, prev_mu_i);
                mu_m.x = __fmaf_rn(dx_i, inv_n, prev_mu_i);

                float prev_mu_j = mu_m.y;
                float dx_j = __fsub_rn(X_j, prev_mu_j);
                mu_m.y = __fmaf_rn(dx_j, inv_n, prev_mu_j);

                float dxj_muj = __fsub_rn(X_j, mu_m.y);
                C = __fmaf_rn(dx_i, dxj_muj, C);

                acc.x = __fmaf_rn(X_i, X_i, acc.x);
                acc.y = __fmaf_rn(X_j, X_j, acc.y);

                float dx_full_i = __fsub_rn(X_full_i, mu_m.z);
                mu_m.z = __fmaf_rn(dx_full_i, inv_n, mu_m.z);

                float dx_full_j = __fsub_rn(X_full_j, mu_m.w);
                mu_m.w = __fmaf_rn(dx_full_j, inv_n, mu_m.w);
            }
        }

        for (; k < nb_samples; ++k) {
            n++;
            float inv_n = __frcp_rn(static_cast<float>(n));
            float X_i      = __powf(d_Y[k + i * nb_samples], a);
            float X_j      = __powf(d_Y[k + j * nb_samples], a);
            float X_full_i = __powf(d_Yfull[k + i * nb_samples], a);
            float X_full_j = __powf(d_Yfull[k + j * nb_samples], a);

            float prev_mu_i = mu_m.x;
            float dx_i = __fsub_rn(X_i, prev_mu_i);
            mu_m.x = __fmaf_rn(dx_i, inv_n, prev_mu_i);

            float prev_mu_j = mu_m.y;
            float dx_j = __fsub_rn(X_j, prev_mu_j);
            mu_m.y = __fmaf_rn(dx_j, inv_n, prev_mu_j);

            float dxj_muj = __fsub_rn(X_j, mu_m.y);
            C = __fmaf_rn(dx_i, dxj_muj, C);
            
            acc.x = __fmaf_rn(X_i, X_i, acc.x);
            acc.y = __fmaf_rn(X_j, X_j, acc.y);

            float dx_full_i = __fsub_rn(X_full_i, mu_m.z);
            mu_m.z = __fmaf_rn(dx_full_i, inv_n, mu_m.z);

            float dx_full_j = __fsub_rn(X_full_j, mu_m.w);
            mu_m.w = __fmaf_rn(dx_full_j, inv_n, mu_m.w);
        }

        float n_mui_sq = __fmul_rn(n, __fmul_rn(mu_m.x, mu_m.x));
        float n_muj_sq = __fmul_rn(n, __fmul_rn(mu_m.y, mu_m.y));
        
        float sum_sq_i = __fsub_rn(acc.x, n_mui_sq);
        float sum_sq_j = __fsub_rn(acc.y, n_muj_sq);

        float a_i = (mu_m.z != 0.0f) ? __frcp_rn(mu_m.z) : 0.0f;
        float a_j = (mu_m.w != 0.0f) ? __frcp_rn(mu_m.w) : 0.0f;

        float ai_sq = __fmul_rn(a_i, a_i);
        float aj_sq = __fmul_rn(a_j, a_j);
        float aij = __fmul_rn(a_i, a_j);

        float term1 = __fmul_rn(sum_sq_i, ai_sq);
        float term2 = __fmul_rn(sum_sq_j, aj_sq);
        float term3 = __fmul_rn(2.0f, __fmul_rn(aij, C));
        
        float S         = __fadd_rn(term1, __fsub_rn(term2, term3));
        float a_sq      = __fmul_rn(a, a);
        float denom     = __fmul_rn(a_sq, static_cast<float>(n - 1));
        float lrv_value = __fdiv_rn(S, denom);

        int pair_index = (i * (i - 1)) / 2 + j;
        d_variances[pair_index] = lrv_value;
    }
}



__global__
void
propr::detail::cuda::lrv_alpha_weighted(float* __restrict__ d_Y,
                             float* __restrict__ d_Yfull,
                             float* __restrict__ d_W,
                             float* __restrict__ d_Wfull,
                             float a,
                             float* __restrict__ d_variances,
                             int nb_samples, int nb_genes) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i < nb_genes && j < i) {
        float sum_w    = 0.0f;
        float sum_w_sq = 0.0f;

        float sum_wX_i = 0.0f;
        float sum_wX_j = 0.0f;
        
        float sum_wX_i_sq = 0.0f;
        float sum_wX_j_sq = 0.0f;

        float sum_wX_iX_j = 0.0f;

        float sum_w_full  = 0.0f;

        float sum_w_full_X_full_i = 0.0f; 
        float sum_w_full_X_full_j = 0.0f;
        int k = 0;

        #pragma unroll
        for (; k <= nb_samples - 4; k += 4) {
            float4 y_i     = *reinterpret_cast<float4*>(&d_Y[k + i * nb_samples]);
            float4 y_j     = *reinterpret_cast<float4*>(&d_Y[k + j * nb_samples]);
            float4 yfull_i = *reinterpret_cast<float4*>(&d_Yfull[k + i * nb_samples]);
            float4 yfull_j = *reinterpret_cast<float4*>(&d_Yfull[k + j * nb_samples]);
            float4 w_i     = *reinterpret_cast<float4*>(&d_W[k + i * nb_samples]);
            float4 w_j     = *reinterpret_cast<float4*>(&d_W[k + j * nb_samples]);
            float4 wfull_i = *reinterpret_cast<float4*>(&d_Wfull[k + i * nb_samples]);
            float4 wfull_j = *reinterpret_cast<float4*>(&d_Wfull[k + j * nb_samples]);

            #pragma unroll
            for (int m = 0; m < 4; ++m) {
                float X_i      = __powf(reinterpret_cast<float*>(&y_i)[m]    , a);
                float X_j      = __powf(reinterpret_cast<float*>(&y_j)[m]    , a);
                float X_full_i = __powf(reinterpret_cast<float*>(&yfull_i)[m], a);
                float X_full_j = __powf(reinterpret_cast<float*>(&yfull_j)[m], a);
                
                float w      = __fmul_rn(reinterpret_cast<float*>(&w_i)[m]    , reinterpret_cast<float*>(&w_j)[m]);
                float w_full = __fmul_rn(reinterpret_cast<float*>(&wfull_i)[m], reinterpret_cast<float*>(&wfull_j)[m]);
                
                sum_w    = __fadd_rn(sum_w, w);
                sum_w_sq = __fmaf_rn(w, w, sum_w_sq);
                
                sum_wX_i = __fmaf_rn(w, X_i, sum_wX_i);
                sum_wX_j = __fmaf_rn(w, X_j, sum_wX_j);
                
                float X_i_sq = __fmul_rn(X_i, X_i);
                float X_j_sq = __fmul_rn(X_j, X_j);
                sum_wX_i_sq  = __fmaf_rn(w, X_i_sq, sum_wX_i_sq);
                sum_wX_j_sq  = __fmaf_rn(w, X_j_sq, sum_wX_j_sq);
                
                sum_wX_iX_j = __fmaf_rn(w, __fmul_rn(X_i, X_j), sum_wX_iX_j);
                
                sum_w_full          = __fadd_rn(sum_w_full, w_full);
                sum_w_full_X_full_i = __fmaf_rn(w_full, X_full_i, sum_w_full_X_full_i);
                sum_w_full_X_full_j = __fmaf_rn(w_full, X_full_j, sum_w_full_X_full_j);
            }
        }

        for (; k < nb_samples; ++k) {
            float X_i      = __powf(d_Y[k + i * nb_samples], a);
            float X_j      = __powf(d_Y[k + j * nb_samples], a);
            float X_full_i = d_Yfull[k + i * nb_samples];
            float X_full_j = d_Yfull[k + j * nb_samples];
            
            float w = __fmul_rn(d_W[k + i * nb_samples], d_W[k + j * nb_samples]);
            float w_full = __fmul_rn(d_Wfull[k + i * nb_samples], d_Wfull[k + j * nb_samples]);
            
            sum_w    = __fadd_rn(sum_w, w);
            sum_w_sq = __fmaf_rn(w, w, sum_w_sq);
            
            sum_wX_i = __fmaf_rn(w, X_i, sum_wX_i);
            sum_wX_j = __fmaf_rn(w, X_j, sum_wX_j);
            
            float X_i_sq = __fmul_rn(X_i, X_i);
            float X_j_sq = __fmul_rn(X_j, X_j);
            sum_wX_i_sq  = __fmaf_rn(w, X_i_sq, sum_wX_i_sq);
            sum_wX_j_sq  = __fmaf_rn(w, X_j_sq, sum_wX_j_sq);
            
            sum_wX_iX_j = __fmaf_rn(w, __fmul_rn(X_i, X_j), sum_wX_iX_j);
            
            sum_w_full = __fadd_rn(sum_w_full, w_full);
            sum_w_full_X_full_i = __fmaf_rn(w_full, X_full_i, sum_w_full_X_full_i);
            sum_w_full_X_full_j = __fmaf_rn(w_full, X_full_j, sum_w_full_X_full_j);
        }

        float inv_sum_w = __frcp_rn(sum_w);
        float inv_sum_w_full = __frcp_rn(sum_w_full);
        
        float w_valid = __fdiv_rn(fminf(fabsf(sum_w), FLT_MAX), __fadd_rn(fabsf(sum_w), FLT_MIN));
        
        // float mu_i = __fmul_rn(sum_wX_i, __fmul_rn(inv_sum_w, w_valid));
        // float mu_j = __fmul_rn(sum_wX_j, __fmul_rn(inv_sum_w, w_valid));
        
        float sum_wX_i_ratio = __fmul_rn(sum_wX_i, sum_wX_i);
        float sum_wX_j_ratio = __fmul_rn(sum_wX_j, sum_wX_j);
        
        float sum_sq_i = __fsub_rn(sum_wX_i_sq, __fmul_rn(sum_wX_i_ratio, inv_sum_w));
        float sum_sq_j = __fsub_rn(sum_wX_j_sq, __fmul_rn(sum_wX_j_ratio, inv_sum_w));
        float C = __fsub_rn(sum_wX_iX_j, __fmul_rn(__fmul_rn(sum_wX_i, sum_wX_j), inv_sum_w));
        
        float wfull_valid = __fdiv_rn(fminf(fabsf(sum_w_full), FLT_MAX), __fadd_rn(fabsf(sum_w_full), FLT_MIN));
        
        float m_i = __fmul_rn(sum_w_full_X_full_i, __fmul_rn(inv_sum_w_full, wfull_valid));
        float m_j = __fmul_rn(sum_w_full_X_full_j, __fmul_rn(inv_sum_w_full, wfull_valid));
        
        float m_i_valid = __fdiv_rn(fminf(fabsf(m_i), FLT_MAX), __fadd_rn(fabsf(m_i), FLT_MIN));
        float m_j_valid = __fdiv_rn(fminf(fabsf(m_j), FLT_MAX), __fadd_rn(fabsf(m_j), FLT_MIN));
        
        float a_i = __fmul_rn(__frcp_rn(m_i), m_i_valid);
        float a_j = __fmul_rn(__frcp_rn(m_j), m_j_valid);
        
        float a_i_sq = __fmul_rn(a_i, a_i);
        float a_j_sq = __fmul_rn(a_j, a_j);
        
        float term1 = __fmul_rn(sum_sq_i, a_i_sq);
        float term2 = __fmul_rn(sum_sq_j, a_j_sq);
        float term3 = __fmul_rn(2.0f, __fmul_rn(__fmul_rn(a_i, a_j), C));
        
        float numerator = __fadd_rn(term1, __fsub_rn(term2, term3));
        float denominator_term = __fsub_rn(sum_w, __fmul_rn(sum_w_sq, inv_sum_w));
        
        float denom_valid = __fdiv_rn(fminf(denominator_term, FLT_MAX),  __fadd_rn(denominator_term, FLT_MIN));
        float denominator = __fmul_rn(__fmul_rn(a, a), denominator_term);
        float lrv_value   = __fmul_rn(__fdiv_rn(numerator, denominator), denom_valid);
        
        int pair_index = (i * (i - 1)) / 2 + j;
        d_variances[pair_index] = lrv_value;
    }
}
