#pragma once

#define IS_POWER_OF_2(x) (((x) != 0) && (((x) & ((x) - 1)) == 0))
#define PROPR_WARP_SIZE 32U

#define HLF_EPSILON __float2half(4.887581E-04)
#define HLF_MIN     __float2half(6.103516E-05)
#define HLF_MAX     __float2half(6.550400E+04)
