#pragma once

namespace propr {
    namespace cuda {
        namespace traits {
           struct ctzRcpp_config {
            static constexpr int PHASE_ONE_BLK_X = 256;
            
            static constexpr int PHASE_TWO_BLK_X = 16;
            static constexpr int PHASE_TWO_BLK_Y = 16;
           };
        }
    }
}