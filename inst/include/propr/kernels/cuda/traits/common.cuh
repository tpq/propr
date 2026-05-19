#pragma once

namespace propr {
    namespace cuda {
        namespace traits {
            
            template<int threads_x = 256>
            struct thread_layout_1d {
                static constexpr int BLK_X = threads_x;
            };

            template<int threads_x = 16, int threads_y = 16>
            struct thread_layout_2d {
                static constexpr int BLK_X = threads_x;
                static constexpr int BLK_Y = threads_y;
            };
        }
    }
}
            