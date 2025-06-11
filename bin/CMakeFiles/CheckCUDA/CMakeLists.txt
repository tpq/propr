cmake_minimum_required(VERSION 3.31.4)
set(CMAKE_MODULE_PATH "/gpfs/projects/bsc99/bsc099260/EuroHPC/PROPR/SRC/BASELINE/propr/cmake")
project(CheckCUDA CUDA)
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/result.cmake"
  "set(CMAKE_CUDA_COMPILER \"${CMAKE_CUDA_COMPILER}\")\n"
  "set(CMAKE_CUDA_HOST_COMPILER \"${CMAKE_CUDA_HOST_COMPILER}\")\n"
  )