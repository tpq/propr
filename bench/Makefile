# Benchmark Makefile
BUILD_DIR ?= build
NVBENCH_DIR = ../prereq/nvbench

# Compiler settings
CXX = g++
NVCC = nvcc

# Flags
CXXFLAGS = -std=c++17 -O3 -Wall -Wextra
NVCCFLAGS = -std=c++17 -O3 -arch=sm_70,sm_75,sm_80,sm_86

# Include paths
INCLUDES = -Iprereq/nvbench -I../src -Iprereq/nvbench/nvbench

# Library paths
LDFLAGS = -Lprereq/nvbench/build -lcuda -lcudart

# Source files
CU_SOURCES = src/getFdrRandcounts_Propd.cu src/lrm.cu src/lrv.cu src/omega.cu
CPP_SOURCES = src/backend.cpp src/comparison.cpp src/ctzRcpp.cpp src/graflex.cpp src/lr2propr.cpp src/lrm.cpp src/lrv.cpp src/omega.cpp src/RcppExports.cpp
CU_OBJECTS = 
CPP_OBJECTS = 
OBJECTS =  

# Target
TARGET = build/benchmark_executable

.PHONY: all clean

all: 

:  | build
	nvcc   -o bench-makefile 

build/%.o: src/%.cu | build
	nvcc   -c  -o bench-makefile

build/%.o: src/%.cpp | build
	g++   -c  -o bench-makefile

build:
	mkdir -p build

clean:
	rm -rf build
