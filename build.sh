#!/bin/bash
# Script to build all components from scratch, using the maximum available CPU power
#
# Given parameters are passed over to CMake.
# Examples:
#    * ./build_all.sh -DCMAKE_BUILD_TYPE=Debug
#    * ./build_all.sh VERBOSE=1
#
# Written by Tiffany Huang, 12/14/2016
#

# Go into the directory where this bash script is contained.
cd `dirname $0`

# Compile code.
mkdir -p build
cd build
cmake -G"Eclipse CDT4 - Unix Makefiles" -D CMAKE_BUILD_TYPE=Debug ../src/
make -j `nproc` $*
