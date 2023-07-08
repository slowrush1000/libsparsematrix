#!/bin/tcsh -f

set CMAKE_C_COMPILER    = "/usr/bin/gcc"
set CMAKE_CXX_COMPILER  = "/usr/bin/g++"
#set CMAKE_C_COMPILER    = "/usr/bin/clang"
#set CMAKE_CXX_COMPILER  = "/usr/bin/clang++"
#set CMAKE_C_COMPILER    = "/opt/intel/oneapi/compiler/2022.1.0/linux/bin/intel64/icc"
#set CMAKE_CXX_COMPILER  = "/opt/intel/oneapi/compiler/2022.1.0/linux/bin/intel64/icpc"

\rm -rf debug release
mkdir debug release

cd debug
cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} ../..
#cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} ../..
cd ..

cd release
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} ../..
#cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} ../..
cd .. 
