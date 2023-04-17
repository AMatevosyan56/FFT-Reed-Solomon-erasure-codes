cmake_minimum_required(VERSION 3.22)
project(fft)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20")

set(SOURCE_FILES RSErasureCode.cpp)
add_executable(fft ${SOURCE_FILES})

install(TARGETS fft DESTINATION bin)
