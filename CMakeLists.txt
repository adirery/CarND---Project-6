project(Extended_and_non-stinky_KF)

cmake_minimum_required (VERSION 3.5)

set(sources
    src/FusionEKF.cpp
		src/ukf.cpp
    src/kalman_filter.cpp
    src/main.cpp
    src/tools.cpp)

add_executable(Extended_and_non-stinky_KF ${sources})
