cmake_minimum_required(VERSION 3.14)
project(obliczenia_naukowe)
set(CMAKE_CXX_STANDARD 14)
find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)
include_directories(${BLAS_INCLUDE_DIR})

add_executable(obliczenia_naukowe main.cpp nellspsv.cpp)

set(LAPACK_LIBRARIES /usr/lib/x86_64-linux-gnu/liblapacke.so ${LAPACK_LIBRARIES})

message(${LAPACK_LIBRARIES})

target_link_libraries(obliczenia_naukowe ${LAPACK_LIBRARIES})

