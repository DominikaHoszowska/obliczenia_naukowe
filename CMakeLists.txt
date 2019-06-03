cmake_minimum_required(VERSION 3.14)
project(obliczenia_naukowe)
set(CMAKE_CXX_STANDARD 14)
find_package( BLAS REQUIRED )
find_package(LAPACK REQUIRED)
include_directories(${BLAS_INCLUDE_DIR})
include_directories(${LAPACK_INCLUDE_DIR})

add_executable(obliczenia_naukowe main.cpp nellspsv.cpp)

target_link_libraries(obliczenia_naukowe ${BLAS_LIBRARIES})

target_link_libraries(obliczenia_naukowe ${LAPACK_LIBRARIES})
