# CMakeList.txt : CMake project for RaytracingEngine, include source and define
# project specific logic here.

cmake_minimum_required (VERSION 3.0)

project(RaytracingEngine)

set(CMAKE_CXX_STANDARD_LIBRARIES "-static-libgcc -static-libstdc++ ${CMAKE_CXX_STANDARD_LIBRARIES}")

include(CheckCXXCompilerFlag)

function(enable_cxx_compiler_flag_if_supported flag)
    string(FIND "${CMAKE_CXX_FLAGS}" "${flag}" flag_already_set)
    if(flag_already_set EQUAL -1)
        check_cxx_compiler_flag("${flag}" flag_supported)
        if(flag_supported)
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${flag}" PARENT_SCOPE)
        endif()
        unset(flag_supported CACHE)
    endif()
endfunction()

if(MSVC)
    enable_cxx_compiler_flag_if_supported("/Yd")
    enable_cxx_compiler_flag_if_supported("/std:c++17")
    enable_cxx_compiler_flag_if_supported("/openmp")
elseif(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
    enable_cxx_compiler_flag_if_supported("-Wextra")
    enable_cxx_compiler_flag_if_supported("-pedantic")
    enable_cxx_compiler_flag_if_supported("-std=c++17")
    enable_cxx_compiler_flag_if_supported("-O1")
    enable_cxx_compiler_flag_if_supported("-O2")
    enable_cxx_compiler_flag_if_supported("-fopenmp")
endif()

set(SRC_DIR ${PROJECT_NAME})

file(GLOB_RECURSE SOURCES ${SRC_DIR}/*.h ${SRC_DIR}/*.cpp)
set(LIB_NAMES Eigen3::Eigen)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

add_executable(${PROJECT_NAME} ${SOURCES})

target_link_libraries(${PROJECT_NAME} ${LIB_NAMES})

target_include_directories(${PROJECT_NAME} PUBLIC
						   "${PROJECT_SOURCE_DIR}/${PROJECT_NAME}")