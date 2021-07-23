# CMakeList.txt : CMake project for DeepRT_M, include source and define
# project specific logic here.

cmake_minimum_required (VERSION 3.8)

project(RaytracingEngine)

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

enable_cxx_compiler_flag_if_supported("-Wall")
enable_cxx_compiler_flag_if_supported("-Wextra")
enable_cxx_compiler_flag_if_supported("-pedantic")
enable_cxx_compiler_flag_if_supported("-std=c++17")
enable_cxx_compiler_flag_if_supported("-O2")
enable_cxx_compiler_flag_if_supported("-fopenmp")

set(CMAKE_CXX_STANDARD_LIBRARIES "-static-libgcc -static-libstdc++ ${CMAKE_CXX_STANDARD_LIBRARIES}")

set(SRC_DIR RaytracingEngine)
set(LIB_DIR lib)

file(GLOB_RECURSE SOURCES ${SRC_DIR}/*.h ${SRC_DIR}/*.cpp)
set(LIB_DIRS ${LIB_DIR}/eigen)
set(LIB_NAMES eigen)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")

add_executable(${PROJECT_NAME} ${SOURCES})

add_subdirectory(${LIB_DIRS})

target_link_libraries(${PROJECT_NAME} PUBLIC ${LIB_NAMES})

target_include_directories(${PROJECT_NAME} PUBLIC
                           "${PROJECT_BINARY_DIR}"
						   "${PROJECT_SOURCE_DIR}/DeepRT"
                           )

# TODO: Add tests and install targets if needed.
