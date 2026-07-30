cmake_minimum_required(VERSION 3.20)
project(resonarsat VERSION 0.1.0 LANGUAGES C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

# ---------------------------------------------------------------------------
# Build types. ASAN is a first-class configuration, not an afterthought: the
# readers parse untrusted external binary formats and ctest runs under it.
# ---------------------------------------------------------------------------
set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG")
set(CMAKE_C_FLAGS_ASAN "-O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer")
set(CMAKE_EXE_LINKER_FLAGS_ASAN "-fsanitize=address,undefined")

# Deliberately NOT -ffast-math: it permits reassociation and flushes denormals,
# which perturbs sub-pixel correlation peaks and interferometric phase -- the
# exact quantities this project measures.
add_compile_options(-Wall -Wextra -Wshadow -Wconversion -Wno-sign-conversion)

option(RS_NATIVE_ARCH "Build with -march=native" OFF)
if(RS_NATIVE_ARCH)
    include(CheckCCompilerFlag)
    check_c_compiler_flag(-march=native RS_HAS_MARCH_NATIVE)
    if(RS_HAS_MARCH_NATIVE)
        add_compile_options(-march=native)
    endif()
endif()

# ---------------------------------------------------------------------------
# OpenMP is optional -- every parallel region is guarded by #ifdef _OPENMP and
# the code is correct single-threaded. But "optional" quietly became "absent"
# here: Apple clang ships without it, Homebrew installs libomp keg-only so it is
# not on any default search path, and find_package then reports failure while
# the build succeeds. The result is a serial binary that looks identical to a
# threaded one until you notice backprojection using a single core.
#
# So look where Homebrew actually puts it before giving up, and say loudly which
# outcome occurred.
# ---------------------------------------------------------------------------
if(APPLE AND NOT DEFINED OpenMP_ROOT)
    find_program(RS_BREW brew)
    if(RS_BREW)
        execute_process(COMMAND ${RS_BREW} --prefix libomp
                        OUTPUT_VARIABLE RS_LIBOMP_PREFIX
                        OUTPUT_STRIP_TRAILING_WHITESPACE
                        ERROR_QUIET RESULT_VARIABLE RS_BREW_RC)
        if(RS_BREW_RC EQUAL 0 AND EXISTS "${RS_LIBOMP_PREFIX}")
            set(OpenMP_ROOT "${RS_LIBOMP_PREFIX}")
            message(STATUS "Using Homebrew libomp at ${OpenMP_ROOT}")
        endif()
    endif()
endif()

find_package(OpenMP COMPONENTS C)

if(OpenMP_C_FOUND)
    message(STATUS "OpenMP enabled: backprojection and tracking will use all cores")
else()
    message(WARNING
        "OpenMP NOT found -- this build will be SINGLE-THREADED.\n"
        "Backprojection is the dominant cost and will run on one core.\n"
        "On macOS:  brew install libomp\n"
        "then reconfigure. The build is correct without it, only slow.")
endif()

# ---------------------------------------------------------------------------
# Vendored PocketFFT (BSD-3-Clause). Nothing outside src/core/fft.c may include
# or name it -- see include/resonarsat/fft.h.
# ---------------------------------------------------------------------------
add_library(pocketfft STATIC third_party/pocketfft/pocketfft.c)
target_include_directories(pocketfft PUBLIC third_party/pocketfft)
target_compile_options(pocketfft PRIVATE -w)  # third-party: not our warnings to fix

# ---------------------------------------------------------------------------
# Core library
# ---------------------------------------------------------------------------
add_library(resonarsat_core STATIC
    src/util/status.c
    src/util/geom.c
    src/util/geocode.c
    src/util/orbit.c
    src/core/fft.c
    src/core/focus.c
    src/core/subaperture.c
    src/core/coreg.c
    src/core/phaselink.c
    src/core/microm.c
    src/core/ccd.c
    src/core/spectrum.c
    src/core/tomo.c
    src/core/simulate.c
    src/readers/uavsar.c
    src/readers/cphd.c
    src/readers/sicd.c
    src/io/raster.c
    src/io/png.c
)
target_include_directories(resonarsat_core PUBLIC include)
target_link_libraries(resonarsat_core PRIVATE pocketfft)
target_link_libraries(resonarsat_core PUBLIC m)
if(OpenMP_C_FOUND)
    target_link_libraries(resonarsat_core PUBLIC OpenMP::OpenMP_C)
endif()

# ---------------------------------------------------------------------------
# CLI driver and tools
# ---------------------------------------------------------------------------
add_executable(resonarsat src/main.c)
target_link_libraries(resonarsat PRIVATE resonarsat_core)

add_executable(sim_cphd tools/sim_cphd.c)
target_link_libraries(sim_cphd PRIVATE resonarsat_core)

add_executable(sim_slc tools/sim_slc.c)
target_link_libraries(sim_slc PRIVATE resonarsat_core)

add_executable(crop_slc tools/crop_slc.c)
target_link_libraries(crop_slc PRIVATE resonarsat_core)

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------
enable_testing()

# CLI parsing must reject option tokens where a numeric value is required.
# Before strict parsing, atof("--clutter") silently produced zero and launched
# a long simulation with a different scene from the one requested.
add_test(NAME sim_cphd_rejects_missing_option_value
         COMMAND sim_cphd ignored.cphd --clutter)
set_tests_properties(sim_cphd_rejects_missing_option_value
                     PROPERTIES WILL_FAIL TRUE)
add_test(NAME sim_cphd_rejects_unknown_option
         COMMAND sim_cphd ignored.cphd --cluter 3000)
set_tests_properties(sim_cphd_rejects_unknown_option
                     PROPERTIES WILL_FAIL TRUE)

set(RS_TESTS
    test_constants
    test_fft
    test_focus
    test_subaperture
    test_coreg
    test_tomo
    test_endtoend
    test_validation
    test_depth
    test_tracking
    test_nullmotion
    test_ccd
    test_phaselink
    test_modelc
    test_cphd
    test_readers
    test_geocode
    test_raster
)
foreach(t IN LISTS RS_TESTS)
    add_executable(${t} tests/${t}.c)
    target_link_libraries(${t} PRIVATE resonarsat_core)
    target_include_directories(${t} PRIVATE tests)
    add_test(NAME ${t} COMMAND ${t})
endforeach()

# Documentation coverage: every function in src/ and tools/ must carry the
# block comment required by section 3a of the implementation plan.
find_package(Python3 COMPONENTS Interpreter)
if(Python3_Interpreter_FOUND)
    add_test(NAME docs_coverage
             COMMAND ${Python3_EXECUTABLE}
                     ${CMAKE_SOURCE_DIR}/scripts/check_docs.py
                     ${CMAKE_SOURCE_DIR}/src ${CMAKE_SOURCE_DIR}/tools)
endif()
