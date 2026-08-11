cmake_minimum_required(VERSION 3.26)

project(fstats Fortran)
enable_language(Fortran)

set(CMAKE_Fortran_STANDARD 2018)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)

set(MT_DIR "${CMAKE_BINARY_DIR}/../../..")

# Release artifacts ship stripped. MT reaches every library here across an FFI
# boundary that resolves through the export table alone, so the symbol table is
# dead weight in binaries that are committed to the repository.
#
# `-s` is a link-time flag: under target_compile_options it is handed to the
# compiler only and silently does nothing, which is why the shipped libraries
# were unstripped despite the config asking for it. It has to travel in the
# linker flags instead. It is set here, ahead of the dependency subdirectories,
# so that it also reaches ferror/linalg/collections and the BLAS/LAPACK targets
# LINALG fetches transitively -- those come from third-party CMakeLists this
# project does not own, so there is no target of ours to hang
# target_link_options on, and a per-target flag would strip fstats alone.
set(CMAKE_SHARED_LINKER_FLAGS_RELEASE "${CMAKE_SHARED_LINKER_FLAGS_RELEASE} -s")

# The win64 artifacts additionally ship with a zeroed build timestamp. GNU ld
# stamps the wall-clock build time into both the PE COFF header and the export
# directory, so two links of byte-identical objects still differ there -- and in
# the PE checksum that covers them -- which is enough to defeat an md5 check of
# a committed DLL against a rebuild. `--no-insert-timestamp` zeroes both copies.
# It is a PE-only option that the native ELF linker rejects outright, so it has
# to stay behind the MINGW guard rather than join the `-s` above; it sits at the
# same scope for the same reason, so that the dependency subdirectories inherit
# it too.
if (MINGW)
    set(CMAKE_SHARED_LINKER_FLAGS_RELEASE "${CMAKE_SHARED_LINKER_FLAGS_RELEASE} -Wl,--no-insert-timestamp")
endif()

add_subdirectory(configure)
add_subdirectory(dependencies)

find_package(OpenMP REQUIRED COMPONENTS Fortran)

set(FSTATS_SOURCES
    src/fstats.f90
    src/fstats_special_functions.f90
    src/fstats_descriptive_statistics.f90
    src/fstats_hypothesis.f90
    src/fstats_anova.f90
    src/fstats_distributions.f90
    src/fstats_helper_routines.f90
    src/fstats_regression.f90
    src/fstats_experimental_design.f90
    src/fstats_allan.f90
    src/fstats_errors.f90
    src/fstats_bootstrap.f90
    src/fstats_sampling.f90
    src/fstats_smoothing.f90
    src/fstats_types.f90
    src/fstats_mcmc.f90
    src/fstats_interp.f90
    src/mt_select.f90
    src/mt_expfit.f90
    src/mt_robust.f90
    src/mt_loess.f90
    src/mt_borovicka.f90
    src/vendored/lbfgsb.f
    src/vendored/linpack_subset.f
    src/vendored/timer_stub.f
)

set_source_files_properties(${FSTATS_SOURCES} PROPERTIES
    Fortran_PREPROCESS ON
)

# Fixed-form F77 cannot meet the -std=f2018 -Werror used for the rest of the
# library, and its comment syntax would choke the C preprocessor. Must follow
# the blanket Fortran_PREPROCESS ON above so the per-file OFF wins.
set_source_files_properties(
    src/vendored/lbfgsb.f src/vendored/linpack_subset.f src/vendored/timer_stub.f
    PROPERTIES COMPILE_OPTIONS "-std=legacy;-w" Fortran_PREPROCESS OFF
)

set(MT_SOURCES
    mt/mt_fstats.f90
)

add_library(fstats SHARED
    ${FSTATS_SOURCES}
    ${MT_SOURCES}
)

# Resolve the sibling dependency libraries (ferror/linalg/collections) from this
# library's own directory via an $ORIGIN runpath, so the shipped libfstats.so
# loads from the repo root on any machine without a post-build patchelf step.
set_target_properties(fstats PROPERTIES
    OUTPUT_NAME fstats
    POSITION_INDEPENDENT_CODE TRUE
    INSTALL_RPATH "$ORIGIN"
    BUILD_WITH_INSTALL_RPATH TRUE
    LIBRARY_OUTPUT_DIRECTORY ${MT_DIR}
    RUNTIME_OUTPUT_DIRECTORY ${MT_DIR}
)

target_include_directories(fstats PUBLIC
    $<BUILD_INTERFACE:${mod_dir}>
    $<BUILD_INTERFACE:${ferror_INCLUDE_DIR}>
    $<BUILD_INTERFACE:${linalg_INCLUDE_DIR}>
    $<BUILD_INTERFACE:${collections_INCLUDE_DIR}>
    $<INSTALL_INTERFACE:${include_install_dir}>
)

target_link_libraries(fstats PUBLIC
    ${ferror_LIBRARY} 
    ${linalg_LIBRARY} 
    ${collections_LIBRARY}
)

if (OpenMP_Fortran_FOUND)
    target_link_libraries(fstats PUBLIC OpenMP::OpenMP_Fortran)
    target_compile_definitions(fstats PRIVATE USEOPENMP=1)
endif()

target_compile_options(fstats PRIVATE
    $<$<CONFIG:DEBUG>:-O0 -g>
    $<$<CONFIG:RELEASE>:-O3>
    -std=f2018
    -Wall -Werror -pedantic -pedantic-errors
    -Wno-error=maybe-uninitialized
    -Wno-error=unused-variable
    -Wno-error=unused-dummy-argument
    -Wno-error=unused-label
    -Wno-error=surprising
    -Wno-error=uninitialized
    -fimplicit-none
)

if (WIN32)
    add_custom_command(TARGET fstats POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_BINARY_DIR}/_deps/lapack-build/bin/libblas.dll
        ${CMAKE_BINARY_DIR}/_deps/lapack-build/bin/liblapack.dll
        ${MT_DIR}
    )
endif()
