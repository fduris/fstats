include(FetchContent)

# Pin reference BLAS/LAPACK. We never declare these ourselves -- LINALG fetches
# them, unpinned, when it cannot find a system BLAS/LAPACK -- but FetchContent
# honours the first declaration of a given name, so declaring them here ahead of
# add_subdirectory(linalg) is what fixes their revision. Keep it above the
# add_subdirectory(linalg) line below or the pin goes silently dead: LINALG's own
# floating declaration would be the first one seen.
#
# This revision is contemporaneous with the ferror/linalg/collections pins, so the
# four together form one coherent snapshot. Unlike those three it does NOT
# reproduce the committed libblas.dll/liblapack.dll, which predate 2025-11-11 and
# so carry roughly nine months less upstream history -- 97 of the 101 BLAS/LAPACK
# routines LINALG imports have changed since, including real numerical edits to
# dtrsm/dtrsv. The committed pair is therefore deliberately NOT rebuilt from this
# pin; it is shipped as the CI-proven machine code with only its symbol table
# removed. Anyone regenerating those two artifacts from here is changing numerics
# and needs to re-validate, not just re-diff sizes.
FetchContent_Declare(
    lapack
    GIT_REPOSITORY "https://github.com/Reference-LAPACK/lapack"
    GIT_TAG 51b349470b3b26c948d068deb45c9b120a47ed32
)

# Get FERROR
add_subdirectory(ferror)
set(ferror_LIBRARY ${ferror_LIBRARY} PARENT_SCOPE)
set(ferror_INCLUDE_DIR ${ferror_INCLUDE_DIR} PARENT_SCOPE)

# Get LINALG
add_subdirectory(linalg)
set(linalg_LIBRARY ${linalg_LIBRARY} PARENT_SCOPE)
set(linalg_INCLUDE_DIR ${linalg_INCLUDE_DIR} PARENT_SCOPE)

# Get COLLECTIONS
add_subdirectory(collections)
set(collections_LIBRARY ${collections_LIBRARY} PARENT_SCOPE)
set(collections_INCLUDE_DIR ${collections_INCLUDE_DIR} PARENT_SCOPE)