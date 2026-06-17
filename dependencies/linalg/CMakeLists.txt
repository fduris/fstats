# Get the macros and functions we'll need
include(FetchContent)

# Fetch the proper content
FetchContent_Declare(
    linalg
    GIT_REPOSITORY "https://github.com/jchristopherson/linalg"
)

FetchContent_MakeAvailable(linalg)

# Keep the default unversioned SONAME (liblinalg.so) so consumers record a
# relative DT_NEEDED, and resolve siblings (libferror.so) from the library's own
# directory via an $ORIGIN runpath. Both let the shipped .so load from the repo
# root on any machine without a post-build patchelf step.
set_target_properties(linalg PROPERTIES
    OUTPUT_NAME "linalg"
    INSTALL_RPATH "$ORIGIN"
    BUILD_WITH_INSTALL_RPATH TRUE
    LIBRARY_OUTPUT_DIRECTORY ${MT_DIR}
    RUNTIME_OUTPUT_DIRECTORY ${MT_DIR}
)

set(linalg_INCLUDE_DIR ${linalg_BINARY_DIR}/include)
set(linalg_INCLUDE_DIR ${linalg_INCLUDE_DIR} PARENT_SCOPE)
configure_file(
    "${linalg_SOURCE_DIR}/include/linalg.h"
    "${linalg_INCLUDE_DIR}/linalg.h"
    COPYONLY
)

# Make a parent-scope variable for the library
set(linalg_LIBRARY linalg)
set(linalg_LIBRARY ${linalg_LIBRARY} PARENT_SCOPE)
