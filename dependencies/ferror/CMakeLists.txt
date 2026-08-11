# Get the macros and functions we'll need
include(FetchContent)

# Fetch the proper content
# Pinned, not floating: MT commits the built artifacts to its repository root,
# so an unpinned fetch would silently rebase the shipped binaries onto whatever
# upstream had pushed that morning and make the build unreproducible. This is
# the revision the shipped libferror.so/.dll were built from.
FetchContent_Declare(
    ferror
    GIT_REPOSITORY "https://github.com/jchristopherson/ferror"
    GIT_TAG 9d6c9750f5c14ef2c806104ea8dfdfc9ce91efed
)

FetchContent_MakeAvailable(ferror)

# Keep the default unversioned SONAME (libferror.so) so consumers record a
# relative DT_NEEDED, and resolve siblings from the library's own directory via
# an $ORIGIN runpath. Both let the shipped .so load from the repo root on any
# machine without a post-build patchelf step.
set_target_properties(ferror PROPERTIES
    OUTPUT_NAME "ferror"
    INSTALL_RPATH "$ORIGIN"
    BUILD_WITH_INSTALL_RPATH TRUE
    LIBRARY_OUTPUT_DIRECTORY ${MT_DIR}
    RUNTIME_OUTPUT_DIRECTORY ${MT_DIR}
)

set(ferror_INCLUDE_DIR ${ferror_BINARY_DIR}/include)
set(ferror_INCLUDE_DIR ${ferror_INCLUDE_DIR} PARENT_SCOPE)
configure_file(
    "${ferror_SOURCE_DIR}/include/ferror.h"
    "${ferror_INCLUDE_DIR}/ferror.h"
    COPYONLY
)

# Make a parent-scope variable for the library
set(ferror_LIBRARY ferror)
set(ferror_LIBRARY ${ferror_LIBRARY} PARENT_SCOPE)
