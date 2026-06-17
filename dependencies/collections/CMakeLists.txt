include(FetchContent)

# Fetch the proper content
FetchContent_Declare(
    collections
    GIT_REPOSITORY "https://github.com/jchristopherson/collections"
    GIT_TAG main
)

FetchContent_MakeAvailable(collections)

# Keep the default unversioned SONAME (libcollections.so) so consumers record a
# relative DT_NEEDED, and resolve siblings (libferror.so) from the library's own
# directory via an $ORIGIN runpath. Both let the shipped .so load from the repo
# root on any machine without a post-build patchelf step.
set_target_properties(collections PROPERTIES
    OUTPUT_NAME "collections"
    INSTALL_RPATH "$ORIGIN"
    BUILD_WITH_INSTALL_RPATH TRUE
    LIBRARY_OUTPUT_DIRECTORY ${MT_DIR}
    RUNTIME_OUTPUT_DIRECTORY ${MT_DIR}
)

set(collections_INCLUDE_DIR ${collections_BINARY_DIR}/src/mod_files)
set(collections_INCLUDE_DIR ${collections_INCLUDE_DIR} PARENT_SCOPE)

# Make a parent-scope variable for the library
set(collections_LIBRARY collections)
set(collections_LIBRARY ${collections_LIBRARY} PARENT_SCOPE)