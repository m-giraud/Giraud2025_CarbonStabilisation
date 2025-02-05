add_subdirectory(dune)

configure_file(setup.py.in setup.py)

if (DUNE_COMMON_VERSION VERSION_GREATER_EQUAL 2.9)
dune_python_configure_bindings(
  PATH "."
  CMAKE_METADATA_FILE dune/foamgrid/metadata.cmake
  CMAKE_METADATA_FLAGS DUNE_OPTS_FILE
)
else()
dune_python_install_package(
  PATH "."
  CMAKE_METADATA_FILE dune/foamgrid/metadata.cmake
  CMAKE_METADATA_FLAGS DUNE_OPTS_FILE
)
endif()

if(POLICY CMP0087)
  cmake_policy(SET CMP0087 NEW)
endif()
