add_definitions(-DDUNE_FOAMGRID_EXAMPLE_GRIDS_PATH=\"${PROJECT_SOURCE_DIR}/doc/grids/\")
add_definitions(-DDUNE_GRID_EXAMPLE_GRIDS_PATH=\"${DUNE_GRID_EXAMPLE_GRIDS_PATH}\")

if(${DUNE_GRID_VERSION} VERSION_LESS 2.7)
  # before version 2.7 these tests fail because of an overly strict requirement
  # in the grid test about the sum over all outer normals of a grid view
  dune_add_test(SOURCES foamgrid-test.cc EXPECT_FAIL)
  dune_add_test(SOURCES local-refine-test.cc EXPECT_FAIL)
  # before version 2.7 this test fails to run because of a bug
  # in the GmshReader that failed to read boundary segments for 1d grids
  # (fixed in dune-grid commit eaadf1b1)
  dune_add_test(SOURCES boundary-segment-test.cc COMPILE_ONLY)
else()
  dune_add_test(SOURCES foamgrid-test.cc)
  dune_add_test(SOURCES local-refine-test.cc)
  dune_add_test(SOURCES boundary-segment-test.cc)
endif()

dune_add_test(SOURCES global-refine-test.cc)
dune_add_test(SOURCES growth-test-1d.cc)
dune_add_test(SOURCES growth-test-2d.cc)
dune_add_test(SOURCES setposition-test.cc)
