#include <config.h>

#include <iostream>
#include <memory>
#include <string>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/exceptions.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/test/gridcheck.hh>
#include <dune/grid/test/checkintersectionit.hh>
#include <dune/foamgrid/foamgrid.hh>

template<class Grid>
void checkBoundarySegments(const Grid& grid,
                           const Dune::GridFactory<Grid>& factory,
                           const std::vector<int>& boundaryMarkers,
                           int boundaryIntersectionsExpected,
                           int insertedBoundariesExpected,
                           int numIntersectionsExpected,
                           int boundaryIdExpected = 42)
{
  std::cout << "  Checking boundary segment indices" << std::endl;

  int numIntersections = 0;
  int boundaryIntersections = 0;
  int insertedBoundaries = 0;

  for ( const auto& element : elements(grid.leafGridView()) )
  {
    for ( const auto& intersection : intersections(grid.leafGridView(), element) )
    {
      numIntersections++;
      if ( intersection.boundary() )
      {
        boundaryIntersections++;
        if ( factory.wasInserted( intersection ) )
        {
          insertedBoundaries++;
          const int bndsegIdx = factory.insertionIndex(intersection);
          const int bndId = boundaryMarkers[bndsegIdx];
          if (bndId != boundaryIdExpected)
            DUNE_THROW(Dune::GridError, "Wrong boundary ID. Expected " << boundaryIdExpected << " got " << bndId);
        }
      }
    }
  }

  // sanity checks
  if (boundaryIntersections != boundaryIntersectionsExpected)
      DUNE_THROW(Dune::GridError, "Wrong number of boundary intersections. Expected " << boundaryIntersectionsExpected << " got " << boundaryIntersections);
  if (insertedBoundaries != insertedBoundariesExpected)
      DUNE_THROW(Dune::GridError, "Wrong number of inserted boundaries. Expected " << insertedBoundariesExpected << " got " << insertedBoundaries);
  if (numIntersections != numIntersectionsExpected)
      DUNE_THROW(Dune::GridError, "Wrong number of intersections. Expected " << numIntersectionsExpected << " got " << numIntersections);
}

int main (int argc, char *argv[]) try
{
    using namespace Dune;

    Dune::MPIHelper::instance(argc, argv);

    // paths to gmsh test files
    const std::string dune_foamgrid_path = std::string(DUNE_FOAMGRID_EXAMPLE_GRIDS_PATH) + "gmsh/";

    {
        std::cout << "\n################################################\n";
        std::cout << "Checking FoamGrid<1, 2> (1d in 2d grid)\n";
        std::cout << "################################################\n\n";

        std::cout << "  Creating grid" << std::endl;
        using Grid = FoamGrid<1, 2>;
        GridFactory<Grid> factory;

        std::vector<int> boundaryMarkers, elementMarkers;
        GmshReader<Grid>::read(factory, dune_foamgrid_path + "line1d2dbseg.msh", boundaryMarkers, elementMarkers, /*verbose*/ true, /*insertBoundarySegments*/ true);
        auto grid = std::shared_ptr<Grid>(factory.createGrid());

        // check the boundary segments
        checkBoundarySegments(*grid, factory, boundaryMarkers, 2, 1, 8);
    }

    {
        std::cout << "\n################################################\n";
        std::cout << "Checking FoamGrid<2, 3> (2d in 3d grid)\n";
        std::cout << "################################################\n\n";

        std::cout << "  Creating grid" << std::endl;
        using Grid = FoamGrid<2, 3>;
        GridFactory<Grid> factory;

        std::vector<int> boundaryMarkers, elementMarkers;
        GmshReader<Grid>::read(factory, dune_foamgrid_path + "junction2d3dbseg.msh", boundaryMarkers, elementMarkers, /*verbose*/ true, /*insertBoundarySegments*/ true);
        auto grid = std::shared_ptr<Grid>(factory.createGrid());

        // check the boundary segments
        checkBoundarySegments(*grid, factory, boundaryMarkers, 6, 3, 12);
    }
}
// //////////////////////////////////
//   Error handler
// /////////////////////////////////
catch (const Dune::Exception& e) {
    std::cout << e << std::endl;
    return 1;
}
