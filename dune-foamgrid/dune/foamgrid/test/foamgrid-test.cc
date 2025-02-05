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

#include "make2din3dgrid.hh"

template<class G>
void traversal (G& grid)
{
  // first we extract the dimensions of the grid
  const int dimgrid = G::dimension;

  // Leaf Traversal
  std::cout << "*** Traverse codim 0 leaves" << std::endl;

  // get the instance of the LeafGridView
  auto leafGridView = grid.leafGridView();

  // iterate through all entities of codim 0 at the leaves
  int count = 0;
  for (const auto& element : elements(leafGridView))
  {
    Dune::GeometryType gt = element.type();
    std::cout << "visiting leaf " << gt
              << " with first vertex at " << element.geometry().corner(0)
              << std::endl;
    count++;
  }

  std::cout << "there are/is " << count << " leaf element(s)" << std::endl;

  // Leafwise traversal of codim dim
  std::cout << std::endl;
  std::cout << "*** Traverse codim " << dimgrid << " leaves" << std::endl;

  // iterate through all entities of codim 0 on the given level
  count = 0;
  for (const auto& vertex : vertices(leafGridView))
  {
    Dune::GeometryType gt = vertex.type();
    std::cout << "visiting " << gt
              << " at " << vertex.geometry().corner(0)
              << std::endl;
    count++;
  }
  std::cout << "there are/is " << count << " leaf vertices(s)"
            << std::endl;

  // Levelwise traversal of codim 0
  std::cout << std::endl;
  std::cout << "*** Traverse codim 0 level-wise" << std::endl;

  // iterate through all entities of codim 0 on the given level
  for (int level=0; level<=grid.maxLevel(); level++)
  {
    // get the instance of the LeafGridView
    auto levelGridView = grid.levelGridView(level);

    count = 0;
    for (const auto& element : elements(levelGridView))
    {
      Dune::GeometryType gt = element.type();
      std::cout << "visiting " << gt
                << " with first vertex at " << element.geometry().corner(0)
                << std::endl;
      count++;
    }
    std::cout << "there are/is " << count << " element(s) on level "
              << level << std::endl;
    std::cout << std::endl;
  }
  // Iterate over all intersections
  std::cout << std::endl;
  std::cout << "*** Traverse intersections with level iterator" << std::endl;
  auto levelGridView = grid.levelGridView(0);

  for(const auto& element : elements(levelGridView))
  {
    Dune::GeometryType gt = element.type();
    std::cout << "visiting leaf " << gt
              << " with first vertex at " << element.geometry().corner(0)
              << " and second vertex at " << element.geometry().corner(1);
    if(dimgrid==2)
        std::cout << " and third vertex at " << element.geometry().corner(2);
    std::cout << std::endl;


    count = 0;
    for (const auto& intersection : intersections(levelGridView, element))
    {
        if(intersection.neighbor()) {
            std::cout << "found neighbor with first vertex at: "
                      << intersection.outside().geometry().corner(0) << " and second vertex at: "
                      << intersection.outside().geometry().corner(1);
            if(dimgrid==2)
                std::cout << " and third vertex at " << intersection.outside().geometry().corner(2);
            std::cout << std::endl;
            ++count;
        } else if(intersection.boundary()) {
            std::cout << "    this is a boundary intersection." << std::endl;
        }
    }
    std::cout << "This element knows about " << count << " neighbors." << std::endl << std::endl;
  }

  // Iterate over all intersections
  std::cout << std::endl;
  std::cout << "*** Traverse intersections with leaf iterator" << std::endl;

  for(const auto& element : elements(leafGridView))
  {
    Dune::GeometryType gt = element.type();
    std::cout << "visiting leaf " << gt
              << " with first vertex at " << element.geometry().corner(0)
              << " and second vertex at " << element.geometry().corner(1);
    if(dimgrid==2)
        std::cout << " and third vertex at " << element.geometry().corner(2);
    std::cout << std::endl;

    count = 0;
    for (const auto& intersection : intersections(leafGridView, element))
    {
        if(intersection.neighbor()){
            std::cout << "    found neighbor with first vertex at: "
                      << intersection.outside().geometry().corner(0)
                      << " and second vertex at: "
                      << intersection.outside().geometry().corner(1);
            if(dimgrid==2)
                std::cout << " and third vertex at " << intersection.outside().geometry().corner(2);
            std::cout << std::endl;
            ++count;
        } else if(intersection.boundary()) {
            std::cout << "    this is a boundary intersection." << std::endl;
        }
    }
    std::cout << "This element knows about " << count << " neighbors." << std::endl << std::endl;
  }
}

int main (int argc, char *argv[]) try
{
    Dune::MPIHelper::instance(argc, argv);

    // paths to gmsh test files
    const std::string dune_grid_path = std::string(DUNE_GRID_EXAMPLE_GRIDS_PATH) + "gmsh/";
    const std::string dune_foamgrid_path = std::string(DUNE_FOAMGRID_EXAMPLE_GRIDS_PATH) + "gmsh/";

    {
        std::cout << "\n################################################\n";
        std::cout << "Checking default-constructed (empty) FoamGrid<1, 3>\n";
        std::cout << "################################################\n\n";

        std::cout << "  Creating grid" << std::endl;
        FoamGrid<1, 3> emptyGrid;

        std::cout << "  Calling gridcheck" << std::endl;
        gridcheck(emptyGrid);

        std::cout << "  Calling checkIntersectionIterator" << std::endl;
        checkIntersectionIterator(emptyGrid);

        std::cout << "  Check if has multiple neighbor functionality" << std::endl;
        traversal(emptyGrid);
    }
    {
        std::cout << "\n################################################\n";
        std::cout << "Checking FoamGrid<1, 3> created from empty grid factory\n";
        std::cout << "################################################\n\n";

        std::cout << "  Creating grid" << std::endl;
        GridFactory<FoamGrid<1, 3>> factory;
        std::shared_ptr<FoamGrid<1, 3> > emptyGrid(factory.createGrid());

        std::cout << "  Calling gridcheck" << std::endl;
        gridcheck(*emptyGrid);

        std::cout << "  Calling checkIntersectionIterator" << std::endl;
        checkIntersectionIterator(*emptyGrid);

        std::cout << "  Check if has multiple neighbor functionality" << std::endl;
        traversal(*emptyGrid);
    }
    {
        std::cout << "\n################################################\n";
        std::cout << "Checking FoamGrid<2, 2> (2d in 2d grid)\n";
        std::cout << "################################################\n\n";

        std::cout << "  Creating grid" << std::endl;
        std::shared_ptr<FoamGrid<2, 2> > grid2d( GmshReader<FoamGrid<2, 2> >::read(dune_grid_path + "curved2d.msh", /*verbose*/ true, false ) );

        std::cout << "  Calling gridcheck" << std::endl;
        gridcheck(*grid2d);

        std::cout << "  Calling checkIntersectionIterator" << std::endl;
        checkIntersectionIterator(*grid2d);
    }
    {
        std::cout << "\n################################################\n";
        std::cout << "Checking FoamGrid<2, 3> (2d in 3d grid)\n";
        std::cout << "################################################\n\n";

        std::cout << "  Creating grid" << std::endl;
        std::shared_ptr<FoamGrid<2, 3>> grid3d( make2Din3DHybridTestGrid<FoamGrid<2, 3> >() );

        std::cout << "  Calling gridcheck" << std::endl;
        gridcheck(*grid3d);

        std::cout << "  Calling checkIntersectionIterator" << std::endl;
        checkIntersectionIterator(*grid3d);
    }
    {
        std::cout << "\n################################################\n";
        std::cout << "Checking FoamGrid<1, 1> (1d in 1d grid)\n";
        std::cout << "################################################\n\n";

        std::cout << "  Creating grid" << std::endl;
        std::shared_ptr<FoamGrid<1, 1> > grid11( GmshReader<FoamGrid<1, 1> >::read(dune_foamgrid_path + "line1d2d.msh", /*verbose*/ true, false ) );

        std::cout << "  Calling gridcheck" << std::endl;
        gridcheck(*grid11);

        std::cout << "  Calling checkIntersectionIterator" << std::endl;
        checkIntersectionIterator(*grid11);
    }
    {
        std::cout << "\n################################################\n";
        std::cout << "Checking FoamGrid<1, 2> (1d in 2d grid)\n";
        std::cout << "################################################\n\n";

        std::cout << "  Creating grid" << std::endl;
        std::shared_ptr<FoamGrid<1, 2> > grid12( GmshReader<FoamGrid<1, 2> >::read(dune_foamgrid_path + "line1d2d.msh", /*verbose*/ true, false ) );

        std::cout << "  Calling gridcheck" << std::endl;
        gridcheck(*grid12);

        std::cout << "  Calling checkIntersectionIterator" << std::endl;
        checkIntersectionIterator(*grid12);
    }
    {
        std::cout << "\n################################################\n";
        std::cout << "Checking FoamGrid<1, 3> (1d in 3d grid)\n";
        std::cout << "################################################\n\n";

        std::cout << "  Creating grid" << std::endl;
        std::shared_ptr<FoamGrid<1, 3> > grid13( GmshReader<FoamGrid<1, 3> >::read(dune_foamgrid_path + "line1d3d.msh", /*verbose*/ true, false ) );

        std::cout << "  Calling gridcheck" << std::endl;
        gridcheck(*grid13);

        std::cout << "  Calling checkIntersectionIterator" << std::endl;
        checkIntersectionIterator(*grid13);
    }
    {
        std::cout << "\n################################################\n";
        std::cout << "Checking FoamGrid<2, 3> (2d in 3d grid)\n";
        std::cout << "################################################\n\n";

        // dimworld == 3,  and a grid containing a T-Junction
        std::cout << "  Creating grid" << std::endl;
        std::shared_ptr<FoamGrid<2, 3> > gridTJunction( GmshReader<FoamGrid<2, 3> >::read(dune_foamgrid_path + "tjunction2d.msh", /*verbose*/ true, false ) );

        std::cout << "  Calling gridcheck" << std::endl;
        gridcheck(*gridTJunction);

        std::cout << "  Calling checkIntersectionIterator" << std::endl;
        checkIntersectionIterator(*gridTJunction);
    }
    {
        std::cout << "\n################################################\n";
        std::cout << "Checking FoamGrid<1, 3> (1d in 3d grid)\n";
        std::cout << "################################################\n\n";

        std::cout << "  Creating grid" << std::endl;
        std::shared_ptr<FoamGrid<1, 3> > gridStar( GmshReader<FoamGrid<1, 3> >::read(dune_foamgrid_path + "bifurcation1d3d.msh", /*verbose*/ true, false ) );

        std::cout << "  Calling gridcheck" << std::endl;
        gridcheck(*gridStar);

        std::cout << "  Calling checkIntersectionIterator" << std::endl;
        checkIntersectionIterator(*gridStar);

        std::cout << "  Check if has multiple neighbor functionality" << std::endl;
        traversal(*gridStar);
    }
    {
        std::cout << "\n################################################\n";
        std::cout << "Checking FoamGrid<2, 3, float> (2d in 3d grid with float coordinates)\n";
        std::cout << "################################################\n\n";

        std::cout << "  Creating grid" << std::endl;
        std::shared_ptr<FoamGrid<2, 3, float>> grid3d( make2Din3DHybridTestGrid<FoamGrid<2, 3, float> >() );

        std::cout << "  Calling gridcheck" << std::endl;
        gridcheck(*grid3d);

        std::cout << "  Calling checkIntersectionIterator" << std::endl;
        checkIntersectionIterator(*grid3d);
    }
}
// //////////////////////////////////
//   Error handler
// /////////////////////////////////
catch (const Exception& e) {
    std::cout << e << std::endl;
    return 1;
}
