from __future__ import absolute_import, unicode_literals

def foamGrid(constructor, dimgrid, dimworld):
    from dune.grid.grid_generator import module

    typeName = "Dune::FoamGrid<{}, {}>".format(dimgrid, dimworld)
    includes = ["dune/foamgrid/foamgrid.hh", "dune/foamgrid/dgffoam.hh"]
    gridModule = module(includes, typeName)

    return gridModule.reader(constructor).leafView
