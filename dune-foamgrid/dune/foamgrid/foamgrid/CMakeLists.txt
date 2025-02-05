# install headers
file(GLOB DUNE_FOAMGRID_HEADERS *.hh)
install(FILES foamgrid.cc ${DUNE_FOAMGRID_HEADERS}
        DESTINATION include/dune/foamgrid/foamgrid)
