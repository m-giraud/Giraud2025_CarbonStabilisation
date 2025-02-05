// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set ts=8 sw=4 et sts=4:
#ifndef DUNE_FOAMGRID_NULLITERATORFACTORY_HH
#define DUNE_FOAMGRID_NULLITERATORFACTORY_HH

/** \file
* \brief The null iterator factory for intersections
*/

#include <dune/foamgrid/foamgrid/foamgridvertex.hh>

namespace Dune {

template <int dimgrid, int dimworld, class ctype>
class FoamGridNullIteratorFactory
{
  public:
    static typename std::vector<const FoamGridEntityImp<dimgrid, dimgrid, dimworld, ctype>*>::const_iterator null()
    { return emptyVector_.end(); }
  private:
    static typename std::vector<const FoamGridEntityImp<dimgrid, dimgrid, dimworld, ctype>*> emptyVector_;
};

template <int dimgrid, int dimworld, class ctype>
typename std::vector<const FoamGridEntityImp<dimgrid, dimgrid, dimworld, ctype>*>
FoamGridNullIteratorFactory<dimgrid, dimworld, ctype>::emptyVector_;
}

#endif
