#ifndef DUNE_FOAMGRID_GEOMETRY_HH
#define DUNE_FOAMGRID_GEOMETRY_HH

/** \file
* \brief The FoamGridGeometry class
*/

#include <vector>
#include <array>

#include <dune/common/version.hh>
#include <dune/geometry/affinegeometry.hh>


namespace Dune {


template<int mydim, int coorddim, class GridImp>
class FoamGridGeometry :
        public AffineGeometry<typename GridImp::ctype, mydim, coorddim>
{
    public:

    /**
     * \brief This is DefaultConstructor
     */
    FoamGridGeometry() {}

    /**
     * \brief Construct geometry from coordinate vector
     */
    FoamGridGeometry(const GeometryType& type, const std::vector<FieldVector<typename GridImp::ctype,coorddim> >& coordinates) :
        AffineGeometry<typename GridImp::ctype, mydim, coorddim>(type, coordinates)
    {}

    /**
     * \brief Construct geometry from coordinate array
     * \note more efficient if corner size is known at compile time, e.g. for simplices
     */
    template<std::size_t size>
    FoamGridGeometry(const GeometryType& type, const std::array<FieldVector<typename GridImp::ctype,coorddim>, size>& coordinates) :
        AffineGeometry<typename GridImp::ctype, mydim, coorddim>(type, coordinates)
    {}
};


}  // namespace Dune

#endif
