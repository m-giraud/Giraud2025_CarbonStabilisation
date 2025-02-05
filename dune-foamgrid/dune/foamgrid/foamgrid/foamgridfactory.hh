// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set ts=8 sw=4 et sts=4:
#ifndef DUNE_FOAMGRID_FACTORY_HH
#define DUNE_FOAMGRID_FACTORY_HH

/** \file
    \brief The specialization of the generic GridFactory for FoamGrid
    \author Oliver Sander
 */

#include <vector>
#include <memory>
#include <unordered_map>
#include <functional>

#include <dune/common/deprecated.hh>
#include <dune/common/version.hh>
// TODO remove header and macro after release Dune 2.8
#define DUNE_FUNCTION_HH_SILENCE_DEPRECATION // silence deprecation warning from <dune/common/function.hh>
#include <dune/common/function.hh>
#include <dune/common/fvector.hh>
#include <dune/common/hash.hh>

#if DUNE_VERSION_EQUAL(DUNE_COMMON, 2, 7)
#include <dune/common/to_unique_ptr.hh>
#endif

#include <dune/grid/common/gridfactory.hh>
#include <dune/foamgrid/foamgrid.hh>

namespace Dune {

/** \brief Specialization of the generic GridFactory for FoamGrid<dimgrid, dimworld>
    */
template <int dimgrid, int dimworld, class ct>
    class GridFactoryBase
        : public GridFactoryInterface<FoamGrid<dimgrid, dimworld, ct> >
    {
    /** \brief Type used by the grid for coordinates */
    using ctype = ct;

    public:

#if DUNE_VERSION_LT(DUNE_COMMON, 2, 7)
    using GridPtrType = FoamGrid<dimgrid, dimworld, ctype>*;
#elif DUNE_VERSION_LT(DUNE_COMMON, 2, 8)
    using GridPtrType = ToUniquePtr<FoamGrid<dimgrid, dimworld, ctype>>;
#else
    using GridPtrType = std::unique_ptr<FoamGrid<dimgrid, dimworld,ctype>>;
#endif

    public:

        /** \brief Default constructor */
        GridFactoryBase()
        : factoryOwnsGrid_(true)
        {
            grid_ = new FoamGrid<dimgrid, dimworld, ctype>;
            grid_->entityImps_.resize(1);
        }

        /** \brief Constructor for a given grid object

        If you already have your grid object constructed you can
        hand it over using this constructor.

        If you construct your factory class using this constructor
        the pointer handed over to you by the method createGrid() is
        the one you supplied here.
         */
        GridFactoryBase(FoamGrid<dimgrid, dimworld, ctype>* grid)
        : factoryOwnsGrid_(false)
        {
            grid_ = grid;
            grid_->entityImps_.resize(1);
        }

        /** \brief Destructor */
        ~GridFactoryBase() override {
            if (grid_ && factoryOwnsGrid_)
                delete grid_;
        }

        /** \brief Insert a vertex into the coarse grid */
        void insertVertex(const FieldVector<ctype,dimworld>& pos) override {
            std::get<0>(grid_->entityImps_[0]).emplace_back(
              0,    // level
              pos,  // position
              grid_->getNextFreeId()
            );
            vertexArray_.push_back(&*std::get<0>(grid_->entityImps_[0]).rbegin());
        }

        /** \brief Obtain an element's insertion index
         */
        unsigned int
        insertionIndex(const typename FoamGrid<dimgrid, dimworld, ctype>::Traits::template Codim<0>::Entity &entity) const override
        {
            return entity.impl().target_->leafIndex_;
        }

        /** \brief Obtain a vertex' insertion index
         */
        unsigned int
        insertionIndex(const typename FoamGrid<dimgrid, dimworld, ctype>::Traits::template Codim<dimgrid>::Entity &vertex) const override
        {
            return vertex.impl().target_->leafIndex_;
        }

        /** \brief Obtain a boundary's insertion index
         */
        unsigned int
        insertionIndex(const typename FoamGrid<dimgrid, dimworld, ctype>::LeafIntersection& intersection) const override
        {
            return intersection.boundarySegmentIndex();
        }

    protected:

        // Pointer to the grid being built
        FoamGrid<dimgrid, dimworld, ctype>* grid_;

        // True if the factory allocated the grid itself, false if the
        // grid was handed over from the outside
        bool factoryOwnsGrid_;

        /** \brief Array containing all vertices */
        std::vector<FoamGridEntityImp<0, dimgrid, dimworld, ctype>*> vertexArray_;

        /** \brief Counter that creates the boundary segment indices */
        unsigned int boundarySegmentCounter_ = 0;
    };

template <int dimgrid, int dimworld, class ct>
    class GridFactory<FoamGrid<dimgrid, dimworld, ct> >
        : public GridFactoryBase<dimgrid, dimworld, ct>
    {};

/** \brief Specialization of GridFactoryBase for 1D-FoamGrid<1, dimworld>
    */
template <int dimworld, class ct>
    class GridFactory<FoamGrid<1, dimworld, ct> >
        : public GridFactoryBase<1, dimworld, ct>
    {
        /** \brief Grid dimension */
        enum {dimgrid = 1};

        /** \brief Type used by the grid for coordinates */
        using ctype = ct;

        /** \brief Alias for FoamGrid entities */
        template<int mydim>
        using EntityImp = FoamGridEntityImp<mydim, dimgrid, dimworld, ctype>;

	using GridPtrType = typename GridFactoryBase<1, dimworld, ct>::GridPtrType;

    public:

        GridFactory() {}

        GridFactory(FoamGrid<1, dimworld, ctype>* grid):
            GridFactoryBase<1,dimworld,ctype>(grid)
        {}

        /** \brief Insert a boundary segment.
            This is only needed if you want to control the numbering of the boundary segments
        */
        void insertBoundarySegment(const std::vector<unsigned int>& vertices) override
        {
            boundarySegmentIndices_[ vertices[0] ] = this->boundarySegmentCounter_++;
        }

        /** \brief Insert a boundary segment and the boundary segment geometry
            This influences the ordering of the boundary segments.
        */
        void insertBoundarySegment(const std::vector<unsigned int>& vertices,
                                   const std::shared_ptr<BoundarySegment<dimgrid, dimworld> >& boundarySegment) override
        {
            DUNE_THROW(Dune::NotImplemented, "Parameterized boundary segments are not implemented");
        }

        /** \brief Return true if leaf intersection was inserted as boundary segment
        */
        bool wasInserted( const typename FoamGrid<dimgrid, dimworld, ctype>::LeafIntersection &intersection ) const override
        {
          if ( !intersection.boundary() || intersection.inside().level() != 0 )
            return false;

          const auto& vertex = intersection.inside().template subEntity<1>(intersection.indexInInside());
          const auto& it = boundarySegmentIndices_.find( this->insertionIndex(vertex) );
          return (it != boundarySegmentIndices_.end());
        }

        /** \brief Insert an element into the coarse grid
            \param type The GeometryType of the new element
            \param vertices The vertices of the new element, using the DUNE numbering
        */
        void insertElement(const GeometryType& type,
                           const std::vector<unsigned int>& vertices) override
        {
            assert(type.isLine());
            // insert new element
            std::get<1>(this->grid_->entityImps_[0]).emplace_back(
                this->vertexArray_[vertices[0]],
                this->vertexArray_[vertices[1]],
                0, // level
                this->grid_->getNextFreeId()
            );
        }

// disable the deprecated interface when using Dune > 2.8
#if DUNE_VERSION_LTE(DUNE_COMMON, 2, 8)
DUNE_NO_DEPRECATED_BEGIN
        /** \brief Insert a parametrized element into the coarse grid
        \param type The GeometryType of the new element
        \param vertices The vertices of the new element, using the DUNE numbering
        \param elementParametrization A function prescribing the shape of this element
        */
        [[deprecated("Signature with VirtualFunction is deprecated and will be removed after Dune 2.8. Use signature with std::function.")]]
        void insertElement(const GeometryType& type,
                           const std::vector<unsigned int>& vertices,
                           const std::shared_ptr<VirtualFunction<FieldVector<ctype,dimgrid>,FieldVector<ctype,dimworld> > >& elementParametrization) override
        {
            assert(type.isLine());
            EntityImp<1> newElement(this->vertexArray_[vertices[0]],
                                    this->vertexArray_[vertices[1]],
                                    0, // level
                                    this->grid_->getNextFreeId());
            // save the pointer to the element parametrization
            newElement.elementParametrization_ =
              [elementParametrization](const FieldVector<ctype,dimgrid>& x){
                FieldVector<ctype,dimworld> y;
                elementParametrization->evaluate(x, y);
                return y;
              };

            std::get<dimgrid>(this->grid_->entityImps_[0]).push_back(newElement);
        }
DUNE_NO_DEPRECATED_END
#endif

        /** \brief Insert a parametrized element into the coarse grid
        \param type The GeometryType of the new element
        \param vertices The vertices of the new element, using the DUNE numbering
        \param elementParametrization A function prescribing the shape of this element
        */
        void insertElement(const GeometryType& type,
                           const std::vector<unsigned int>& vertices,
                           std::function<FieldVector<ctype,dimworld>(FieldVector<ctype,dimgrid>)> elementParametrization)
#if DUNE_VERSION_GT(DUNE_COMMON, 2, 7)
        override
#endif
        {
            assert(type.isLine());
            EntityImp<1> newElement(this->vertexArray_[vertices[0]],
                                    this->vertexArray_[vertices[1]],
                                    0, // level
                                    this->grid_->getNextFreeId());

            newElement.elementParametrization_ = elementParametrization;
            std::get<dimgrid>(this->grid_->entityImps_[0]).push_back(newElement);
        }

        /** \brief Finalize grid creation and hand over the grid

        The receiver takes responsibility of the memory allocated for the grid
        */
        GridPtrType createGrid() override
        {
            // Prevent a crash when this method is called twice in a row
            // You never know who may do this...
            if (this->grid_==nullptr)
                return nullptr;

            // make facets (vertices) know about the element
            for (auto& element : std::get<dimgrid>(this->grid_->entityImps_[0]))
            {
                element.vertex_[0]->elements_.push_back(&element);
                element.vertex_[1]->elements_.push_back(&element);
            }

            // Create the index sets
            this->grid_->setIndices();

            // ////////////////////////////////////////////////
            //   Set the boundary ids
            // ////////////////////////////////////////////////

            for (auto& facet : std::get<0>(this->grid_->entityImps_[0]))
            {
                if (facet.elements_.size() == 1) // if boundary facet
                {
                    const auto it = boundarySegmentIndices_.find( facet.vertex_[0]->leafIndex_ );
                    if (it != boundarySegmentIndices_.end())
                        facet.boundarySegmentIndex_ = it->second;
                    else
                        facet.boundarySegmentIndex_ = this->boundarySegmentCounter_++;
                }
            }

            // ////////////////////////////////////////////////
            //   Hand over the new grid
            // ////////////////////////////////////////////////

            Dune::FoamGrid<dimgrid, dimworld, ctype>* tmp = this->grid_;
            tmp->numBoundarySegments_ = this->boundarySegmentCounter_;
            this->grid_ = nullptr;
            return GridPtrType(tmp);
        }

      private:
        std::unordered_map<unsigned int, unsigned int> boundarySegmentIndices_;
    };

    /** \brief Specialization of GridFactoryBase for 2D-FoamGrid<2, dimworld>
    */
template <int dimworld, class ct>
    class GridFactory<FoamGrid<2, dimworld, ct> >
        : public GridFactoryBase<2, dimworld, ct>
    {
        /** \brief Grid dimension */
        enum {dimgrid = 2};

        /** \brief Type used by the grid for coordinates */
        using ctype = ct;

        /** \brief Alias for FoamGrid entities */
        template<int mydim>
        using EntityImp = FoamGridEntityImp<mydim, dimgrid, dimworld, ctype>;

	using GridPtrType = typename GridFactoryBase<2, dimworld, ct>::GridPtrType;

    public:

        GridFactory() {}

        GridFactory(FoamGrid<2, dimworld, ctype>* grid):
            GridFactoryBase<2, dimworld, ctype>(grid)
        {}

        /** \brief Insert a boundary segment.
            This is only needed if you want to control the numbering of the boundary segments
        */
        void insertBoundarySegment(const std::vector<unsigned int>& vertices) override
        {
            std::array<unsigned int, 2> vertexIndices {{ vertices[0], vertices[1] }};

            // sort the indices
            if ( vertexIndices[0] > vertexIndices[1] )
              std::swap( vertexIndices[0], vertexIndices[1] );

            boundarySegmentIndices_[ vertexIndices ] = this->boundarySegmentCounter_++;
        }

        /** \brief Insert a boundary segment (== a line) and the boundary segment geometry
            This influences the ordering of the boundary segments.
        */
        void insertBoundarySegment(const std::vector<unsigned int>& vertices,
                                   const std::shared_ptr<BoundarySegment<dimgrid, dimworld> >& boundarySegment) override
        {
            DUNE_THROW(Dune::NotImplemented, "Parameterized boundary segments are not implemented");
        }

        /** \brief Return true if leaf intersection was inserted as boundary segment
        */
        bool wasInserted( const typename FoamGrid<dimgrid, dimworld, ctype>::LeafIntersection &intersection ) const override
        {
          if ( !intersection.boundary() || intersection.inside().level() != 0 )
            return false;

          // Get the vertices of the intersection
          const auto refElement = ReferenceElements<ctype, dimgrid>::general(intersection.inside().type());

          const int subIdx0 = refElement.subEntity(intersection.indexInInside(), 1, /*idx*/0, dimgrid);
          const auto vertex0 = intersection.inside().template subEntity<2>( subIdx0 );
          const int subIdx1 = refElement.subEntity(intersection.indexInInside(), 1, /*idx*/1, dimgrid);
          const auto vertex1 = intersection.inside().template subEntity<2>( subIdx1 );

          std::array<unsigned int, 2> vertexIndices {{
            this->insertionIndex( vertex0 ),
            this->insertionIndex( vertex1 )
          }};

          // sort the indices
          if ( vertexIndices[0] > vertexIndices[1] )
            std::swap( vertexIndices[0], vertexIndices[1] );

          const auto& it = boundarySegmentIndices_.find( vertexIndices );
          return (it != boundarySegmentIndices_.end());
        }

        /** \brief Insert an element into the coarse grid
            \param type The GeometryType of the new element
            \param vertices The vertices of the new element, using the DUNE numbering
        */
        void insertElement(const GeometryType& type,
                           const std::vector<unsigned int>& vertices) override {

            assert(type.isTriangle());

            EntityImp<dimgrid> newElement(/*level=*/0, this->grid_->getNextFreeId());
            newElement.vertex_[0] = this->vertexArray_[vertices[0]];
            newElement.vertex_[1] = this->vertexArray_[vertices[1]];
            newElement.vertex_[2] = this->vertexArray_[vertices[2]];

            std::get<dimgrid>(this->grid_->entityImps_[0]).push_back(newElement);
        }

// disable the deprecated interface when using Dune > 2.8
#if DUNE_VERSION_LTE(DUNE_COMMON, 2, 8)
DUNE_NO_DEPRECATED_BEGIN
        /** \brief Insert a parametrized element into the coarse grid
        \param type The GeometryType of the new element
        \param vertices The vertices of the new element, using the DUNE numbering
        \param elementParametrization A function prescribing the shape of this element
        */
        [[deprecated("Signature with VirtualFunction is deprecated and will be removed after Dune 2.8. Use signature with std::function.")]]
        void insertElement(const GeometryType& type,
                           const std::vector<unsigned int>& vertices,
                           const std::shared_ptr<VirtualFunction<FieldVector<ctype,dimgrid>,FieldVector<ctype,dimworld> > >& elementParametrization) override
        {
            assert(type.isTriangle());
            EntityImp<dimgrid> newElement(/*level=*/0, this->grid_->getNextFreeId());
            newElement.vertex_[0] = this->vertexArray_[vertices[0]];
            newElement.vertex_[1] = this->vertexArray_[vertices[1]];
            newElement.vertex_[2] = this->vertexArray_[vertices[2]];
            // save the pointer to the element parametrization
            newElement.elementParametrization_ =
              [elementParametrization](const FieldVector<ctype,dimgrid>& x){
                FieldVector<ctype,dimworld> y;
                elementParametrization->evaluate(x, y);
                return y;
              };

            std::get<dimgrid>(this->grid_->entityImps_[0]).push_back(newElement);
        }
DUNE_NO_DEPRECATED_END
#endif

        /** \brief Insert a parametrized element into the coarse grid
        \param type The GeometryType of the new element
        \param vertices The vertices of the new element, using the DUNE numbering
        \param elementParametrization A function prescribing the shape of this element
        */
        void insertElement(const GeometryType& type,
                           const std::vector<unsigned int>& vertices,
                           std::function<FieldVector<ctype,dimworld>(FieldVector<ctype,dimgrid>)> elementParametrization)
#if DUNE_VERSION_GT(DUNE_COMMON, 2, 7)
        override
#endif
        {
            assert(type.isTriangle());
            EntityImp<dimgrid> newElement(/*level=*/0, this->grid_->getNextFreeId());
            newElement.vertex_[0] = this->vertexArray_[vertices[0]];
            newElement.vertex_[1] = this->vertexArray_[vertices[1]];
            newElement.vertex_[2] = this->vertexArray_[vertices[2]];
            newElement.elementParametrization_ = elementParametrization;
            std::get<dimgrid>(this->grid_->entityImps_[0]).push_back(newElement);
        }

        /** \brief Finalize grid creation and hand over the grid
        The receiver takes responsibility of the memory allocated for the grid
        */
        GridPtrType createGrid() override
        {
            // Prevent a crash when this method is called twice in a row
            // You never know who may do this...
            if (this->grid_==nullptr)
                return nullptr;

            // ////////////////////////////////////////////////////
            //   Create the edges
            // ////////////////////////////////////////////////////

            // For fast retrieval: a map from pairs of vertices to the edge that connects them
            std::unordered_map<std::pair<const EntityImp<0>*, const EntityImp<0>*>, EntityImp<1>*,
                               HashPair<const EntityImp<0>*>> edgeMap;
            for (auto& element : std::get<dimgrid>(this->grid_->entityImps_[0]))
            {
                const auto refElement = ReferenceElements<ctype, dimgrid>::general(element.type());

                // Loop over all edges of this element
                for (std::size_t i=0; i<element.facet_.size(); ++i)
                {
                    // get pointer to the edge (insert edge if it's not in the entity list yet)
                    auto edge = [&](){
                      // Get two vertices of the potential edge
                      auto v0 = element.vertex_[refElement.subEntity(i, 1, 0, 2)];
                      auto v1 = element.vertex_[refElement.subEntity(i, 1, 1, 2)];
                      // sort pointers
                      if ( v0 > v1 )
                        std::swap( v0, v1 );

                      // see if the edge was already inserted
                      // the pointer pair hash is symmetric so we don't have to check for pair(v1,v0)
                      auto e = edgeMap.find(std::make_pair(v0, v1));
                      if (e != edgeMap.end())
                        return e->second;

                      // edge was not inserted yet: insert the edge now
                      std::get<1>(this->grid_->entityImps_[0]).emplace_back(v0, v1, /*level=*/0, this->grid_->getNextFreeId());
                      // insert it into the map of inserted edges
                      auto newEdge = &*std::get<1>(this->grid_->entityImps_[0]).rbegin();
                      edgeMap.insert(std::make_pair(std::make_pair(v0, v1), newEdge));
                      return newEdge;
                    }();

                    // make element know about the edge
                    element.facet_[i] = edge;

                    // make edge know about the element
                    edge->elements_.push_back(&element);
                }
            }

            // Create the index sets
            this->grid_->setIndices();


            // ////////////////////////////////////////////////
            //   Set the boundary ids
            // ////////////////////////////////////////////////

            for (auto& facet : std::get<1>(this->grid_->entityImps_[0]))
            {
                if (facet.elements_.size() == 1) // if boundary facet
                {
                    std::array<unsigned int, 2> vertexIndices {{ facet.vertex_[0]->leafIndex_, facet.vertex_[1]->leafIndex_ }};
                    // sort the indices
                    if ( vertexIndices[0] > vertexIndices[1] )
                      std::swap( vertexIndices[0], vertexIndices[1] );

                    const auto it = boundarySegmentIndices_.find( vertexIndices );
                    if (it != boundarySegmentIndices_.end())
                        facet.boundarySegmentIndex_ = it->second;
                    else
                        facet.boundarySegmentIndex_ = this->boundarySegmentCounter_++;
                }
            }

            // ////////////////////////////////////////////////
            //   Hand over the new grid
            // ////////////////////////////////////////////////

            Dune::FoamGrid<dimgrid, dimworld, ctype>* tmp = this->grid_;
            tmp->numBoundarySegments_ = this->boundarySegmentCounter_;
            this->grid_ = nullptr;
            return GridPtrType(tmp);
        }

      private:
        template<class T, class U = T>
        struct HashPair {
          std::size_t operator() (const std::pair<T, U>& a) const {
            std::size_t seed = 0;
            hash_combine(seed, a.first);
            hash_combine(seed, a.second);
            return seed;
          }
        };

        struct HashUIntArray {
          std::size_t operator() (const std::array<unsigned int, 2>& a) const {
            return hash_range(a.begin(), a.end());
          }
        };

        std::unordered_map<std::array<unsigned int, 2>, unsigned int, HashUIntArray> boundarySegmentIndices_;
    };

} // end namespace Dune

#endif
