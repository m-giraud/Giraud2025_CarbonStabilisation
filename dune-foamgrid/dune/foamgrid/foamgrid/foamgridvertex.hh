// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set ts=8 sw=4 et sts=4:
#ifndef DUNE_FOAMGRID_VERTEX_HH
#define DUNE_FOAMGRID_VERTEX_HH

#include <dune/common/fvector.hh>
#include <dune/geometry/type.hh>
#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/exceptions.hh>


namespace Dune {

    /** \brief Base class for FoamGrid entity implementation classes */
    class FoamGridEntityBase
    {
    public:
        FoamGridEntityBase(int level, unsigned int id)
            : level_(level), id_(id), willVanish_()
        {}

        int level() const {
            return level_;
        }

        //! level
        int level_;

        //! entity number
        unsigned int levelIndex_;

        unsigned int leafIndex_;

        unsigned int id_;
        //! \brief Whether this entity will vanish due to coarsening.
        bool willVanish_;
    };

    /**
     * \brief The actual entity implementation
     *
     * \tparam dimentity The dimension of this entity
     * \tparam dimgrid The dimension of the grid
     * \tparam dimworld The world diemnsion
     */
    template <int dimentity, int dimgrid, int dimworld, class ctype>
    class FoamGridEntityImp {};

    /** \brief Vertex specialization of FoamGridEntityImp */
    template <int dimgrid, int dimworld, class ctype>
    class FoamGridEntityImp<0, dimgrid, dimworld, ctype>
        : public FoamGridEntityBase
    {
    public:

        FoamGridEntityImp(int level, const FieldVector<ctype, dimworld>& pos, unsigned int id)
            : FoamGridEntityBase(level, id), pos_(pos), nSons_(0)
            , elements_(), vertex_{{this}}, father_(nullptr), isNew_(false)
            , growthInsertionIndex_(-1)
        {
            sons_[0] = nullptr;
        }

        bool isLeaf() const {
            return sons_[0]==nullptr;
        }

        GeometryType type() const {
            return GeometryTypes::vertex;
        }

        bool hasFather() const
        {
            return father_!=nullptr;
        }

        //! This has no function yet in Foamgrid
        unsigned int boundarySegmentIndex() const {
            return boundarySegmentIndex_;
        }

        //! This has no function yet in Foamgrid
        unsigned int boundaryId() const {
            return boundaryId_;
        }

        /** \brief Number of corners (==1) */
        int corners() const {
            return 1;
        }

        FieldVector<ctype, dimworld> corner(int i) const {
            assert(i<this->corners());
            return pos_;
        }

        PartitionType partitionType() const {
            return InteriorEntity;
        }

        /** \brief Return level index of sub entity with codim = cc and local number i
         */
        int subLevelIndex (int i, unsigned int codim) const {
            assert(codim==dimgrid);
            return this->levelIndex_;
            DUNE_THROW(GridError, "Non-existing codimension requested!");
        }

        /** \brief Return leaf index of sub entity with codim = cc and local number i
         */
        int subLeafIndex (int i, unsigned int codim) const {
            assert(codim==dimgrid);
            return this->leafIndex_;
            DUNE_THROW(GridError, "Non-existing codimension requested!");
        }

        //! Position vector of this vertex
        FieldVector<ctype, dimworld> pos_;

         //! The number of refined vertices */
        unsigned int nSons_;

        //! Elements the vertex is related to
        std::vector<const FoamGridEntityImp<dimgrid, dimgrid ,dimworld, ctype>*> elements_;

        //! A vertex array for compatibility reasons with edges. Initialized with the this pointer.
        std::array<const FoamGridEntityImp<0, dimgrid ,dimworld, ctype>*, 1> vertex_;

        //! Boundary index if vertex is on boundary
        //  only used if the vertex is a boundary vertex
        unsigned int boundarySegmentIndex_;
        unsigned int boundaryId_;

        //! Pointer to father vertex on next coarser grid */
        FoamGridEntityImp<0, dimgrid, dimworld, ctype>* father_;

        //! Son vertex on the next finer grid
        std::array<FoamGridEntityImp<0, dimgrid, dimworld, ctype>*, 1> sons_;

        //! If the vertex was newly inserted (at run-time)
        bool isNew_;

        /** \brief If this vertex was created in a growth step this will be the index of insertion
         *         So if this is the first vertex added to the growth queue by calling insertVertex the index is 0
         *         The index will be valid until postGrow is called.
         *         To all other times this shall be -1
         */
        int growthInsertionIndex_;
    };

}

#endif
