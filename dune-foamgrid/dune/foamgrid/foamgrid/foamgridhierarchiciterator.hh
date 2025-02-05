#ifndef DUNE_FOAMGRID_HIERARCHIC_ITERATOR_HH
#define DUNE_FOAMGRID_HIERARCHIC_ITERATOR_HH

/** \file
* \brief The FoamGridHierarchicIterator class
*/

#include <stack>

namespace Dune {


//**********************************************************************
//
/** \brief Iterator over the descendants of an entity.
* \ingroup FoamGrid
Mesh entities of codimension 0 ("elements") allow to visit all entities of
codimension 0 obtained through nested, hierarchic refinement of the entity.
Iteration over this set of entities is provided by the HierarchicIterator,
starting from a given entity.
*/
template<class GridImp>
class FoamGridHierarchicIterator
{
    static constexpr int dimworld = GridImp::dimensionworld;
    static constexpr int dimgrid = GridImp::dimension;

    using EntityImpPointer = const FoamGridEntityImp<
        dimgrid, dimgrid, dimworld, typename GridImp::ctype
    >*;

public:
    using Entity = typename GridImp::template Codim<0>::Entity;

    //! We only iterate over elements with this iterator
    enum { codimension = 0 };

    //! Constructor with element impl (begin iterator)
    FoamGridHierarchicIterator(EntityImpPointer target, int maxLevel)
    : maxLevel_(maxLevel)
    {
        // Load sons of target onto the iterator stack
        stackChildren_(target);

        // Set entity target to the next child if exists
        resetEntity_();
    }

    //! Constructor without valid element (end iterator)
    FoamGridHierarchicIterator(int maxLevel)
    : maxLevel_(maxLevel)
    {
        resetEntity_();
    }

    //! \todo Please doc me !
    void increment()
    {
        if (elemStack_.empty())
            return;

        auto target = elemStack_.top();
        elemStack_.pop();

        // Load sons of previous target onto the iterator stack
        stackChildren_(target);

        // Set entity target to the next stacked element if exists
        resetEntity_();
    }

    //! dereferencing
    const Entity& dereference() const { return virtualEntity_; }

    //! equality
    bool equals(const FoamGridHierarchicIterator<GridImp>& other) const
    { return virtualEntity_ == other.virtualEntity_; }

private:
    void stackChildren_(EntityImpPointer target)
    {
        // Load sons of target onto the iterator stack
        if (target->level() < maxLevel_ && !target->isLeaf())
            for (std::size_t i = 0; i < target->nSons(); i++)
                elemStack_.push(target->sons_[i]);
    }

    void resetEntity_()
    {
        virtualEntity_.impl().setToTarget(
            elemStack_.empty() ? nullptr : elemStack_.top()
        );
    }

    //! The entity that the iterator is pointing to
    Entity virtualEntity_;

    //! max level to iterate over
    int maxLevel_;

    /** \brief For depth-first search */
    std::stack<EntityImpPointer> elemStack_;
};


}  // end namespace Dune

#endif
