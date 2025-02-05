#ifndef DUNE_FOAMGRID_LEAFITERATOR_HH
#define DUNE_FOAMGRID_LEAFITERATOR_HH

/** \file
* \brief The FoamGridLeafIterator class
*/

namespace Dune {

/** \brief Iterator over all entities of a given codimension and level of a grid.
*  \ingroup FoamGrid
*/
template<int codim, PartitionIteratorType pitype, class GridImp>
class FoamGridLeafIterator
{
    enum {dimgrid  = GridImp::dimension};
    enum {dimworld = GridImp::dimensionworld};

public:

    using Entity = typename GridImp::template Codim<codim>::Entity;
    enum { codimension = codim };

    FoamGridLeafIterator(const GridImp& grid)
    : grid_(&grid)
    {

        /** \todo Can a make the fullRefineLevel work somehow? */
        const int fullRefineLevel = 0;

        const auto& entities = std::get<dimgrid-codim>(grid_->entityImps_[fullRefineLevel]);
        levelIterator_ = entities.begin();

        if (levelIterator_ != entities.end())
        {

          // The &* turns an iterator into a plain pointer
          virtualEntity_.impl().setToTarget(&*entities.begin());

          if (!virtualEntity_.impl().target_->isLeaf())
            increment();
        }

        // grid doesn't contain entities of this codimension
        else
            virtualEntity_.impl().setToTarget(nullptr);
    }

    //! Default constructor
    FoamGridLeafIterator()
    : grid_(nullptr)
    {
        virtualEntity_.impl().setToTarget(nullptr);
    }

    //! prefix increment
    void increment() {
        // Increment until you find a leaf entity
        do {
            globalIncrement();

        } while (levelIterator_!=std::get<dimgrid-codim>(grid_->entityImps_[grid_->maxLevel()]).end()
                 && !virtualEntity_.impl().target_->isLeaf());
    }

    //! dereferencing
    const Entity& dereference() const { return virtualEntity_; }

    //! equality
    bool equals(const FoamGridLeafIterator<codim, pitype, GridImp>& other) const {
      return virtualEntity_ == other.virtualEntity_;
    }

private:

    /** \brief This increment makes the iterator wander over all entities on all levels */
    void globalIncrement() {

        // Backup current level because it may not be accessible anymore after
        // setting the pointer to the next entity.
        const int oldLevel = virtualEntity_.level();

        // Increment on this level
        ++levelIterator_;
        virtualEntity_.impl().setToTarget(&(*levelIterator_));
        if (levelIterator_==std::get<dimgrid-codim>(grid_->entityImps_[oldLevel]).end())
            virtualEntity_.impl().setToTarget(nullptr);

        // If beyond the end of this level set to first of next level
        if (levelIterator_==std::get<dimgrid-codim>(grid_->entityImps_[oldLevel]).end() && oldLevel < grid_->maxLevel()) {

            const std::list<FoamGridEntityImp<dimgrid-codim, dimgrid, dimworld, typename GridImp::ctype> >& entities = std::get<dimgrid-codim>(grid_->entityImps_[oldLevel+1]);
            levelIterator_ = entities.begin();
            virtualEntity_.impl().setToTarget(&*entities.begin());

        }

    }

    // /////////////////////////////////////
    //   Data members
    // /////////////////////////////////////
    const GridImp* grid_;

    //! The entity that the iterator is pointing to
    Entity virtualEntity_;

    typename std::list<FoamGridEntityImp<dimgrid-codim, dimgrid, dimworld, typename GridImp::ctype> >::const_iterator levelIterator_;
};


}  // namespace Dune

#endif
