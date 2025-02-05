// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set ts=8 sw=4 et sts=4:
#ifndef DUNE_FOAMGRID_INTERSECTIONITERATORS_HH
#define DUNE_FOAMGRID_INTERSECTIONITERATORS_HH

#include <memory>
#include <dune/foamgrid/foamgrid/foamgridintersections.hh>
#include <dune/foamgrid/foamgrid/foamgridentity.hh>
#include <dune/foamgrid/foamgrid/foamgridvertex.hh>
#include <dune/foamgrid/foamgrid/foamgridnulliteratorfactory.hh>

/** \file
* \brief The FoamGridLeafIntersectionIterator and FoamGridLevelIntersectionIterator classes
*/

namespace Dune {

/** \brief Iterator over all element neighbors
* \ingroup FoamGrid
* Mesh entities of codimension 0 ("elements") allow to visit all neighbors, where
* a neighbor is an entity of codimension 0 which has a common entity of codimension 1
* These neighbors are accessed via a IntersectionIterator. This allows the implement
* non-matching meshes. The number of neighbors may be different from the number
* of corners of an element!
*/
template<class GridImp>
class FoamGridLeafIntersectionIterator
{

    enum {dimworld = GridImp::dimensionworld};
    enum {dimgrid  = GridImp::dimension};

    typedef typename GridImp::ctype ctype;

    typedef std::vector<typename std::vector<const FoamGridEntityImp<dimgrid, dimgrid, dimworld, ctype>*>::const_iterator> ElementVector;
    typedef typename ElementVector::const_iterator ElementVectorIterator;

    // Only the codim-0 entity is allowed to call the constructors
    friend class FoamGridEntity<0,dimgrid,GridImp>;

    template<typename, typename, typename>
    friend class Dune::IntersectionIterator;

    FoamGridLeafIntersectionIterator()
    {}

    //! Constructor for a given grid entity and a given neighbor
    FoamGridLeafIntersectionIterator(const FoamGridEntityImp<dimgrid, dimgrid, dimworld, ctype>* center, int facet)
        : intersection_(FoamGridLeafIntersection<GridImp>(center,facet)), leafNeighbors_(std::make_shared<ElementVector>())
    {
        if(facet==center->corners())
        {
            // This is the end iterator
            intersection_.impl().neighbor_=FoamGridNullIteratorFactory<dimgrid, dimworld, ctype>::null();
            return;
        }

        pushBackLeafNeighbours_(intersection_.impl().center_->facet_[intersection_.impl().facetIndex_], leafNeighbors_);
        if(leafNeighbors_->size()==1)
        {
            // This is a boundary facet.
            intersection_.impl().neighbor_ = FoamGridNullIteratorFactory<dimgrid, dimworld, ctype>::null();
            return;
        }

        // Search for the first intersection
        while(intersection_.impl().facetIndex_ != center->corners()) // not an  end iterator
        {
            leafNeighborIterator_ = leafNeighbors_->begin();
            intersection_.impl().neighbor_=*leafNeighborIterator_;
            while(leafNeighborIterator_!=leafNeighbors_->end() &&
                   intersection_.impl().center_==**leafNeighborIterator_)
            {
                ++leafNeighborIterator_;
                if(leafNeighborIterator_ != leafNeighbors_->end())
                  intersection_.impl().neighbor_=*leafNeighborIterator_;
            }
            if(leafNeighborIterator_==leafNeighbors_->end())
            {
                if(leafNeighbors_->size()==1)
                {
                    // This is a boundary intersection.
                    return;
                }else
                {
                    // No valid intersection found on this facet, move to next one.
                    ++intersection_.impl().facetIndex_;
                    if(intersection_.impl().facetIndex_ != center->corners())
                    {
                      leafNeighbors_->clear();
                      pushBackLeafNeighbours_(intersection_.impl().center_->facet_[intersection_.impl().facetIndex_], leafNeighbors_);
                    }
                }
            }else
                // intersection with another element found!
                break;
        }

        if(intersection_.impl().facetIndex_ == center->corners())
        {
            // This is an end iterator
            intersection_.impl().neighbor_=FoamGridNullIteratorFactory<dimgrid, dimworld, ctype>::null();
        }
    }

    /** \brief Constructor creating the 'one-after-last'-iterator */
    FoamGridLeafIntersectionIterator(const FoamGridEntityImp<dimgrid, dimgrid, dimworld, ctype>* center)
        : intersection_(FoamGridLeafIntersection<GridImp>(center,center->corners()))
    {
    }

    //! fill the element vector with all leaf neighbors. leafNeighbours is the element vector we need for the iterator.
    void pushBackLeafNeighbours_(const FoamGridEntityImp<dimgrid-1, dimgrid, dimworld, ctype>* facet, std::shared_ptr<ElementVector> leafNeighbours)
    {
      if (facet->isLeaf())
        for(auto it = facet->elements_.begin(); it != facet->elements_.end(); ++it)
          leafNeighbours->push_back(it);
      else
      {
        // do it recursively until the leaf level
        for(auto&& son : facet->sons_)
          pushBackLeafNeighbours_(son, leafNeighbours);
      }
    }

public:

    typedef Dune::Intersection<const GridImp, typename Dune::FoamGridLeafIntersection<GridImp> > Intersection;

    //! equality
    bool equals(const FoamGridLeafIntersectionIterator<GridImp>& other) const {
        return intersection_.impl().equals(other.intersection_.impl());
    }

    //! prefix increment
    void increment()
    {
        if(intersection_.impl().facetIndex_ == intersection_.impl().center_->corners())
        {
            // This is already the end iterator
            DUNE_THROW(InvalidStateException, "Cannot increment a one past the end iterator");
            return;
        }
        if(leafNeighbors_->size()==1)
        {
            // This was a boundary intersection go to the next facet
            ++intersection_.impl().facetIndex_;
            if(intersection_.impl().facetIndex_ < intersection_.impl().center_->corners()){
                // There is another facet, initialize neighbor_ iterator.
                leafNeighbors_->clear();
                pushBackLeafNeighbours_(intersection_.impl().center_->facet_[intersection_.impl().facetIndex_], leafNeighbors_);
                leafNeighborIterator_ = leafNeighbors_->begin();
                intersection_.impl().neighbor_=*leafNeighborIterator_;
            }
            else
            {
                // This is the end iterator
                intersection_.impl().neighbor_=FoamGridNullIteratorFactory<dimgrid, dimworld, ctype>::null();
                return;
            }
        }else
        {
            // Move to the next intersection of this facet
            ++leafNeighborIterator_;
            if(leafNeighborIterator_ != leafNeighbors_->end())
              intersection_.impl().neighbor_=*leafNeighborIterator_;
        }

        // Search for the first intersection with a neighbor
        while(intersection_.impl().facetIndex_ != intersection_.impl().center_->corners()) // still a valid facet
        {
            if(leafNeighbors_->size()==1)
            {
                // This is a boundary intersection.
                intersection_.impl().neighbor_=FoamGridNullIteratorFactory<dimgrid, dimworld, ctype>::null();
                return;
            }

            while(leafNeighborIterator_ != leafNeighbors_->end() &&
                   intersection_.impl().center_==**leafNeighborIterator_)
            {
                // Wrong level or neighbor points to center. In both cases this intersection is invalid.
                ++leafNeighborIterator_;
                if(leafNeighborIterator_ != leafNeighbors_->end())
                  intersection_.impl().neighbor_=*leafNeighborIterator_;
            }
            if(leafNeighborIterator_==leafNeighbors_->end())
            {
                if(leafNeighbors_->size()==1)
                {
                    // This is a boundary intersection.
                    intersection_.impl().neighbor_=FoamGridNullIteratorFactory<dimgrid, dimworld, ctype>::null();
                    return;
                }
                else
                {
                    // No valid intersection found on this facet, move to next facet.
                    ++intersection_.impl().facetIndex_;
                    if(intersection_.impl().facetIndex_ < intersection_.impl().center_->corners())
                    {
                        leafNeighbors_->clear();
                        pushBackLeafNeighbours_(intersection_.impl().center_->facet_[intersection_.impl().facetIndex_], leafNeighbors_);
                        // There is another facet, initialize neighbor_ iterator.
                        leafNeighborIterator_ = leafNeighbors_->begin();
                        intersection_.impl().neighbor_=*leafNeighborIterator_;
                    }
                }
            }else
                // intersection with another element found!
                break;
        }

        if(intersection_.impl().facetIndex_ == intersection_.impl().center_->corners())
        {
            // This is an end iterator
            intersection_.impl().neighbor_=FoamGridNullIteratorFactory<dimgrid, dimworld, ctype>::null();
        }
    }


    //! \brief dereferencing
    const Intersection & dereference() const {
        return intersection_;
    }

private:
    Intersection intersection_;
    //! The neighbor elements on the leaf level. Shared pointer to prevent copying when copying the iterator
    std::shared_ptr<ElementVector> leafNeighbors_;
    ElementVectorIterator leafNeighborIterator_;
};




//! \todo Please doc me !
template<class GridImp>
class FoamGridLevelIntersectionIterator
{

    enum { dimgrid  = GridImp::dimension };
    enum { dimworld = GridImp::dimensionworld };

    typedef typename GridImp::ctype ctype;

    // Only the codim-0 entity is allowed to call the constructors
    friend class FoamGridEntity<0, dimgrid, GridImp>;

    template<typename, typename, typename>
    friend class Dune::IntersectionIterator;

    FoamGridLevelIntersectionIterator()
    {}

    //! \brief Constructor for a given grid entity and a given neighbor
    //! \param center Pointer to the element where the iterator was created.
    //! \param facet The index of the facet to start the investigation.
    FoamGridLevelIntersectionIterator(const FoamGridEntityImp<dimgrid, dimgrid, dimworld, ctype>* center, int facet)
        : intersection_(FoamGridLevelIntersection<GridImp>(center,facet))
    {
        if(facet==center->corners())
        {
            // This is the end iterator
            intersection_.impl().neighbor_=FoamGridNullIteratorFactory<dimgrid, dimworld, ctype>::null();
            return;
        }

        if(center->facet_[intersection_.impl().facetIndex_]->elements_.size()==1)
        {
            // This is a boundary facet.
            intersection_.impl().neighbor_ =
              intersection_.impl().neighborEnd_=
              intersection_.impl().center_->facet_[intersection_.impl().facetIndex_]->elements_.end();
            return;
        }

        // Search for the first intersection with a same level neighbor
        while(intersection_.impl().facetIndex_ != center->corners()) // not an  end iterator
        {
            intersection_.impl().neighbor_=intersection_.impl().center_->facet_[intersection_.impl().facetIndex_]->elements_.begin();
            intersection_.impl().neighborEnd_=intersection_.impl().center_->facet_[intersection_.impl().facetIndex_]->elements_.end();
            while(intersection_.impl().neighbor_!=intersection_.impl().neighborEnd_ &&
                  (intersection_.impl().center_==*intersection_.impl().neighbor_
                   ||intersection_.impl().center_->level()!=(*intersection_.impl().neighbor_)->level()))
            {
                ++intersection_.impl().neighbor_;
            }
            if(intersection_.impl().neighbor_==intersection_.impl().neighborEnd_)
            {
                if(intersection_.impl().center_->facet_[intersection_.impl().facetIndex_]->elements_.size()==1)
                {
                    // This is a boundary intersection.
                    return;
                }else
                    // No valid intersection found on this facet, move to next one.
                    ++intersection_.impl().facetIndex_;
            }else
                // intersection with another element found!
                break;
        }

        if(intersection_.impl().facetIndex_ == center->corners())
        {
            // This is an end iterator
            intersection_.impl().neighbor_=FoamGridNullIteratorFactory<dimgrid, dimworld, ctype>::null();
        }

    }

    /** \brief Constructor creating the 'one-after-last'-iterator */
    FoamGridLevelIntersectionIterator(const FoamGridEntityImp<dimgrid, dimgrid, dimworld, ctype>* center)
        : intersection_(FoamGridLevelIntersection<GridImp>(center,center->corners()))
    {
    }

public:

    typedef Dune::Intersection<const GridImp, typename Dune::FoamGridLevelIntersection<GridImp> > Intersection;

  //! equality
  bool equals(const FoamGridLevelIntersectionIterator<GridImp>& other) const {
      return intersection_.impl().equals(other.intersection_.impl());
  }

    //! prefix increment
    void increment() {
        if(intersection_.impl().facetIndex_==
           intersection_.impl().center_->corners())
        {
            // This is already the end iterator
            DUNE_THROW(InvalidStateException, "Cannot increment a one past the end iterator");
            return;
        }
        if(intersection_.impl().center_->facet_[intersection_.impl().facetIndex_]->elements_.size()==1)
        {
            // This was a boundary intersection.
            ++intersection_.impl().facetIndex_;
            if(intersection_.impl().facetIndex_ < intersection_.impl().center_->corners()){
                // There is another facet, initialize neighbor_ iterator.
                intersection_.impl().neighbor_=intersection_.impl().center_->facet_[intersection_.impl().facetIndex_]->elements_.begin();
                intersection_.impl().neighborEnd_=intersection_.impl().center_->facet_[intersection_.impl().facetIndex_]->elements_.end();;
            }
            else
            {
                // This is the end iterator
                intersection_.impl().neighbor_=FoamGridNullIteratorFactory<dimgrid, dimworld, ctype>::null();
                return;
            }
        }else
        {
            // Move to the next intersection of this facet
            ++intersection_.impl().neighbor_;
        }

        // Search for the first intersection with a same level neighbor
        while(intersection_.impl().facetIndex_ != intersection_.impl().center_->corners()) // still a valid facet
        {
            if(intersection_.impl().center_->facet_[intersection_.impl().facetIndex_]->elements_.size()==1)
            {
                // This is a boundary intersection.
                intersection_.impl().neighbor_=intersection_.impl().neighborEnd_;
                return;
            }

            while(intersection_.impl().neighbor_!=intersection_.impl().neighborEnd_ &&
                  (intersection_.impl().center_==*intersection_.impl().neighbor_
                   ||intersection_.impl().center_->level()!=(*intersection_.impl().neighbor_)->level()))
            {
                // Wrong level or neighbor points to center. In both cases this intersection is invalid.
                ++intersection_.impl().neighbor_;
            }
            if(intersection_.impl().neighbor_==
               intersection_.impl().neighborEnd_)
            {
                if(intersection_.impl().center_->facet_[intersection_.impl().facetIndex_]->elements_.size()==1)
                {
                    // This is a boundary intersection.
                    return;
                }
                else
                {
                    // No valid intersection found on this facet, move to next facet.
                    ++intersection_.impl().facetIndex_;
                    if(intersection_.impl().facetIndex_ < intersection_.impl().center_->corners())
                    {
                        // There is another facet, initialize neighbor_ iterator.
                        intersection_.impl().neighbor_=intersection_.impl().center_->facet_[intersection_.impl().facetIndex_]->elements_.begin();
                        intersection_.impl().neighborEnd_=intersection_.impl().center_->facet_[intersection_.impl().facetIndex_]->elements_.end();
                    }
                }
            }else
                // intersection with another element found!
                break;
        }

        if(intersection_.impl().facetIndex_ == intersection_.impl().center_->corners())
        {
            // This is an end iterator
            intersection_.impl().neighbor_=FoamGridNullIteratorFactory<dimgrid, dimworld, ctype>::null();
        }
    }


    //! \brief dereferencing
    const Intersection & dereference() const
    {
        return intersection_;
    }
private:
    /** \brief The actual intersection */
    Intersection intersection_;
};


}  // namespace Dune

#endif
