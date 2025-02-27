#pragma once

#include "SumOfDensities.h"
#include <limits>

template<class TF>
class SumOfDensitiesDensityIterator : public DensityIterator<TF> {
public:
    virtual void barycenter( TF &pint, TF &area, TF l0, TF l1 ) const override {
        for( const RcPtr<DensityIterator<TF>> &iterator : iterators )
            if ( iterator->x0 < l1 && iterator->x1 > l0 )
                iterator->barycenter( pint, area, max( l0, iterator->x0 ), min( l1, iterator->x1 ) );
    }
    
    virtual TF integral( TF l0, TF l1 ) const override {
        TF res = 0;
        for( const RcPtr<DensityIterator<TF>> &iterator : iterators )
            if ( iterator->x0 < l1 && iterator->x1 > l0 )
                res += iterator->integral( max( l0, iterator->x0 ), min( l1, iterator->x1 ) );
        return res;
    }

    virtual TF value( TF pos ) const override {
        TF res = 0;
        for( const RcPtr<DensityIterator<TF>> &iterator : iterators )
            if ( iterator->x0 <= pos && iterator->x1 > pos )
                res += iterator->value( pos );
        return res;
    }

    virtual void display( Displayer &ds ) const override {
        DS_OBJECT( SumOfDensitiesDensityIterator, this->x0, this->x1, iterators );
    }
    
    virtual bool move_backward() override {
        using namespace std;

        // try find the previous "event" (i.e. the previous bnd of interval < this->x0)
        TF new_x0 = numeric_limits<TF>::lowest();
        for( const RcPtr<DensityIterator<TF>> &iterator : iterators ) {
            if ( iterator->x0 < this->x0 && iterator->x0 > new_x0 )
                new_x0 = iterator->x0;
            if ( iterator->x1 < this->x0 && iterator->x1 > new_x0 )
                new_x0 = iterator->x1;
        }

        // no new event -> impossible to move forward
        if ( new_x0 == this->x0 )
            return false;

        // else, use the new interval
        this->x1 = this->x0;
        this->x0 = new_x0;
        return true;
    }

    virtual bool move_forward() override {
        using namespace std;

        // try find the next "event" (i.e. the next bnd of interval > this->x1)
        TF new_x1 = numeric_limits<TF>::max();
        for( const RcPtr<DensityIterator<TF>> &iterator : iterators ) {
            if ( iterator->x0 > this->x1 && iterator->x0 < new_x1 )
                new_x1 = iterator->x0;
            if ( iterator->x1 > this->x1 && iterator->x1 < new_x1 )
                new_x1 = iterator->x1;
        }

        // no new event -> impossible to move forward
        if ( new_x1 == this->x1 )
            return false;

        // else, use the new interval
        this->x0 = this->x1;
        this->x1 = new_x1;
        return true;
    }

    Vec<RcPtr<DensityIterator<TF>>> iterators; ///< beware: these iterators may be "valid" but without any intersection with *this
};

#define DTP template<class TF>
#define UTP SumOfDensities<TF>

DTP UTP::SumOfDensities( const SubDensities &sub_densities ) : sub_densities( sub_densities ) {
}

DTP RcPtr<Density<TF>> UTP::convoluted( TF std ) const {
    SubDensities ns;
    for( const RcPtr<Density<TF>> &sd : sub_densities )
        ns << sd->convoluted( std );
    return new SumOfDensities<TF>( ns );
}

DTP RcPtr<DensityIterator<TF>> UTP::iterator() const {
    using namespace std;

    // try to make a SumOfDensitiesDensityIterator
    auto res = RcPtr<SumOfDensitiesDensityIterator<TF>>::New();

    // get valid iterators + x0 (first event)
    res->x0 = numeric_limits<TF>::max();
    for( const RcPtr<Density<TF>> &sd : sub_densities ) {
        if ( RcPtr<DensityIterator<TF>> iterator = sd->iterator() ) {
            res->x0 = min( res->x0, iterator->x0 );
            res->iterators << iterator;
        }
    }

    // no valid sub-iterator ?
    if ( res->iterators.empty() )
        return {};

    // get x1
    res->x1 = numeric_limits<TF>::max();
    for( const RcPtr<DensityIterator<TF>> &iterator : res->iterators ) {
        if ( iterator->x0 > res->x0 )
            res->x1 = min( res->x1, iterator->x0 );
        res->x1 = min( res->x1, iterator->x1 );
    }

    return res;
}

DTP void UTP::display( Displayer &ds ) const {
    DS_OBJECT( SumOfDensities, sub_densities );
}

DTP TF UTP::mass() const {
    TF res = 0;
    for( const RcPtr<Density<TF>> &sd : sub_densities )
        res += sd->mass();
    return res;
}

#undef DTP
#undef UTP
