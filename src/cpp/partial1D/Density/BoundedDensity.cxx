#pragma once

#include <tl/support/ASSERT.h>
#include "BoundedDensity.h"

namespace usdot {

template<class TF,class Cl>
class BoundedDensityDensityIterator : public DensityIterator<TF> {
public:
    virtual void barycenter( TF &pint, TF &area, TF l0, TF l1 ) const override {
        iterator->barycenter( pint, area, l0, l1 );
    }
    
    virtual TF integral( TF l0, TF l1 ) const override {
        return iterator->integral( l0, l1 );
    }

    virtual void display( Displayer &ds ) const override {
        DS_OBJECT( BoundedDensityDensityIterator, iterator, min_x, max_x, this->x0, this->x1 );
    }
    
    virtual TF value( TF pos ) const override {
        return iterator->value( pos );
    }

    virtual bool move_backward() override {
        // if there's nothing intersecting before the current iterator position
        if ( iterator->x0 <= min_x )
            return false;

        //
        if ( ! iterator->move_backward() )
            return false;

        // if we've gone too far
        if ( iterator->x1 <= min_x ) {
            iterator->move_forward();
            return false;
        }

        //
        this->x0 = max( min_x, iterator->x0 ); 
        this->x1 = min( max_x, iterator->x1 );
        return true;
    }

    virtual bool move_forward() override {
        // if there's nothing intersecting after the current iterator position
        if ( iterator->x1 >= max_x )
            return false;

        //
        if ( ! iterator->move_forward() )
            return false;

        // if we've gone too far
        if ( iterator->x0 >= max_x ) {
            iterator->move_backward();
            return false;
        }

        //
        this->x0 = max( min_x, iterator->x0 ); 
        this->x1 = min( max_x, iterator->x1 );
        return true;
    };

    RcPtr<DensityIterator<TF>> iterator;
    TF min_x; 
    TF max_x;
};

#define DTP template<class TF>
#define UTP BoundedDensity<TF>

DTP UTP::BoundedDensity( const RcPtr<Density<TF>> &density, TF x0, TF x1 ) : density( density ), x0( x0 ), x1( x1 ) {
    ASSERT( x0 <= x1 );
}

DTP RcPtr<Density<TF>> UTP::convoluted( RcPtr<Convolution<TF>> convolution ) const {
    TODO;
    return const_cast<BoundedDensity *>( this );
}

DTP RcPtr<DensityIterator<TF>> UTP::iterator() const {
    using namespace std;

    // look if possible to find an intersecting iterator
    auto iterator = density->iterator();
    if ( ! iterator )
        return {};

    while ( iterator->x1 <= x0 )
        if ( ! iterator->move_forward() )
            return {};

    if ( iterator->x0 >= x1 )
        return {};

    //
    auto res = RcPtr<BoundedDensityDensityIterator<TF,UTP>>::New();
    res->x0 = max( x0, iterator->x0 ); 
    res->x1 = min( x1, iterator->x1 );
    res->iterator = iterator;
    res->min_x = x0; 
    res->max_x = x1;

    return res;
}

DTP void UTP::display( Displayer &ds ) const {
    DS_OBJECT( BoundedDensity, density, x0, x1 );
}

DTP TF UTP::min_x( TF eps ) const {
    return x0;
}

DTP TF UTP::max_x( TF eps ) const {
    return x1;
}

#undef DTP
#undef UTP

} // namespace usdot
