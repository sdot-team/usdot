#pragma once

#include "PiecewiseAffineDensity.h"
#include <tl/support/ASSERT.h>

template<class TF>
class PiecewiseAffineDensityIterator : public DensityIterator<TF> {
public:
    virtual void barycenter( TF &pint, TF &area, TF l0, TF l1 ) const override {
        TF a = ( l1 - l0 ) * cl->ys[ i ];
        pint += a * ( l1 + l0 ) / 2;
        area += a;
    }
    
    virtual TF integral( TF l0, TF l1 ) const override {
        return ( l1 - l0 ) * cl->ys[ i ];
    }

    virtual void display( Displayer &ds ) const override {
        DS_OBJECT( PiecewiseAffineDensityIterator, cl, this->x0, this->x1, i );
    }
    
    virtual TF value( TF pos ) const override {
        return cl->ys[ i ];
    }

    virtual bool move_backward() override {
        if ( i == 0 )
            return false;
        --i;

        this->x0 = cl->xs[ i + 0 ];
        this->x1 = cl->xs[ i + 1 ];
        return true;
    }

    virtual bool move_forward() override {
        if ( i + 1 >= cl->ys.size() )
            return false;
        ++i;

        this->x0 = cl->xs[ i + 0 ];
        this->x1 = cl->xs[ i + 1 ];
        return true;
    }

    const PiecewiseAffineDensity<TF> *cl;
    PI i;
};

#define DTP template<class TF>
#define UTP PiecewiseAffineDensity<TF>

DTP UTP::PiecewiseAffineDensity( const Vec<TF> &xs, const Vec<TF> &ys ) : xs( xs ), ys( ys ) {
    ASSERT( xs.size() == ys.size() + 1 );
}

DTP RcPtr<Density<TF>> UTP::convoluted( TF std ) const {
    //return new PiecewiseAffineDensity<TF>{ x0, x1, h, this->std + std };
    TODO;
    return const_cast<PiecewiseAffineDensity *>( this );
}

DTP RcPtr<DensityIterator<TF>> UTP::iterator() const {
    auto *res = new PiecewiseAffineDensityIterator<TF>;
    ASSERT( xs.size() >= 2 );
    res->x0 = xs[ 0 ];
    res->x1 = xs[ 1 ];
    res->cl = this;
    res->i = 0;
    return res;
}

DTP void UTP::display( Displayer &ds ) const {
    DS_OBJECT( PiecewiseAffineDensity, xs, ys );
}

DTP TF UTP::mass() const {
    TF res = 0;
    for( PI i = 0; i < ys.size(); ++i )
        res += ( xs[ i + 1 ] - xs[ i + 0 ] ) * ys[ i ];
    return res;
}

#undef DTP
#undef UTP
