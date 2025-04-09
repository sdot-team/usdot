#pragma once

#include "ConvolutedLebesgue.h"

namespace usdot {


template<class TF,class Cl>
class ConvolutedLebesgueDensityIterator : public DensityIterator<TF> {
public:
    // static TF clamped( TF a ) {
    //     using namespace std;
    //     if ( isinf( a ) )
    //         return a < 0 ? numeric_limits<TF>::lowest() : numeric_limits<TF>::max();
    //     return a;
    // }

    virtual void barycenter( TF &pint, TF &area, TF l0, TF l1 ) const override {
        // pint += cl->h * cl->int_conv.px_primitive( l0 - cl->x1 );
        // pint -= cl->h * cl->int_conv.px_primitive( l0 - cl->x0 );
        // pint += cl->h * cl->int_conv.px_primitive( l1 - cl->x0 );
        // pint -= cl->h * cl->int_conv.px_primitive( l1 - cl->x1 );

        // area += cl->h * cl->int_conv.primitive( l0 - cl->x1 );
        // area -= cl->h * cl->int_conv.primitive( l0 - cl->x0 );
        // area += cl->h * cl->int_conv.primitive( l1 - cl->x0 );
        // area -= cl->h * cl->int_conv.primitive( l1 - cl->x1 );
        TODO;
    }
    
    virtual TF integral( TF l0, TF l1 ) const override {
        return cl->h * cl->convolution->step_integral( l0, l1, cl->x0, cl->x1 );
        // using namespace std;

        // const TF sqrt_pi = sqrt( /*pi*/ atan( TF( 1 ) ) * 4 );
        // const TF s = cl->std;
    
        // // compute a * b - c * d with diffs
        // //     ( a - c ) * ( b + d ) = ( a * b - c * d ) + a * d  - c * b
        // //     ( a + c ) * ( b - d ) = ( a * b - c * d ) - a * d  + c * b
        // auto dm = [&]( TF x, TF c0, TF c1 ) {
        //     const TF a = clamped( x - c0 );
        //     const TF c = clamped( x - c1 );
        //     const TF b = clamped( erf( a / s ) );
        //     const TF d = clamped( erf( c / s ) );
        //     return ( 
        //         + clamped( a + c ) * ( b - d )
        //         + ( c1 - c0 ) * ( b + d )
        //     ) / 2;
        // };
    
        // auto ga = [&]( TF x, TF c0, TF c1 ) {
        //     return (
        //         + exp( - clamped( pow( clamped( ( x - c0 ) / s ), 2 ) ) )
        //         - exp( - clamped( pow( clamped( ( x - c1 ) / s ), 2 ) ) )
        //     );
        // };
    
        // return cl->h / 2 * (
        //     + dm( l1, cl->x0, cl->x1 )
        //     + dm( l0, cl->x1, cl->x0 )
        //     + s / sqrt_pi * (
        //         + ga( l1, cl->x0, cl->x1 )
        //         + ga( l0, cl->x1, cl->x0 )
        //     )
        // );
    }
    
    virtual TF value( TF pos ) const override {
        return cl->h * cl->convolution->step_value( pos, cl->x0, cl->x1 );
    }

    virtual void display( Displayer &ds ) const override {
        DS_OBJECT( ConvolutedLebesgueDensityIterator, cl );
    }

    virtual bool move_backward() override {
        return false;
    }

    virtual bool move_forward() override {
        return false;
    }

    const Cl *cl;
};

#define DTP template<class TF>
#define UTP ConvolutedLebesgue<TF>

DTP UTP::ConvolutedLebesgue( TF x0, TF x1, TF h, RcPtr<Convo> convolution ) : convolution( convolution ), x0( x0 ), x1( x1 ), h( h ) {
    if ( x0 > x1 ) {
        std::swap( this->x0, this->x1 );
        this->h = - this->h;
    }
}

DTP RcPtr<Density<TF>> UTP::convoluted( RcPtr<Convo> convolution ) const {
    return new ConvolutedLebesgue<TF>{ x0, x1, h, this->convolution->convoluted( convolution ) };
}

DTP RcPtr<DensityIterator<TF>> UTP::iterator() const {
    if ( x0 == x1 )
        return {};

    auto *res = new ConvolutedLebesgueDensityIterator<TF,UTP>;
    const TF d = convolution->numerical_width();
    res->cl = this;
    if ( x0 > x1 ) {
        res->x0 = x1 - d;
        res->x1 = x0 + d;
    } else {
        res->x0 = x0 - d;
        res->x1 = x1 + d;
    }
    return res;
}

DTP void UTP::display( Displayer &ds ) const {
    DS_OBJECT( ConvolutedLebesgue, convolution, x0, x1, h );
}

DTP TF UTP::min_x( TF eps ) const {
    return x0 - convolution->numerical_width( eps );
}

DTP TF UTP::max_x( TF eps ) const {
    return x1 + convolution->numerical_width( eps );
}

DTP TF UTP::mass() const {
    return h * ( x1 - x0 );
}

#undef DTP
#undef UTP

} // namespace usdot
