#pragma once

// #include <tl/support/P.h>
#include "Density.h"

namespace usdot {

#define DTP template<class TF>
#define UTP Density<TF>

DTP TF UTP::integral( TF l0, TF l1 ) const {
    RcPtr<DensityIterator<TF>> p = iterator();
    if ( ! p )
        return 0;

    while ( p->x1 <= l0 )
        if ( ! p->move_forward() )
            return 0;
    
    TF res = 0;
    while ( p->x0 < l1 ) {
        using namespace std;
        res += p->integral( max( l0, p->x0 ), min( l1, p->x1 ) );
        if ( ! p->move_forward() )
            break;
    }
    return res;
}

DTP TF UTP::value( TF x ) const {
    RcPtr<DensityIterator<TF>> p = iterator();
    if ( ! p )
        return 0;

    while ( p->x1 <= x )
        if ( ! p->move_forward() )
            return 0;
    
    if ( p->x0 > x )
        return 0;
    
    return p->value( x );
}

DTP TF UTP::mass() const { 
    RcPtr<DensityIterator<TF>> p = iterator();
    if ( ! p )
        return 0;

    TF res = 0; 
    do {
        res += p->integral( p->x0, p->x1 ); 
    } while ( p->move_forward() );
    return res;
}

#undef DTP
#undef UTP

} // namespace usdot
