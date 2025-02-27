#pragma once

#include "Convolution.h"

namespace usdot {
    
#define DTP template<class TF>
#define UTP Convolution<TF>

DTP TF UTP::dirac_integral( TF l0, TF l1, TF dirac_pos ) const {
    return p1( l1 - dirac_pos ) - p1( l0 - dirac_pos );
}

DTP TF UTP::step_integral( TF l0, TF l1, TF x0, TF x1 ) const {
    return p2( l1 - x0 ) - p2( l0 - x0 ) - p2( l1 - x1 ) + p2( l0 - x1 );
}

DTP TF UTP::dirac_value( TF l, TF dirac_pos ) const {
    return p0( l - dirac_pos );
}

DTP TF UTP::step_value( TF l, TF x0, TF x1 ) const {
    return p1( x1 - l ) - p1( x0 - l );
}

#undef DTP
#undef UTP

} // namespace usdot
