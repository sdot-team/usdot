#pragma once

#include "GaussianConvolution.h"

namespace usdot {
    
#define DTP template<class TF>
#define UTP GaussianConvolution<TF>

DTP void UTP::display( Displayer &ds ) const {
    DS_OBJECT( GaussianConvolution, width );
}

DTP TF UTP::numerical_width( TF epsilon ) const {
    using namespace std;
    const TF sqrt_pi = sqrt( 4 * atan( TF( 1 ) ) );
    return width * sqrt( - log( sqrt_pi * epsilon / 2 * width ) );
}

DTP TF UTP::p0( TF x ) const {
    using namespace std;
    const TF sqrt_pi = sqrt( 4 * atan( TF( 1 ) ) );
    return exp( - pow( x / width, 2 ) ) / ( width * sqrt_pi );
}

DTP TF UTP::p1( TF x ) const {
    using namespace std;
    return erf( x / width ) / 2;
}

DTP TF UTP::p2( TF x ) const {
    using namespace std;
    const TF sqrt_pi = sqrt( 4 * atan( TF( 1 ) ) );
    return ( width / sqrt_pi * exp( - pow( x / width, 2 ) ) + x * erf( x / width ) ) / 2;
}

// DTP TF UTP::step_integral( TF l0, TF l1, TF x0, TF x1 ) const {
//     using namespace std;
//     const TF sqrt_pi = sqrt( 4 * atan( TF( 1 ) ) );
//     auto primitive = [&]( TF x ) {
//         return width / sqrt_pi * exp( pow( x / width, 2 ) ) + x * erf( x / width );
//     };
//     return ( 
//         + primitive( l1 - x0 ) - primitive( l0 - x0 )
//         - primitive( l1 - x1 ) + primitive( l0 - x1 )
//     ) / 2;
// }

// DTP TF UTP::dirac_value( TF l, TF dirac_pos ) const {
//     using namespace std;
//     const TF sqrt_pi = sqrt( 4 * atan( TF( 1 ) ) );
// }

// DTP TF UTP::step_value( TF l, TF x0, TF x1 ) const {
//     using namespace std;
//     return ( erf( ( x1 - l ) / width ) - erf( ( x0 - l ) / width ) )/ 2;
// }


#undef DTP
#undef UTP

} // namespace usdot
