#pragma once

#include "InvX2Convolution.h"

namespace usdot {
    
#define DTP template<class TF>
#define UTP InvX2Convolution<TF>

DTP void UTP::display( Displayer &ds ) const {
    DS_OBJECT( InvX2Convolution, width );
}

DTP TF UTP::numerical_width( TF epsilon ) const {
    using namespace std;
    const TF pi = 4 * atan( TF( 1 ) );
    return width * sqrt( 1 / ( width * pi * epsilon ) - 1 );
}

DTP TF UTP::p2( TF x ) const {
    using namespace std;
    const TF pi = 4 * atan( TF( 1 ) );
    return ( x * atan( x / width ) - width / 2 * log( width * width + x * x ) ) / pi;
}

DTP TF UTP::p1( TF x ) const {
    using namespace std;
    const TF pi = 4 * atan( TF( 1 ) );
    return atan( x / width ) / pi;
}

DTP TF UTP::p0( TF x ) const {
    using namespace std;
    const TF pi = 4 * atan( TF( 1 ) );
    return 1 / ( pow( x / width, 2 ) + 1 ) / ( width * pi );
}

#undef DTP
#undef UTP

} // namespace usdot
