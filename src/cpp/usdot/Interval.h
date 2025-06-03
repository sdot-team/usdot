#pragma once

#include "BndType.h"

namespace usdot {

template<class TF>
struct Interval {
    bool     is_disk() const { return ( PI8( t0 ) & ( 4 * 3 ) ) == 0; }
    
    #ifdef   TL_DISPLAYER_IS_DEFINED
    void     display( Displayer &ds ) const { ds.start_array(); ds << x0 << x1 << t0 << t1; ds.end_array(); }
    #endif

    BndType  t0;    ///< left cut type
    BndType  t1;    ///< right cut type
    TF       x0; 
    TF       x1; 
};

} // namespace usdot
