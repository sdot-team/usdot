#pragma once

#include "BndType.h"

namespace usdot {

template<class TF>
struct Interval {
    void    display( Displayer &ds ) const { DS_OBJECT( Interval, t0, t1, x0, x1 ); }
    TF      mid    () const { return ( x0 + x1 ) / 2; }

    BndType t0; ///< left cut type
    BndType t1; ///< right cut type
    TF      x0; 
    TF      x1; 
};

} // namespace usdot
