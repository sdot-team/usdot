#pragma once

#include <tl/support/containers/Vec.h>
#include <tl/support/Displayer.h>

namespace usdot {

template<class TF>
struct CdfApproximation {
    Vec<TF> xs; ///< stop coordinates
    Vec<TF> ys; ///< height for each x
    Vec<TF> zs; ///< additional value at the middle of each interval
};

} // namespace usdot
