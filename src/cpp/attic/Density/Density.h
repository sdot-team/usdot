#pragma once

#include "../Convolution/Convolution.h"
#include <tl/support/containers/Vec.h>
#include <tl/support/memory/RcPtr.h>
#include "../CdfApproximation.h"
#include "DensityIterator.h"
#include <limits>

namespace usdot {

/**
*/
template<class TF>
class Density : public WithRefCount {
public:
    virtual     ~Density          () {}

    virtual auto cdf_approximation( TF epsilon ) const -> CdfApproximation<TF>; ///< 
    virtual auto convoluted       ( RcPtr<Convolution<TF>> convolution ) const -> RcPtr<Density<TF>> = 0; 
    virtual auto iterator         () const -> RcPtr<DensityIterator<TF>> = 0; ///< fist sub-interval where it's "easy" to get integral of value. return nullptr if not possible to find such interval
    virtual TF   integral         ( TF l0, TF l1 ) const;
    virtual void display          ( Displayer &ds ) const = 0;
    virtual TF   value            ( TF x ) const;
    virtual TF   min_x            ( TF eps = std::numeric_limits<TF>::epsilon() ) const = 0;
    virtual TF   max_x            ( TF eps = std::numeric_limits<TF>::epsilon() ) const = 0;
    virtual TF   mass             () const;
};

} // namespace usdot

#include "Density.cxx"