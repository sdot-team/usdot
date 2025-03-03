#pragma once

#include <tl/support/containers/Vec.h>
#include "Density.h"

namespace usdot {

/**
*/
template<class TF>
class PiecewiseConstantDensity : public Density<TF> {
public:
    /**/         PiecewiseConstantDensity( const Vec<TF> &xs, const Vec<TF> &ys ); ///< xs.size() must be equal to ys.size() + 1
       
    virtual auto cdf_approximation       ( TF epsilon ) const -> CdfApproximation<TF> override;
    virtual auto convoluted              ( RcPtr<Convolution<TF>> convolution ) const -> RcPtr<Density<TF>> override;
    virtual auto iterator                () const -> RcPtr<DensityIterator<TF>> override;
    virtual void display                 ( Displayer &ds ) const override;
    virtual TF   min_x                   ( TF eps = std::numeric_limits<TF>::epsilon() ) const override;
    virtual TF   max_x                   ( TF eps = std::numeric_limits<TF>::epsilon() ) const override;
    virtual TF   mass                    () const override;

    Vec<TF>      xs;
    Vec<TF>      ys;
};

} // namespace usdot

#include "PiecewiseConstantDensity.cxx"
