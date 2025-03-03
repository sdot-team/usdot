#pragma once

#include "Density.h"

namespace usdot {

/**
    For each interval of xs, ys prodives
    - value at the beginning at the interval
    - slope for x = 0 at the beginning of the interval, x = 1 at the end of the interval
*/
template<class TF>
class PiecewiseAffineDensity : public Density<TF> {
public:
    /**/           PiecewiseAffineDensity( const Vec<TF> &xs, const Vec<Vec<TF,2>> &ys ); ///< xs.size() must be equal to ys.size() + 1
         
    virtual auto   cdf_approximation     ( TF epsilon ) const -> CdfApproximation<TF> override;
    virtual auto   convoluted            ( RcPtr<Convolution<TF>> convolution ) const -> RcPtr<Density<TF>> override;
    virtual auto   iterator              () const -> RcPtr<DensityIterator<TF>> override;
    virtual void   display               ( Displayer &ds ) const override;
    virtual TF     min_x                 ( TF eps = std::numeric_limits<TF>::epsilon() ) const override;
    virtual TF     max_x                 ( TF eps = std::numeric_limits<TF>::epsilon() ) const override;
    virtual TF     mass                  () const override;

    Vec<TF>        xs;
    Vec<Vec<TF,2>> ys;
};

} // namespace usdot

#include "PiecewiseAffineDensity.cxx"
