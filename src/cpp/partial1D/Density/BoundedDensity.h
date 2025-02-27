#pragma once

#include "Density.h"

namespace usdot {

/**
*/
template<class TF>
class BoundedDensity : public Density<TF> {
public:
    /**/               BoundedDensity( const RcPtr<Density<TF>> &density, TF x0, TF x1 );
             
    virtual auto       convoluted    ( RcPtr<Convolution<TF>> convolution ) const -> RcPtr<Density<TF>> override;
    virtual auto       iterator      () const -> RcPtr<DensityIterator<TF>> override;
    virtual void       display       ( Displayer &ds ) const override;
    virtual TF         min_x         ( TF eps = std::numeric_limits<TF>::epsilon() ) const override;
    virtual TF         max_x         ( TF eps = std::numeric_limits<TF>::epsilon() ) const override;

    RcPtr<Density<TF>> density;
    TF                 x0;
    TF                 x1;
};

} // namespace usdot

#include "BoundedDensity.cxx"
