#pragma once

#include <tl/support/containers/Vec.h>
#include "Density.h"

namespace usdot {

/**
*/
template<class TF>
class SumOfDensities : public Density<TF> {
public:
    using        SubDensities     = Vec<std::pair<TF,RcPtr<Density<TF>>>>;
   
    /**/         SumOfDensities   ( const SubDensities &sub_densities = {} );

    virtual auto cdf_approximation( TF epsilon ) const -> CdfApproximation<TF> override;
    virtual auto convoluted       ( RcPtr<Convolution<TF>> convolution ) const -> RcPtr<Density<TF>> override;
    virtual auto iterator         () const -> RcPtr<DensityIterator<TF>> override;
    virtual void display          ( Displayer &ds ) const override;
    virtual TF   min_x            ( TF eps = std::numeric_limits<TF>::epsilon() ) const override;
    virtual TF   max_x            ( TF eps = std::numeric_limits<TF>::epsilon() ) const override;
    virtual TF   mass             () const override;
  
    SubDensities sub_densities;   ///<
};

} // namespace usdot

#include "SumOfDensities.cxx"
