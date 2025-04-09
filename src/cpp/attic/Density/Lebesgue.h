#pragma once

#include "Density.h"

namespace usdot {


/**
*/
template<class TF>
class Lebesgue : public Density<TF> {
public:
    /**/         Lebesgue         ( TF x0 = 0, TF x1 = 1, TF h = 1 );
       
    virtual auto cdf_approximation( TF epsilon ) const -> CdfApproximation<TF> override; ///< 
    virtual auto convoluted       ( RcPtr<Convolution<TF>> convolution ) const -> RcPtr<Density<TF>> override;
    virtual auto iterator         () const -> RcPtr<DensityIterator<TF>> override;
    virtual void display          ( Displayer &ds ) const override;
    virtual TF   min_x            ( TF eps = std::numeric_limits<TF>::epsilon() ) const override;
    virtual TF   max_x            ( TF eps = std::numeric_limits<TF>::epsilon() ) const override;
    virtual TF   mass             () const override;

    TF           x0;
    TF           x1;
    TF           h;
};


} // namespace usdot

#include "Lebesgue.cxx"
