#pragma once

#include "Density.h"

namespace usdot {

/**
*/
template<class TF>
class ConvolutedLebesgue : public Density<TF> {
public:
    using        Convo             = Convolution<TF>;

    /**/         ConvolutedLebesgue( TF x0, TF x1, TF h, RcPtr<Convo> convolution );
       
    virtual auto convoluted        ( RcPtr<Convo> convolution ) const -> RcPtr<Density<TF>> override;
    virtual auto iterator          () const -> RcPtr<DensityIterator<TF>> override;
    virtual void display           ( Displayer &ds ) const override;
    virtual TF   min_x             ( TF eps = std::numeric_limits<TF>::epsilon() ) const override;
    virtual TF   max_x             ( TF eps = std::numeric_limits<TF>::epsilon() ) const override;
    virtual TF   mass              () const override;

    RcPtr<Convo> convolution;
    TF           x0;
    TF           x1;
    TF           h;
};

} // namespace usdot

#include "ConvolutedLebesgue.cxx"
