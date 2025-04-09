#pragma once

#include <tl/support/containers/Vec.h>
#include "../Convolution/Convolution.h"
#include "Density.h"

namespace usdot {

/**

*/
template<class TF>
class ConvolutedDiracs : public Density<TF> {
public:
    using        Convo           = Convolution<TF>;

    /**/         ConvolutedDiracs( Vec<TF> positions, Vec<TF> weights, RcPtr<Convo> convolution );
       
    virtual auto convoluted      ( RcPtr<Convolution<TF>> convolution ) const -> RcPtr<Density<TF>> override;
    virtual auto iterator        () const -> RcPtr<DensityIterator<TF>> override;
    virtual void display         ( Displayer &ds ) const override;
    virtual TF   min_x           ( TF eps = std::numeric_limits<TF>::epsilon() ) const override;
    virtual TF   max_x           ( TF eps = std::numeric_limits<TF>::epsilon() ) const override;
    virtual TF   mass            () const override;

    RcPtr<Convo> convolution;
    Vec<TF>      positions;
    Vec<TF>      weights;
};

} // namespace usdot

#include "ConvolutedDiracs.cxx"
