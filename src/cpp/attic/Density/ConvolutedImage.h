#pragma once

#include <tl/support/containers/Vec.h>
#include "ConvolutedLebesgue.h"

/**
*/
template<class TF,class IntConv=IntGauss<TF>>
class ConvolutedImage : public Density<TF> {
public:
    /**/         ConvolutedImage( Vec<TF> h, IntConv int_conv = {} );
       
    virtual auto convoluted     ( TF std ) const -> RcPtr<Density<TF>> override;
    virtual auto iterator       () const -> RcPtr<DensityIterator<TF>> override;
    virtual void display        ( Displayer &ds ) const override;
    virtual TF   mass           () const override;

    IntConv      int_conv;
    Vec<TF>      h;
};

#include "ConvolutedImage.cxx"
