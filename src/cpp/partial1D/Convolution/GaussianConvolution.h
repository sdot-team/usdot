#pragma once

#include "Convolution.h"

namespace usdot {

/**
    exp( - ( x / width )**2 ) * 2 / ( sqrt( pi ) * width )
*/
template<class TF>
class GaussianConvolution : public Convolution<TF> {
public:
    /* */        GaussianConvolution( TF width ) : width( width ) {}

    virtual TF   numerical_width    ( TF epsilon = std::numeric_limits<TF>::epsilon() ) const override;
    virtual void display            ( Displayer &ds ) const override;
  
    virtual TF   p0                 ( TF x ) const override;
    virtual TF   p1                 ( TF x ) const override;
    virtual TF   p2                 ( TF x ) const override;

    // virtual TF   dirac_integral     ( TF l0, TF l1, TF dirac_pos ) const override;
    // virtual TF   step_integral      ( TF l0, TF l1, TF x0, TF x1 ) const override;
    // virtual TF   dirac_value        ( TF l, TF dirac_pos ) const override;
    // virtual TF   step_value         ( TF l, TF x0, TF x1 ) const override;

    TF           width;
};

} // namespace usdot

#include "GaussianConvolution.cxx"
