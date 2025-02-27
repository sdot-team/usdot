#pragma once

#include "Convolution.h"

namespace usdot {

/**
    1 / ( ( x / width )**2 + 1 ) / ( width * pi )
*/
template<class TF>
class InvX2Convolution : public Convolution<TF> {
public:
    /* */        InvX2Convolution( TF width ) : width( width ) {}

    virtual TF   numerical_width    ( TF epsilon = std::numeric_limits<TF>::epsilon() ) const override;
    virtual void display            ( Displayer &ds ) const override;

    virtual TF   p0                 ( TF x ) const override;
    virtual TF   p1                 ( TF x ) const override;
    virtual TF   p2                 ( TF x ) const override;

    TF           width;
};

} // namespace usdot

#include "InvX2Convolution.cxx"
