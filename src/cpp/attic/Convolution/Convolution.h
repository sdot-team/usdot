#pragma once

#include <tl/support/memory/RcPtr.h>
#include <limits>

namespace usdot {

/**
*/
template<class TF>
class Convolution : public WithRefCount {
public:
    virtual     ~Convolution        () {}
    
    virtual TF   numerical_width    ( TF epsilon = std::numeric_limits<TF>::epsilon() ) const = 0;
    virtual auto convoluted         ( RcPtr<Convolution> that ) const -> RcPtr<Convolution> { TODO; return {}; }
    virtual void display            ( Displayer &ds ) const = 0;
      
    virtual TF   p0                 ( TF x ) const = 0; ///< value
    virtual TF   p1                 ( TF x ) const = 0; ///< primitive
    virtual TF   p2                 ( TF x ) const = 0; ///< primitive of the primitive

    // helper functions
    virtual TF   dirac_integral     ( TF l0, TF l1, TF dirac_pos ) const;
    virtual TF   step_integral      ( TF l0, TF l1, TF x0, TF x1 ) const;
    virtual TF   dirac_value        ( TF l, TF dirac_pos ) const;
    virtual TF   step_value         ( TF l, TF x0, TF x1 ) const;
};

} // namespace usdot

#include "Convolution.cxx"