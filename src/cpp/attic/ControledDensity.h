#pragma once

#include <memory>
#include <vector>

namespace usdot {

/**
    values on a regular grid. Interpolation order = 1
*/
template<class TF>
class ControledDensity {
public:
    using Vec              = std::vector<TF>;
    using PI               = std::size_t;

    /**/  ControledDensity ( const Vec &values, TF ratio_at_end = 0, TF target_mass = -1 );

    TF    x_primitive      ( TF x ) const;
    TF    x_integral       ( TF x0, TF x1 ) const;

    TF    primitive        ( TF x ) const;
    TF    integral         ( TF x0, TF x1 ) const;
    TF    value            ( TF x ) const;
    
    TF    contrast         ( PI beg, PI end ) const;

    Vec   x_primitives;    ///<
    Vec   primitives;      ///<
    Vec   values;          ///<
    PI    offset;          ///<
};

} // namespace usdot

#include "ControledDensity.cxx"
