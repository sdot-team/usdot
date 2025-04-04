#pragma once

#include <vector>

namespace usdot {

/**
    values on a regular grid. Interpolation order = 1
*/
template<class TF>
class GridDensity {
public:
    using VF               = std::vector<TF>;
    using TI               = std::size_t;

    /**/  GridDensity      ( VF &&values );

    TF    x_primitive      ( TF x ) const;
    TF    x_integral       ( TF x0, TF x1 ) const;

    TF    primitive        ( TF x ) const;
    TF    integral         ( TF x0, TF x1 ) const;
    TF    value            ( TF x ) const;

    TF    mass             () const;

private:
    VF    x_primitives;    ///<
    VF    primitives;      ///<
    VF    values;          ///<
};

} // namespace usdot

#include "GridDensity.cxx"
