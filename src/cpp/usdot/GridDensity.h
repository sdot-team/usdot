#pragma once

#include "utility/common_types.h"
#include <ostream>
#include <vector>

namespace usdot {

/**
    values on a regular grid. Interpolation order = 1
*/
template<class TF>
class GridDensity {
public:
    using VF               = std::vector<TF>;

    /**/  GridDensity      ( VF &&values );

    TF    x_primitive      ( TF x ) const;
    TF    x_integral       ( TF x0, TF x1 ) const;

    TF    primitive        ( TF x ) const;
    TF    integral         ( TF x0, TF x1 ) const;
    TF    value            ( TF x ) const;

    TF    width            () const;
    TF    mass             () const;

    void  get_inv_cdf      ( VF &inv_cdf_values, TF &mul_coeff, PI nb_bins ) const;
    void  plot             ( std::ostream &fs ) const;

private:
    VF    x_primitives;    ///<
    VF    primitives;      ///<
    VF    values;          ///<
};

} // namespace usdot

#include "GridDensity.cxx"
