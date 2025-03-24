#pragma once

#include <vector>

namespace usdot {

/**
    values on a regular grid. Interpolation order = 1
*/
template<class TF>
class GridDensity {
public:
    using                   Vec              = std::vector<TF>;
    using                   PI               = std::size_t;

    /**/                    GridDensity      ( const Vec &values, TF filter = 0, PI mul_x = 1, TF cut_ratio = 1e-6, TF target_mass = 1 );
    
    TF                      der_primitive    ( TF x ) const;
    TF                      der_integral     ( TF x0, TF x1 ) const;
    TF                      der_value        ( TF x ) const;

    TF                      x_primitive      ( TF x ) const;
    TF                      x_integral       ( TF x0, TF x1 ) const;

    TF                      primitive        ( TF x ) const;
    TF                      integral         ( TF x0, TF x1 ) const;
    TF                      value            ( TF x ) const;
    

    Vec                     der_primitives;  ///<
    Vec                     der_values;      ///<

    Vec                     x_primitives;    ///<
    Vec                     primitives;      ///<
    Vec                     values;          ///<
    PI                      offset;          ///<
};

} // namespace usdot

#include "GridDensity.cxx"
