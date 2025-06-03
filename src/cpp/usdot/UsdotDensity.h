#pragma once

#include "utility/IndexedPositions.h"
#include <ostream>

namespace usdot {

/**
    Continuous affine density.

    Values are normalized so that integral( -inf, +inf ) = 1    

    BEWARE: `set_lag_ratio` must be called before any call to functions like primitive, value, ...

    lag_ratio = 0 -> use the original value
    lag_ratio = 1 -> flat curve
*/
template<class TF>
class UsdotDensity {
public:
    using IP                        = IndexedPositions<TF>;
    using VF                        = std::vector<TF>;
    using VI                        = std::vector<PI>;
         
    /**/  UsdotDensity              ( TF beg_positions, TF end_positions, const VF &values );
         
    TF    x2_integral               ( TF x0, TF x1, TF di = 0 ) const; ///< integral( (x-di)^2 d\rho, x0, x1 )
    TF    x_primitive               ( TF x ) const; ///< primitive( x d\rho )
    TF    x_integral                ( TF x0, TF x1 ) const; ///< integral( x d\rho, x0, x1 )
       
    TF    derivative                ( TF x ) const;
    TF    primitive                 ( TF x ) const;
    TF    integral                  ( TF x0, TF x1 ) const;
    TF    value                     ( TF x ) const;
      
    TF    inv_primitive             ( TF x ) const;
    TF    min_x                     () const;
    TF    max_x                     () const;
    TF    ptp_x                     () const;
             
    void  plot                      ( std::ostream &fs, std::string linestyle = "-", double linewidth = 1 ) const;
      
    void  _compute_values_for       ( TF flattening_ratio );
    void  _compute_primitives       ();
    VF    _primitive_of             ( const VF &values ) const;


    TF    beg_positions;            ///<
    TF    end_positions;            ///<
    VF    values;                   ///<
 
    IP    primitive_positions;      ///<
    VF    x_primitives;             ///<
    VF    primitives;               ///<
};

} // namespace usdot

#include "UsdotDensity.cxx"
