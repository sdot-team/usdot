#pragma once

#include "utility/TridiagonalSymmetricMatrix.h"
#include "utility/common_types.h"
#include <ostream>
#include <vector>

namespace usdot {

/**
    Continuous affine density.

    Values are normalized so that integral( -inf, +inf ) = 1    

    BEWARE: `set_lag_ratio` must be called before any call to functions like primitive, value, ...

    lag_ratio = 0 -> use the original value
    lag_ratio = 1 -> flat curve
*/
template<class TF>
class DiffusionDensity {
public:
    using SY                        = TridiagonalSymmetricMatrix<TF>;
    using VF                        = std::vector<TF>;
    using VI                        = std::vector<PI>;
    using MF                        = std::vector<VF>;
         
    /**/  DiffusionDensity          ( TF beg_original_positions, TF end_original_positions, const VF &original_values );
    /**/  DiffusionDensity          ( const VF &original_positions, const VF &original_values );
    /**/  DiffusionDensity          ( const VF &original_values );
         
    void  set_flattening_ratio      ( TF flattening_ratio );
    void  compute_derivatives       ( PI nb_derivatives );
      
    TF    x_primitive               ( TF x ) const; ///< primitive( x d\rho )
    TF    x_integral                ( TF x0, TF x1 ) const; ///< integral( x d\rho, x0, x1 )
       
    TF    derivative                ( TF x ) const;
    TF    primitive                 ( TF x ) const;
    TF    integral                  ( TF x0, TF x1 ) const;
    TF    value                     ( TF x ) const;
      
    TF    derivative                ( TF x, PI num_der_lag_ratio );
    TF    primitive                 ( TF x, PI num_der_lag_ratio );
    TF    integral                  ( TF x0, TF x1, PI num_der_lag_ratio );
    TF    value                     ( TF x, PI num_der_lag_ratio );
      
    TF    min_x                     () const;
    TF    max_x                     () const;
    TF    ptp_x                     () const;
             
    void  plot                      ( std::ostream &fs ) const;
      
    void  _compute_values_for       ( TF flattening_ratio );
    VF    _primitive_of             ( const VF &values ) const;
  
    VF    original_values;          ///<
    VF    positions;                ///<
    VF    values;                   ///<
  
    VF    x_primitives;             ///<
    VF    primitives;               ///<

    TF    current_flattening_ratio; ///<
    MF    der_primitives;           ///<
    MF    der_values;               ///<
    SY    system;                   ///<
     
    VI    opt_pos_beg_inds;         ///<
    TF    opt_pos_beg_x;            ///<
    TF    opt_pos_mul_x;            ///<
};

} // namespace usdot

#include "DiffusionDensity.cxx"
