#pragma once

#include "utility/TridiagonalSymmetricMatrix.h"
#include "utility/IndexedPositions.h"
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
    using IP                        = IndexedPositions<TF>;
    using VF                        = std::vector<TF>;
    using VI                        = std::vector<PI>;
    using MF                        = std::vector<VF>;
         
    /**/  DiffusionDensity          ( TF beg_original_positions, TF end_original_positions, const VF &original_values, TF start_flattening_ratio = 0 );
    /**/  DiffusionDensity          ( const VF &original_positions, const VF &original_values, TF start_flattening_ratio = 0 );
    /**/  DiffusionDensity          ( const VF &original_values, TF start_flattening_ratio = 0 );
         
    void  set_flattening_ratio      ( TF flattening_ratio );
    void  compute_derivatives       ( PI nb_derivatives );
      
    TF    x2_integral               ( TF x0, TF x1, TF di = 0 ) const; ///< integral( (x-di)^2 d\rho, x0, x1 )
    TF    x_primitive               ( TF x ) const; ///< primitive( x d\rho )
    TF    x_integral                ( TF x0, TF x1 ) const; ///< integral( x d\rho, x0, x1 )
       
    void  disk_integral             ( TF x0, TF x1, TF c, TF r, TF &mass, TF &dmass ) const;
    TF    disk_integral             ( TF x0, TF x1, TF c, TF r ) const;
    TF    derivative                ( TF x ) const;
    TF    primitive                 ( TF x ) const;
    TF    integral                  ( TF x0, TF x1 ) const;
    TF    value                     ( TF x ) const;
      
    TF    derivative                ( TF x, PI num_der_lag_ratio );
    TF    primitive                 ( TF x, PI num_der_lag_ratio ) const;
    TF    integral                  ( TF x0, TF x1, PI num_der_lag_ratio ) const;
    TF    value                     ( TF x, PI num_der_lag_ratio ) const;
      
    TF    min_x                     () const;
    TF    max_x                     () const;
    TF    ptp_x                     () const;
             
    void  plot                      ( std::ostream &fs, std::string linestyle = "-", double linewidth = 1 ) const;
      
    void  _compute_values_for       ( TF flattening_ratio );
    void  _compute_primitives       ();
    VF    _primitive_of             ( const VF &values ) const;
  
    VF    extended_positions;       ///<
    IP    indexed_positions;        ///<
    VF    original_values;          ///<
    TF    max_value;                ///<
    VF    values;                   ///<
  
    VF    x_primitives;             ///<
    VF    primitives;               ///<

    TF    current_flattening_ratio; ///<
    TF    coeff_flattening_ratio;   ///<
    MF    der_primitives;           ///<
    MF    der_values;               ///<
    TF    sys_div;                  ///<
    SY    sys_mat;                  ///<
    VF    sys_vec;                  ///<

    bool  normalization = 1;        ///<
};

} // namespace usdot

#include "DiffusionDensity.cxx"
