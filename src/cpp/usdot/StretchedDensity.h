#pragma once

#include "utility/TridiagonalPlus1LineSymmetricMatrix.h"
#include "utility/common_types.h"
#include <ostream>
#include <vector>

namespace usdot {

/**
    Continuous affine density.

*/
template<class TF>
class StretchedDensity {
public:
    using VF                      = std::vector<TF>;
    using VI                      = std::vector<PI>;
    using MF                      = std::vector<VF>;
       
    /**/  StretchedDensity        ( TF beg_original_positions, TF end_original_positions, const VF &original_values );
    /**/  StretchedDensity        ( const VF &original_positions, const VF &original_values );
    /**/  StretchedDensity        ( const VF &original_values );
       
    void  set_flattening_ratio    ( TF flattening_ratio );
    
    TF    x_primitive             ( TF x ) const; ///< primitive( x d\rho ). Works only for flattening_ratio == 0.
    TF    x_integral              ( TF x0, TF x1 ) const; ///< integral( x d\rho, x0, x1 ) Works only for flattening_ratio == 0.
       
    TF    primitive               ( TF x ) const;
    TF    integral                ( TF x0, TF x1 ) const;
    TF    value                   ( TF x ) const;
    
    TF    primitive               ( TF x, PI num_der_lag_ratio );
    TF    integral                ( TF x0, TF x1, PI num_der_lag_ratio );
    TF    value                   ( TF x, PI num_der_lag_ratio );
    
    TF    min_x                   () const;
    TF    max_x                   () const;
    TF    ptp_x                   () const;
           
    void  plot                    ( std::ostream &fs ) const;
    
    void  _append_der_value       ();
    VF    _primitive_of           ( const VF &positions, const VF &values ) const;
    void  _set_values             ( TF flattening_ratio );

    VF    original_x_primitives;  ///<
    VF    original_positions;     ///<
    VF    original_values;        ///<

    VF    flat_positions;         ///<
    TF    flat_value;             ///<

    TF    flattening_ratio;       ///<
    VF    primitives;             ///<
    VF    positions;              ///<
    VF    values;                 ///<

    MF    der_primitives;         ///<
    MF    der_values;             ///<
   
    VI    opt_pos_beg_inds;       ///<
    TF    opt_pos_beg_x;          ///<
    TF    opt_pos_mul_x;          ///<
};

} // namespace usdot

#include "StretchedDensity.cxx"
