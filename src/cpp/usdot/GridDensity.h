#pragma once

#include "utility/TridiagonalPlus1LineSymmetricMatrix.h"
#include "utility/common_types.h"
#include <ostream>
#include <vector>

namespace usdot {

/**
    Continuous affine density.

    Values are normalized so that integral( -inf, +inf ) = 1    

    regular_positions = 0 if positions have no specific properties (at least sorted :) )
    regular_positions = 1 if positions are equally spaced 
    regular_positions = 2 if positions are 0, 1, 2, ... 

    BEWARE: `set_lag_ratio` must be called before any call to functions like primitive, value, ...

    lag_ratio = 0 -> use the original value
    lag_ratio = 1 -> flat curve
*/
template<class TF,int regular_positions=0>
class GridDensity {
public:
    using SY                  = TridiagonalPlus1LineSymmetricMatrix<TF>;
    using VF                  = std::vector<TF>;
    using VI                  = std::vector<PI>;
    using MF                  = std::vector<VF>;
   
    /**/  GridDensity         ( TF beg_original_positions, TF end_original_positions, const VF &original_values );
    /**/  GridDensity         ( const VF &original_positions, const VF &original_values );
    /**/  GridDensity         ( const VF &original_values );
   
    void  set_lag_ratio       ( TF t );

    TF    x_primitive         ( TF x ) const; ///< primitive( x d\rho )
    TF    x_integral          ( TF x0, TF x1 ) const; ///< integral( x d\rho, x0, x1 )
   
    TF    derivative          ( TF x ) const;
    TF    primitive           ( TF x ) const;
    TF    integral            ( TF x0, TF x1 ) const;
    TF    value               ( TF x ) const;

    TF    derivative          ( TF x, PI num_der_lag_ratio );
    TF    value               ( TF x, PI num_der_lag_ratio );

    TF    position            ( PI i ) const;
    TF    min_x               () const;
    TF    max_x               () const;
    TF    ptp_x               () const;
       
    void  plot                ( std::ostream &fs ) const;

    void  _append_der_value   ();
    void  _set_values         ( TF t );

    TF    beg_original_positions; ///<
    TF    end_original_positions; ///<
    VF    original_positions;     ///<
    VF    original_values;        ///<

    VF    x_primitives;           ///<
    VF    primitives;             ///<
    VF    positions;              ///<
    VF    values;                 ///<

    TF    current_lag_ratio;      ///<
    MF    der_values;             ///<
    SY    system;                 ///<
   
    TF    opt_pos_coeff;          ///<
    VI    opt_pos_begs;           ///<
    TF    opt_pos_beg;            ///<
    TF    opt_pos_end;            ///<
};

} // namespace usdot

#include "GridDensity.cxx"
