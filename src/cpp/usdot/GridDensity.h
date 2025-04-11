#pragma once

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
*/
template<class TF,int regular_positions=0>
class GridDensity {
public:
    using VF                  = std::vector<TF>;
    using VI                  = std::vector<PI>;
   
    /**/  GridDensity         ( TF beg_positions, TF end_positions, const VF &values );
    /**/  GridDensity         ( const VF &positions, const VF &values );
    /**/  GridDensity         ( const VF &values );
   
    TF    x_primitive         ( TF x ) const; ///< primitive( x d\rho )
    TF    x_integral          ( TF x0, TF x1 ) const; ///< integral( x d\rho, x0, x1 )
   
    TF    primitive           ( TF x ) const;
    TF    integral            ( TF x0, TF x1 ) const;
    TF    value               ( TF x ) const;
   
    TF    position            ( PI i ) const;
    TF    min_x               () const;
    TF    max_x               () const;
    TF    ptp_x               () const;
       
    TF    inv_cdf_no_check    ( TF u ) const;
    TF    inv_cdf             ( TF u, TF ext_density ) const;
    TF    inv_cdf             ( TF u ) const;
    TF    cdf                 ( TF x, TF ext_density ) const;
    TF    cdf                 ( TF x ) const { return primitive( x ); }

    TF    der_inv_cdf_no_check( TF u ) const;
    TF    der_inv_cdf         ( TF u, TF ext_density ) const; ///< der inv cdf
   
    void  plot                ( std::ostream &fs ) const;
   
    VF    x_primitives;       ///<
    VF    primitives;         ///<
    VF    positions;          ///<
    VF    values;             ///<
   
    TF    opt_pos_coeff;      ///<
    VI    opt_pos_begs;       ///<
    TF    opt_pos_beg;        ///<
    TF    opt_pos_end;        ///<
   
    VI    opt_cdf_begs;       ///<
};

} // namespace usdot

#include "GridDensity.cxx"
