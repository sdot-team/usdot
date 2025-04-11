#pragma once

#include "utility/BumpPointerPool.h"
#include "System.h"

namespace usdot {
    
/**
*/
template<class TF,class Density>
class WeightInitializer {
public:
    using            Sys                  = System<TF,Density>;
    using            VF                   = std::vector<TF>;
     
    /**/             WeightInitializer    ( Sys &sys );
     
    void             run                  (); ///<

    PI               nb_iterations;       ///<
    
private:
    struct Ag { ///< aggregate
        TF  len_n                         () const { return end_n - beg_n; }

        TF  beg_u, end_u;                 ///< end position, and total length
        PI  beg_n, end_n;                 ///< cells indices
        TF  test_width;                   ///< 
        Ag* prev_opt;                     ///< prev Ag to optimize
        Ag* prev;
    };
    
    void             merge_touching_aggregates(); ///< 
    void             make_isolated_aggregates ();
    void             optimize_aggregate       ( Ag *item, TF x_tol );
    void             set_the_weights          ();
    TF               der_inv_cdf              ( TF u ) const;
    TF               inv_cdf                  ( TF u ) const;
    TF               max_u                    () const { return inv_cdf_values.size() - 1; }
    TF               cdf                      ( TF x ) const;
         
    Ag*              last_ag_to_optimize;     ///<
    TF               coeff_ext_inv_cdf;       ///<
    VF               inv_cdf_values;          ///<     
    VF               dirac_masses;            ///< normalized
    TF               mul_coeff;               ///<
    Ag*              last_ag;
    BumpPointerPool  pool;
    Sys&             sys;
};


} // namespace usdot

#include "WeightInitializer.cxx"
