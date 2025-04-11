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
        TF  total_mass;                   ///< 
        Ag* prev_opt;                     ///< prev Ag to optimize
        Ag* prev;
    };
      
    void             optimize_the_new_aggregates();
    void             merge_touching_aggregates  (); ///< 
    void             make_isolated_aggregates   ();
    void             set_the_weights            ();
           
    Ag*              last_ag_to_optimize;       ///<
    TF               coeff_ext_density;         ///<
    Ag*              last_ag;
    BumpPointerPool  pool;
    Sys&             sys;
};


} // namespace usdot

#include "WeightInitializer.cxx"
