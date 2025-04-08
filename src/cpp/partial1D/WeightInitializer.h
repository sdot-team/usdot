#pragma once

#include <tl/support/memory/BumpPointerPool.h>
#include "System.h"

namespace usdot {
    
/**
*/
template<class TF,class Density>
class WeightInitializer {
public:
    using            Sys                 = System<TF,Density>;
    using            TI                  = std::size_t;
     
    /**/             WeightInitializer   ( Sys &sys );
     
    void             run                 (); ///<

private:
    struct Ag {
        TF  len_n() const { return end_n - beg_n; }

        TF  beg_u, end_u; ///< end position, and total length
        PI  beg_n, end_n; ///< cells indices
        TF  test_width;   ///< 
        Ag* prev_opt;     ///< prev Ag to optimize
        Ag* prev;
    };

    void             optimize_agglomerate( Ag *item );
    TF               inv_cdf             ( TF u ) const;
    TF               cdf                 ( TF x ) const;
     
    Vec<TF>          inv_cdf_values;
    Vec<TF>          dirac_masses;       ///< normalized
    TF               mul_coeff;          ///<
    Ag*              last_ag;
    BumpPointerPool  pool;
    Sys&             sys;
};


} // namespace usdot

#include "WeightInitializer.cxx"
