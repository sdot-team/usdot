#pragma once

#include <tl/support/memory/BumpPointerPool.h>
#include "System.h"

namespace usdot {
    
/**
*/
template<class TF,class Density>
class WeightInitializer {
public:
    using            Sys              = System<TF,Density>;
    using            TI               = std::size_t;
  
    /**/             WeightInitializer( Sys &sys );
  
    void             run              (); ///<

private:
    struct Ag {
        TF  len_n() const { return end_n - beg_n; }

        TF  beg_x, end_x; ///< end position, and total length
        PI  beg_n, end_n; ///< cells indices
        Ag* prev;
    };

    TF               optimal_radius( TF pos, TF mass ) const;
    void             append_cell   ( TI n, TF beg_x, TF end_x );
     
    Ag*              last_ag;
    BumpPointerPool  pool;
    Sys&             sys;
};


} // namespace usdot

#include "WeightInitializer.cxx"
