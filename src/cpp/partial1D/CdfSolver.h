#pragma once

#include <tl/support/memory/BumpPointerPool.h>
#include "CdfApproximation.h"

namespace usdot {
    
/**
*/
template<class TF>
class CdfSolver {
public:
    /**/           CdfSolver    ( CdfApproximation<TF> &&cdf, const Vec<TF> &seed_coords, const Vec<TF> &mass_ratios );

    PI             nb_cells     () const { return seed_coords.size(); }
    void           solve        ( Vec<TF> &seed_weights ); ///<

private:
    struct Agglomeration {
        TF end_y() const { return beg_y + len_y; }
        TF len_n() const { return end_n - beg_n; }

        Agglomeration *prev;
        TF beg_y, len_y; ///< end position, and total length
        PI beg_n, end_n; ///< cells indices
    };

    void            append_cells( PI beg_n, PI end_n, TF beg_y, TF len_y );

    const Vec<TF>&  seed_coords;
    const Vec<TF>&  mass_ratios;
    Vec<TF>         cx;
    Vec<TF>         cy;
    Vec<TF>         cz;
    
    BumpPointerPool pool;
    Agglomeration*  last;
};


} // namespace usdot

#include "CdfSolver.cxx"
