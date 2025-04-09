#pragma once

#include "PowerDiagram.h"

struct GluedCells {
    using            TF               = PowerDiagram::TF;
    using            PI               = TL_NAMESPACE::PI;
                         
    PI               nb_cells         () const { return end_ind - beg_ind; }
    void             display          ( Displayer &ds ) const;

    PI               beg_density_ind;  ///<

    PI               beg_ind;          ///< dirac index 
    PI               end_ind;          ///< dirac index 

    TF               beg_p;            ///<

    PI               num_moving_phase; ///<
    GluedCells*      next_gc_to_test;  ///<
    TF               off_update;       ///<
    GluedCells*      prev_gc;          ///<
    GluedCells*      next_gc;          ///<
};

inline void GluedCells::display( Displayer &ds ) const {
    DS_OJBECT( beg_density_ind, beg_ind, end_ind, beg_p );
}
