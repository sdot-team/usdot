#pragma once

#include "SymmetricBandMatrix.h"
#include "SimpleSolverInput.h"

namespace usdot {

/**
*/
template<class TF>
class SimpleSolver {
public:
    /**/                    SimpleSolver                    ( const SimpleSolverInput<TF> &input );

    void                    set_density_contrast            ( TF max_ratio ); ///< 
    PI                      nb_diracs                       () const;

    void                    solve                           (); ///< all in one

    // helpers
    TF                      normalized_dirac_convolution    ( TF normalized_pos, TF convolution_width = 2e-2 ) const;
    TF                      normalized_density_value        ( TF normalized_pos ) const;

    // directly modifiable inputs
    TF                      max_mass_ratio_error_target;
    bool                    throw_if_error;
    int                     verbosity;

    // outputs
    Vec<TF>                 max_mass_ratio_error_history;
    Vec<TF>                 norm_2_residual_history;
    Vec<TF>                 norm_2_rhs_history;

private:
    // diracs
    Vec<TF>                 sorted_dirac_positions;         ///<
    Vec<TF>                 sorted_dirac_weights;           ///<
    Vec<TF>                 sorted_dirac_masses;            ///<
    TF                      sum_of_dirac_masses;            ///<
    Vec<PI>                 sorted_dirac_nums;              ///<

    // density
    TF                      max_of_original_density_values; ///<
    Vec<TF>                 original_density_values;        ///<

    bool                    need_lower_contrast_ratio;      ///<
    TF                      dirac_contrast_ratio;           ///<
    Vec<TF>                 density_integrals;              ///<
    Vec<TF>                 density_values;                 ///<
    TF                      add_right;                      ///<
    TF                      add_left;                       ///<
    TF                      beg_x;                          ///<
    TF                      end_x;                          ///<

    // parameters
    TF                      global_mass_ratio;
};


} // namespace usdot

#include "SimpleSolver.cxx"
