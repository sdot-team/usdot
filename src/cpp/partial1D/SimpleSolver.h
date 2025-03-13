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

    TF                      normalized_dirac_convolution    ( TF normalized_pos, TF convolution_width = 2e-2 ) const;

    TF                      normalized_density_primitive    ( TF normalized_pos ) const;
    TF                      normalized_density_integral     ( TF x0, TF x1 ) const;
    TF                      normalized_density_value        ( TF normalized_pos ) const;

    Vec<TF>                 normalized_cell_masses          (); ///<
    TF                      normalized_error                (); ///<

    auto                    newton_system_ap                ( TF eps = 1e-6 ) -> std::tuple<SymmetricBandMatrix<TF>,Vec<TF>,TF,bool>; ///<
    auto                    newton_system                   () -> std::tuple<SymmetricBandMatrix<TF>,Vec<TF>,TF,bool>; ///<
    Vec<TF>                 newton_dir                      (); ///<

    bool                    update_weights                  (); ///<
    void                    solve                           (); ///< all in one

    // directly modifiable inputs
    TF                      max_mass_ratio_error_target;
    bool                    throw_if_error;
    int                     verbosity;

    // outputs
    Vec<TF>                 max_mass_ratio_error_history;
    Vec<TF>                 norm_2_residual_history;
    Vec<TF>                 norm_2_rhs_history;

private:
    enum {                  CUT                             = 0 };
    enum {                  BALL                            = 1 };
    enum {                  DENSITY                         = 2 };

    void                    for_each_cell_mt                ( auto &&func ) const; ///< func( dirac_position, dirac_weight, rdist, b0, b1, num_thread )
    void                    for_each_cell                   ( auto &&func ) const; ///< func( dirac_position, dirac_weight, rdist, b0, b1 )

    // diracs
    Vec<TF>                 sorted_dirac_positions;         ///<
    Vec<TF>                 sorted_dirac_weights;           ///<
    Vec<TF>                 sorted_dirac_masses;            ///<
    Vec<PI>                 sorted_dirac_nums;              ///<

    TF                      sum_of_dirac_masses;            ///<

    // density
    TF                      max_of_original_density_values; ///<
    bool                    need_lower_contrast_ratio;      ///<
    TF                      min_density_value;              ///<

    Vec<TF>                 original_density_values;        ///<

    Vec<TF>                 normalized_density_integrals;   ///<
    Vec<TF>                 normalized_density_values;      ///<
    TF                      add_right;                      ///<
    TF                      add_left;                       ///<
    TF                      beg_x;                          ///<
    TF                      end_x;                          ///<

    // parameters
    TF                      global_mass_ratio;
};


} // namespace usdot

#include "SimpleSolver.cxx"
