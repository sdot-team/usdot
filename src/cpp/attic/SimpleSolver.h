#pragma once

#include "SymmetricBandMatrix.h"
#include "SimpleSolverInput.h"

namespace usdot {

/**
*/
template<class TF>
class SimpleSolver {
public:
    struct                  Errors                           { TF n2_residual, ni_residual; bool bad_cell = false; };

    /**/                    SimpleSolver                     ( const SimpleSolverInput<TF> &input );
 
    void                    set_density_contrast             ( TF max_ratio ); ///< 
    PI                      nb_diracs                        () const;
 
    TF                      normalized_dirac_convolution     ( TF normalized_pos, TF convolution_width = 2e-2 ) const;
     
    TF                      normalized_density_x_primitive   ( TF normalized_pos ) const;
    TF                      normalized_density_x_integral    ( TF x0, TF x1 ) const;
    TF                      normalized_density_primitive     ( TF normalized_pos ) const;
    TF                      normalized_density_integral      ( TF x0, TF x1 ) const;
    TF                      normalized_density_value         ( TF normalized_pos ) const;
 
    Vec<TF>                 normalized_cell_barycenters      () const;
    Vec<TF>                 normalized_cell_masses           () const; ///<
    
    auto                    newton_system_ap                 ( TF eps = 1e-6 ) -> std::tuple<SymmetricBandMatrix<TF>,Vec<TF>,Errors,bool>; ///<
    auto                    newton_system                    () -> std::tuple<SymmetricBandMatrix<TF>,Vec<TF>,Errors,bool>; ///<
    Vec<TF>                 newton_dir                       (); ///<
    Errors                  errors                           () const; ///<
 
    bool                    update_weights                   (); ///<
    bool                    converged                        ( const Errors &er ) const;
    void                    solve                            (); ///< all in one

    void                    for_each_cell_mt                 ( auto &&func ) const; ///< func( dirac_position, dirac_weight, ldist, rdist, b0, b1, num_thread )
    void                    for_each_cell                    ( auto &&func ) const; ///< func( dirac_position, dirac_weight, ldist, rdist, b0, b1 )
    void                    plot                             ( Str filename = "glot.py" ) const;

    // directly modifiable inputs
    TF                      max_mass_ratio_error_target      = 1e-4;
    TF                      current_contrast_ratio           = 1e-9;
    TF                      target_contrast_ratio            = 1e-9;
    Vec<TF>                 sorted_dirac_weights;            ///<
    Vec<TF>                 sorted_dirac_masses;             ///<
    bool                    throw_if_error;
    bool                    multithread;                     ///<
    int                     verbosity;

    // outputs
    Vec<TF>                 max_mass_ratio_error_history;
    Vec<TF>                 norm_2_residual_history;
    Vec<TF>                 norm_2_rhs_history;

private:
    enum {                  CUT                              = 0 };
    enum {                  BALL                             = 1 };
    enum {                  DENSITY                          = 2 };
 
 
    // diracs 
    Vec<TF>                 sorted_dirac_positions;          ///<
    Vec<PI>                 sorted_dirac_nums;               ///<
 
    TF                      sum_of_dirac_masses;             ///<
 
    // density 
    TF                      max_of_original_density_values;  ///<
    bool                    need_lower_contrast_ratio;       ///<
    TF                      min_density_value;               ///<
 
    Vec<TF>                 original_density_values;         ///<

    Vec<TF>                 normalized_density_x_primitives; ///<
    Vec<TF>                 normalized_density_primitives;   ///<
    Vec<TF>                 normalized_density_values;       ///<
    TF                      add_right;                       ///<
    TF                      add_left;                        ///<
    TF                      beg_x;                           ///<
    TF                      end_x;                           ///<

    // parameters
    TF                      global_mass_ratio;
};


} // namespace usdot

#include "SimpleSolver.cxx"
