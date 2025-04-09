#pragma once

#include "LogGridSolverInput.h"
#include "TridiagonalMatrix.h"

namespace usdot {

/**
*/
template<class TF>
class LogGridSolver {
public:
    struct                  System                           { TridiagonalMatrix<TF> M; Vec<TF> V; };

    /**/                    LogGridSolver                    ( LogGridSolverInput<TF> &&input );
 
    void                    set_density_contrast             ( TF max_ratio ); ///<
    void                    update_weights                   ();
    auto                    newton_dir                       () const -> std::function<Vec<TF>( TF a )>;
    void                    solve                            ();

    TF                      normalized_error                 () const; ///< norm_2( log( current_mass / target_mass ) )
    Vec<TF>                 cell_barycenters                 () const;
    Vec<TF>                 dirac_positions                  () const;
    PI                      nb_diracs                        () const;
    void                    plot                             ( Str filename = "glot.py" ) const;

    TF                      normalized_density_x_integral_ap ( TF x0, TF x1, PI nb_steps = 100000 ) const;
    TF                      normalized_density_x_primitive   ( TF normalized_pos ) const;
    TF                      normalized_density_x_integral    ( TF x0, TF x1 ) const;
    TF                      normalized_density_primitive     ( TF normalized_pos ) const;
    TF                      normalized_density_integral      ( TF x0, TF x1 ) const;
    TF                      normalized_density_value         ( TF normalized_pos ) const;
 
    Vec<TF>                 normalized_cell_barycenters      () const;
    Vec<Vec<TF,2>>          normalized_cell_boundaries       () const;
    Vec<TF>                 normalized_cell_masses           () const; ///<

    void                    for_each_normalized_system_item  ( auto &&func ) const; ///< func( PI index, TF m0, TF m1, TF v, bool bad_cell )
    void                    for_each_normalized_cell_mass    ( auto &&func ) const; ///< func( PI index, TF mass, bool bad_cell )
    
    void                    for_each_normalized_cell_mt      ( auto &&func ) const; ///< func( dirac_position, dirac_weight, ldist, rdist, b0, b1, num_thread )
    void                    for_each_normalized_cell         ( auto &&func ) const; ///< func( dirac_position, dirac_weight, ldist, rdist, b0, b1 )

    // directly modifiable inputs
    TF                      current_contrast_ratio           = 1e-9;
    TF                      target_contrast_ratio            = 1e-9;
    TF                      target_l2_error                  = 1e-5;
    bool                    throw_if_error;                  ///<
    bool                    multithread;                     ///<
    int                     verbosity;                       ///<
    
    Vec<TF>                 sorted_dirac_weights;            ///<
    Vec<TF>                 sorted_dirac_masses;             ///<

private:
    // diracs 
    Vec<TF>                 sorted_dirac_positions;            ///<
    Vec<PI>                 sorted_dirac_nums;                 ///<
   
    TF                      sum_of_dirac_masses;               ///<
   
    // density   
    TF                      max_of_original_density_values;    ///<
    bool                    need_lower_contrast_ratio;         ///<
    Vec<TF>                 original_density_values;           ///<
    TF                      min_density_value;                 ///<
    
    
    TF                      normalized_density_gaussian_width; ///<
    Vec<TF>                 normalized_density_x_primitives;   ///<
    Vec<TF>                 normalized_density_primitives;     ///<
    Vec<TF>                 normalized_density_values;         ///<
    TF                      beg_x_density;                     ///<
    TF                      end_x_density;                     ///<

    // parameters
    TF                      global_mass_ratio;
};


} // namespace usdot

#include "LogGridSolver.cxx"
