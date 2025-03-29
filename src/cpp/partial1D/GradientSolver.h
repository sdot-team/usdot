#pragma once

#include <tl/support/containers/Vec.h>
#include "GradientSolverInput.h"
#include "ControledDensity.h"

namespace usdot {

/**
*/
template<class TF>
class GradientSolver {
public:
    using                   DensityPtr                       = std::unique_ptr<ControledDensity<TF>>;
    using                   PI                               = std::size_t;
    using                   TV                               = Vec<TF>;

    /**/                    GradientSolver                   ( GradientSolverInput<TF> &&input );
 
    int                     update_weights                   ();
    void                    solve                            ();
    void                    plot                             ( std::string filename = "glot.py" ) const;

    TV                      cell_barycenters                 () const;
    Vec<Vec<TF,2>>          cell_boundaries                  () const;
    Vec<TF>                 cell_masses                      () const;
    TF                      error                            () const;
    
    Vec<TF>                 dirac_positions                  () const;
    TF                      density_value                    ( TF pos, PI n = 0 ) const;
    Vec<Vec<TF>>            newton_dirs                      () const;
    PI                      nb_cells                         () const;

    TV                      normalized_cell_barycenters      () const;
    Vec<Vec<TF,2>>          normalized_cell_boundaries       () const;

    void                    for_each_normalized_system_item  ( auto &&func ) const;
    void                    for_each_normalized_cell         ( auto &&get_weight, auto &&func ) const;
    void                    for_each_normalized_cell         ( auto &&func ) const;

    // behavior
    bool                    throw_if_error;                  ///<
    bool                    multithread;                     ///<
    int                     verbosity;                       ///<
    
    // diracs 
    TV                      sorted_dirac_positions;          ///<
    TV                      sorted_dirac_weights;            ///<
    TV                      sorted_dirac_masses;             ///<
    Vec<PI>                 sorted_dirac_nums;               ///<
    
    // 
    TF                      sum_of_dirac_masses;             ///<
    TF                      global_mass_ratio;               ///<
    TF                      target_l2_error                  = 1e-5;
   
    // density
    TV                      original_density_values;         ///<
    TF                      beg_x_density;                   ///<
    TF                      end_x_density;                   ///<
    Vec<DensityPtr>         densities;                       ///<
};

} // namespace usdot

#include "GradientSolver.cxx" // IWYU pragma: export
