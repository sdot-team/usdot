#pragma once

#include "ConvGridSolverInput.h"
#include "GridDensity.h"

namespace usdot {

/**
*/
template<class TF>
class ConvGridSolver {
public:
    using                   DensityPtr                       = std::unique_ptr<GridDensity<TF>>;
    using                   Vec                              = std::vector<TF>;
    using                   PI                               = std::size_t;

    /**/                    ConvGridSolver                   ( ConvGridSolverInput<TF> &&input );
 
    void                    initialize_filter_value          ( TF filter_value ); ///<
    void                    go_to_filter_value               ( TF filter_value ); ///<
    int                     update_weights                   ();
    TF                      line_search                      ( const Vec &dir );
    std::pair<Vec,TF>       newton_dir                       () const;
    void                    solve                            ();

    Vec                     cell_barycenters                 () const;
    Vec                     dirac_positions                  () const;
    TF                      density_value                    ( TF pos ) const;
    PI                      nb_diracs                        () const;
    void                    plot                             ( std::string filename = "glot.py" ) const;
    
    Vec                     normalized_cell_barycenters      () const;
    auto                    normalized_cell_boundaries       () const -> std::vector<std::array<TF,2>>;
    Vec                     normalized_cell_masses           () const; ///<
    TF                      normalized_error                 () const; ///< norm_2( log( current_mass / target_mass ) )

    void                    for_each_normalized_system_item  ( auto &&func ) const; ///< func( PI index, TF m0, TF m1, TF v, bool bad_cell )
    void                    for_each_normalized_cell_mass    ( auto &&func ) const; ///< func( PI index, TF mass, bool bad_cell )
    
    void                    for_each_normalized_cell_mt      ( auto &&func ) const; ///< func( dirac_position, dirac_weight, ldist, rdist, b0, b1, num_thread )
    void                    for_each_normalized_cell         ( auto &&func ) const; ///< func( dirac_position, dirac_weight, ldist, rdist, b0, b1 )

    // directly modifiable inputs
    TF                      current_filter_value             = 1e-9;
    TF                      target_filter_value              = 1e-9;
    TF                      target_l2_error                  = 1e-5;
    bool                    throw_if_error;                  ///<
    bool                    multithread;                     ///<
    int                     verbosity;                       ///<
    
    Vec                     sorted_dirac_weights;            ///<
    Vec                     sorted_dirac_masses;             ///<

private:
    // diracs 
    Vec                     sorted_dirac_positions;            ///<
    std::vector<PI>         sorted_dirac_nums;                 ///<
   
    TF                      sum_of_dirac_masses;               ///<
   
    // density
    Vec                     original_density_values;           ///<
    TF                      beg_x_density;                     ///<
    TF                      end_x_density;                     ///<
    DensityPtr              density;                           ///<

    // parameters
    TF                      global_mass_ratio;
};


} // namespace usdot

#include "ConvGridSolver.cxx"
