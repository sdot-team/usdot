#pragma once

#include <tl/support/containers/Vec.h>
#include "ConvGridSolverInput.h"
#include "GridDensity.h"

namespace usdot {

/**
*/
template<class TF>
class ConvGridSolver {
public:
    using                   DensityPtr                       = std::unique_ptr<GridDensity<TF>>;
    struct                  Poly                             { TF c0, c1, c2, ca; void display( auto &ds) const { DS_OBJECT(Poly,c0,c1,c2,ca);} }; ///< poly to get optimal weight wrt radius of the first cell and blend ratio for target masses
    using                   PI                               = std::size_t;
    using                   TV                               = Vec<TF>;

    /**/                    ConvGridSolver                   ( ConvGridSolverInput<TF> &&input );
 
    void                    initialize_filter_value          ( TF filter_value ); ///<
    void                    go_to_filter_value               ( TF filter_value ); ///<
    void                    solve                            ();

    TV                      cell_barycenters                 () const;
    TV                      dirac_positions                  () const;
    TF                      density_value                    ( TF pos ) const;
    PI                      nb_diracs                        () const;
    void                    plot                             ( std::string filename = "glot.py" ) const;
    
    TV                      normalized_cell_barycenters      () const;
    auto                    normalized_cell_boundaries       () const -> std::vector<std::array<TF,2>>;
    TV                      normalized_cell_masses           () const; ///<
    TF                      normalized_error                 () const; ///< norm_2( log( current_mass / target_mass ) )

    void                    for_each_normalized_system_item  ( auto &&bad_cell, auto &&bb, auto &&bc, auto &&cb, auto &&cc ) const; ///< bb( PI index ), bc( PI index ), cb( PI index ), cc( PI index, TF m0, TF m1, TF v ) 
    void                    for_each_normalized_cell_mass    ( auto &&func ) const; ///< func( PI index, TF mass, bool bad_cell )
    
    void                    for_each_normalized_cell_mt      ( auto &&func ) const; ///< func( dirac_position, dirac_weight, ldist, rdist, b0, b1, num_thread )
    void                    for_each_normalized_cell         ( auto &&func ) const; ///< func( dirac_position, dirac_weight, ldist, rdist, b0, b1 )

    // directly modifiable inputs
    TF                      current_filter_value             = 1e-9;
    TF                      target_filter_value              = 1e-9;
    TF                      target_mass_error                = 1e-5;
    bool                    throw_if_error;                  ///<
    bool                    multithread;                     ///<
    int                     verbosity;                       ///<
    
    TV                      sorted_dirac_weights;            ///<
    TV                      sorted_dirac_masses;             ///<
    
private:
    void                    initialize_weights               ();
    int                     get_weights_for                  ( Vec<TF> &new_dirac_weights, const Vec<Vec<PI,2>> &connected_cells, const Vec<Poly> &polys, TF a );
    TF                      best_r_for_ib                    ( const Vec<Poly> &polys, PI n, TF a ) const;
    void                    get_system                       ( Vec<Vec<PI,2>> &connected_cells, Vec<Poly> &polys, TF &max_a, TF &cell_error, int &has_bad_cell, const TV &sorted_dirac_weight ) const;
    Vec<TF>                 ce_errors                        ( const Vec<Poly> &polys, PI nb, PI ne, TF a, TF r ) const;
    TF                      bi_error                         ( const Vec<Poly> &polys, PI n, TF a, TF r ) const;
    TF                      ii_error                         ( const Vec<Poly> &polys, PI n, TF a, TF r ) const;


    // diracs 
    TV                      sorted_dirac_positions;          ///<
    std::vector<PI>         sorted_dirac_nums;               ///<
   
    TF                      sum_of_dirac_masses;             ///<
   
    // density
    TV                      original_density_values;         ///<
    TF                      beg_x_density;                   ///<
    TF                      end_x_density;                   ///<
    DensityPtr              density;                         ///<

    // parameters
    TF                      global_mass_ratio;
};


} // namespace usdot

#include "ConvGridSolver.cxx"
