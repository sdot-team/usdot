#pragma once

#include "DiffusionDensity.h"
#include "Interval.h"
#include <limits>

namespace usdot {

/**
*/
template<class TF,class Density=DiffusionDensity<TF>>
class SdotSystem {
public:
    using              VB                         = std::vector<std::array<TF,2>>;
    using              VC                         = std::vector<Interval<TF>>;
    using              VS                         = std::vector<VC>;
    using              VF                         = std::vector<TF>;
    using              VI                         = std::vector<PI>;
  
    /**/               SdotSystem                 ( Density *density, const VF &dirac_positions, TF global_mass_ratio = 1, TF sep_equ_coords = 1e-6 );
  
    //   
    void               set_original_weights       ( const VF &weights ) const; ///<
    VF                 get_original_weights       () const; ///<
  
    void               solve_with_diffusion       ();
    void               initialize_weights         ();
    int                newton_iterations          ( const TF min_relax = 1e-6 );
    void               solve                      ();

    // visitors  
    T_T int            for_each_cell_polyline     ( T &&func, PI nb_divs = 100 ) const; ///< return a value != 0 is if void of negative cell
    T_T int            for_each_interval          ( T &&func ) const; ///< return a value != 0 is if void of negative cell
    void               plot                       ( Str filename = "glot.py" ) const;
  
    // basic information  
    PI                 nb_original_diracs         () const { return sorted_dirac_num_values.size(); }
    PI                 nb_sorted_diracs           () const { return sorted_dirac_positions.size(); }
    
    #ifdef TL_DISPLAYER_IS_DEFINED
    void               display                    ( Displayer &ds ) const;
    #endif
       
    // computations  
    int                get_sorted_newton_system   ( TridiagonalSymmetricMatrix<TF> &M, VF &V, TF &max_error_ratio ) const;
    VF                 sorted_cell_barycenters    () const;
    VB                 sorted_cell_boundaries     () const;
    VF                 sorted_cell_masses         () const;
    VS                 sorted_cells               () const; ///< intervals for each cell
      
    VF                 original_cell_barycenters  () const;
    VF                 original_cell_masses       () const;
    VS                 original_cells             () const; ///< intervals for each cell
    
    // approximated versions (for debug purpose)
    int                get_sorted_newton_system_ap( TridiagonalSymmetricMatrix<TF> &M, VF &V, TF &max_error_ratio, TF eps = 1e-20 );
    VF                 sorted_cell_barycenters_ap ( PI ni = 100000 ) const;
    VF                 sorted_cell_masses_ap      ( PI ni = 100000 ) const;
    void               check_newton_system        ( TF eps = 1e-6 );
    static void        plot_bnds_evolution        ( const std::vector<VB> &bnds );
    static bool        no_nan                     ( const VF &v );

    std::vector<VB> bnds;

    // values
    VI                 sorted_dirac_num_offsets;  ///< offsets in sorted_dirac_num_values. size = nb_sorted_diracs + 1
    VI                 sorted_dirac_num_values;   ///<
    VF                 sorted_dirac_positions;    ///<
    VF                 sorted_dirac_weights;      ///<
    VF                 sorted_dirac_masses;       ///<
    Density*           density;                   ///<
    
    mutable std::vector<VF> mass_history; ///< history of masses for each cell, indexed by time


    // parameters
    TF                 target_max_error_ratio     = 1e-4;
    bool               allow_ball_cut             = 1;
    int                verbosity                  = 0;
    TF                 epsilon                    = std::numeric_limits<TF>::max(); ///<

    private:
    void               initialize_dirac_positions ( const VF &dirac_positions, TF sep_equ_coords ); ///<
    void               initialize_dirac_weights   (); ///<
    void               initialize_dirac_masses    ( const VF &relative_dirac_masses, TF global_mass_ratio ); ///<

    PI                 nb_newton_iterations;
    TF                 global_mass_ratio;
};

    
} // namespace usdot

#include "SdotSystem.cxx" // IWYU pragma: keep
