#pragma once

#include "utility/TridiagonalSymmetricMatrix.h"
#include "UsdotDensity.h"

namespace usdot {

/**
*/
template<class TF,class Density=UsdotDensity<TF>>
class UsdotSystem {
public:
    using              VB                         = std::vector<std::array<TF,2>>;
    using              VF                         = std::vector<TF>;
    using              VI                         = std::vector<PI>;
  
    /**/               UsdotSystem                ( Density *density, const VF &dirac_positions, TF global_mass_ratio = 1, TF sep_equ_coords = 1e-6 );
  
    //   
    void               set_original_weights       ( const VF &weights ) const; ///<
    VF                 get_original_weights       () const; ///<
  
    int                newton_iterations          ( const TF min_relax = 1e-6 );
    void               solve                      ();
    
    // basic information  
    PI                 nb_original_diracs         () const { return sorted_dirac_num_values.size(); }
    PI                 nb_sorted_diracs           () const { return sorted_dirac_positions.size(); }
    void               plot                       ( Str filename = "glot.py" ) const;
    
    #ifdef TL_DISPLAYER_IS_DEFINED
    void               display                    ( Displayer &ds ) const;
    #endif
       
    // computations  
    int                get_sorted_newton_system   ( TridiagonalSymmetricMatrix<TF> &M, VF &V, TF &max_error_ratio ) const;
    VF                 sorted_cell_barycenters    () const;
    VB                 sorted_cell_boundaries     () const;
    VF                 sorted_cell_masses         () const;
    int                newton_dir                 ( VF &dir, TF &max_error_ratio ) const;
      
    void               for_each_cell              ( auto &&func ) const;

    VF                 original_cell_barycenters  () const;
    VF                 original_cell_masses       () const;
    
    // approximated versions (for debug purpose)
    int                get_sorted_newton_system_ap( TridiagonalSymmetricMatrix<TF> &M, VF &V, TF &max_error_ratio, TF eps = 1e-20 );
    VF                 sorted_cell_barycenters_ap ( PI ni = 100000 ) const;
    VF                 sorted_cell_masses_ap      ( PI ni = 100000 ) const;
    void               check_newton_system        ( TF eps = 1e-6 );
    static void        plot_bnds_evolution        ( const std::vector<VB> &bnds );
    static bool        no_nan                     ( const VF &v );

    // values
    VI                 sorted_dirac_num_offsets;  ///< offsets in sorted_dirac_num_values. size = nb_sorted_diracs + 1
    VI                 sorted_dirac_num_values;   ///<
    VF                 sorted_dirac_positions;    ///<
    VF                 sorted_dirac_weights;      ///<
    VF                 sorted_dirac_masses;       ///<
    PI                 beg_unknowns;              ///<
    PI                 end_unknowns;              ///<
    Density*           density;                   ///<
    
    // parameters
    TF                 target_max_error_ratio     = 1e-4;
    int                verbosity                  = 0;

    // outputs
    PI                 nb_newton_iterations       = 0;

private:
    void               initialize_dirac_positions ( const VF &dirac_positions, TF sep_equ_coords ); ///<
    void               initialize_dirac_weights   (); ///<
    void               initialize_dirac_masses    ( const VF &relative_dirac_masses, TF global_mass_ratio ); ///<

    TF                 global_mass_ratio;
};

    
} // namespace usdot

#include "UsdotSystem.cxx" // IWYU pragma: keep
