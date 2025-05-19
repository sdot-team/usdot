#pragma once

#include "utility/TridiagonalSymmetricMatrix.h"
#include "DiffusionDensity.h"

namespace usdot {
    
/**
*/
template<class TF,class Density=DiffusionDensity<TF>>
class System {
public:
    using          ConnectedCells                   = std::vector<std::array<PI,2>>;
    struct         Poly                             { TF c0, c1, c2, ca; }; ///< poly to get optimal weight wrt radius of the first cell and blend ratio for target masses
    using          TSM                              = TridiagonalSymmetricMatrix<TF>;
    using          VI                               = std::vector<PI>;
    using          VF                               = std::vector<TF>;
    using          MF                               = std::vector<VF>;
    using          VB                               = std::vector<std::array<TF,2>>;

    /**/           System                           ();

    void           set_relative_dirac_masses        ( const VF &relative_mass_ratios );
    void           set_global_mass_ratio            ( const TF &global_mass_ratio );
    void           set_dirac_positions              ( const VF &dirac_positions, TF min_dirac_separation = 1e-6 );
    void           set_mass_ratio                   ( const TF &mass_ratio );
    void           set_density                      ( Density *density );

    void           initialize_with_flat_density     ();
    MF             der_weights_wrt_flat_ratio       ( PI nb_ders = 2, bool use_approx_for_ders = false ); ///< assuming we're on a solution
    int            newton_iterations                ( TF min_relax = 1e-6 ); ///< return 0 if OK
    void           solve_using_cdf                  ();
    void           solve                            ( bool use_approx_for_ders = false );

    VF             dirac_positions                  () const;
    VF             dirac_weights                    () const;

    VF             cell_barycenters                 () const;
    VB             cell_boundaries                  () const;
    VF             cell_masses                      () const;
    TF             cost                             () const;

    TF             density_value                    ( TF pos ) const;
    
    PI             nb_original_diracs               () const;
    PI             nb_sorted_diracs                 () const;
    TF             max_mass_error                   () const; ///< 
    TF             l2_mass_error                    ( bool max_if_bad_cell = false ) const; ///< 
    TF             x_tol                            () const; ///<
    void           plot                             ( std::string filename = "glot.py" ) const;

    int            get_weights_for                  ( VF &new_dirac_weights, const ConnectedCells &connected_cells, const std::vector<Poly> &polys, TF a );
    void           get_system_con                   ( ConnectedCells &connected_cells, std::vector<Poly> &opt_weights, TF &max_a, TF &cell_error, int &has_bad_cell, const VF &sorted_dirac_weights ) const;
    TF             best_r_for_ib                    ( const std::vector<Poly> &polys, PI n, TF a ) const;
    int            newton_con                       ( TF min_a );
    T_T static T   sgn                              ( const T &a ) { return ( a > 0 ) - ( a < 0 ); }

    // directly modifiable inputs
    TF             target_max_mass_error            = 1e-4;
    TF             target_newton_error              = std::numeric_limits<TF>::epsilon();
    bool           multithread                      = false; ///<
    int            verbosity                        = 0; ///<
    std::ostream*  stream                           = nullptr; ///<
    
    // output
    PI             nb_newton_iterations             = 0; ///<
    PI             nb_iterations_init               = 0; ///<
    TF             time_in_update                   = 0; ///<
    TF             time_in_init                     = 0; ///<

// private:
    T_T void       _for_each_unintersected_cell     ( const T &func ) const;
    T_T void       _for_each_newton_item            ( PI num_der, const T &func );
    T_T void       _for_each_newton_item            ( const T &func ) const;
    int            _make_newton_system              ();
    void           _update_system                   ( bool need_weights = true ) const;

    static void    plot_bnds_evolution              ( const std::vector<VB> &bnds );
    T_T void       for_each_cell                    ( const T &func ) const;

    // 
    VF             relative_mass_ratios;            ///<
    TF             global_mass_ratio;               ///<

    // diracs 
    VI             sorted_dirac_num_offsets;        ///< offsets in sorted_dirac_num_values. size = nb_sorted_diracs + 1
    VI             sorted_dirac_num_values;         ///<
    VF             sorted_dirac_positions;          ///<
    mutable VF     sorted_dirac_weights;            ///<
    mutable VF     sorted_dirac_masses;             ///<
   
    // density
    Density*       density;                         ///<

    //
    TSM            newton_matrix_ldlt;              ///<
    VF             newton_vector;                   ///<
    TF             newton_error;                    ///<
};


} // namespace usdot

#include "System.cxx" // IWYU pragma: export
