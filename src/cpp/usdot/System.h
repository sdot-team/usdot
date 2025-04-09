#pragma once

#include "utility/common_macros.h"
#include "GridDensity.h"

namespace usdot {
template<class TF,class Density> class WeightInitializer;
template<class TF,class Density> class WeightUpdater;
    
/**
*/
template<class TF,class Density=GridDensity<TF>>
class System {
public:
    using          VI                               = std::vector<PI>;
    using          VF                               = std::vector<TF>;
    using          VB                               = std::vector<std::array<TF,2>>;

    /**/           System                           ();

    void           set_relative_dirac_masses        ( const VF &relative_mass_ratios );
    void           set_global_mass_ratio            ( const TF &global_mass_ratio );
    void           set_dirac_positions              ( const VF &dirac_positions, TF min_dirac_separation = 1e-6 );
    void           set_mass_ratio                   ( const TF &mass_ratio );
    void           set_density                      ( const Density *density );

    void           initialize_weights               ();
    void           update_weights                   ();
    void           solve                            ();

    VF             dirac_positions                  () const;
    VF             dirac_weights                    () const;

    VF             cell_barycenters                 () const;
    VB             cell_boundaries                  () const;
    VF             cell_masses                      () const;

    TF             density_value                    ( TF pos ) const;
    
    PI             nb_original_diracs               () const;
    PI             nb_sorted_diracs                 () const;
    TF             l2_mass_error                    () const; ///< 
    TF             x_tol                            () const; ///<
    void           plot                             ( std::string filename = "glot.py" ) const;

    // directly modifiable inputs
    TF             target_max_mass_error            = 1e-4;
    bool           multithread                      = false; ///<
    int            verbosity                        = 0; ///<
    std::ostream*  stream                           = nullptr; ///<
    
    // output
    PI             nb_iterations_update             = 0; ///<
    PI             nb_iterations_init               = 0; ///<

private:
    friend class   WeightInitializer                <TF,Density>;
    friend class   WeightUpdater                    <TF,Density>;

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
    mutable TF     total_dirac_mass;                ///<
   
    // density
    const Density* density;                         ///<
};


} // namespace usdot

#include "System.cxx" // IWYU pragma: export
