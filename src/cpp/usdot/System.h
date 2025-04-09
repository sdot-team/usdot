#pragma once

#include "utility/common_macros.h"
#include "GridDensity.h" // IWYU pragma: export

namespace usdot {
template<class TF,class Density> class WeightInitializer;
template<class TF,class Density> class WeightUpdater;
    
/**
*/
template<class TF,class Density>
class System {
public:
    using          VI                               = std::vector<PI>;
    using          VF                               = std::vector<TF>;
    using          VB                               = std::vector<std::array<TF,2>>;

    /**/           System                           ();

    void           set_relative_dirac_masses        ( const VF &relative_dirac_masses );
    void           set_global_mass_ratio            ( TF mass_ratio );
    void           set_dirac_positions              ( const VF &dirac_positions );
    void           set_mass_ratio                   ( TF mass_ratio );
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
    
    TF             l2_mass_error                    () const; ///< 
    PI             nb_diracs                        () const;
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
    VF             sorted_dirac_positions;          ///<
    mutable VF     sorted_dirac_weights;            ///<
    mutable VF     sorted_dirac_masses;             ///<
    VI             sorted_dirac_nums;               ///<
    mutable TF     total_dirac_mass;                ///<
   
    // density
    const Density* density;                         ///<
};


} // namespace usdot

#include "System.cxx" // IWYU pragma: export
