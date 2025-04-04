#pragma once

#include <tl/support/containers/Vec.h>
#include "GridDensity.h" // IWYU pragma: export

namespace usdot {
template<class TF,class Density>
class WeightInitializer;
    
/**
*/
template<class TF,class Density>
class System {
public:
    using          TI                               = std::size_t;
    using          VI                               = Vec<TI>;
    using          VF                               = Vec<TF>;

    /**/           System                           ();

    void           set_relative_dirac_masses        ( const VF &relative_dirac_masses );
    void           set_global_mass_ratio            ( TF mass_ratio );
    void           set_dirac_positions              ( const VF &dirac_positions );
    void           set_mass_ratio                   ( TF mass_ratio );
    void           set_density                      ( const Density *density );

    void           initialize_weights               ();
    void           update_weights                   ();

    VF             dirac_barycenters                () const;

    VF             cell_barycenters                 () const;
    VF             cell_masses                      () const;

    TF             density_value                    ( TF pos ) const;
    
    TF             l2_mass_error                    () const; ///< 
    TI             nb_diracs                        () const;
    void           plot                             ( std::string filename = "glot.py" ) const;

    // directly modifiable inputs
    TF             target_max_mass_error            = 1e-4;
    bool           multithread                      = false; ///<
    int            verbosity                        = 0; ///<
    
private:
    friend class   WeightInitializer<TF,Density>;

    void           _update_system                   ( bool need_weights = true ) const;

    void           for_each_cell                    ( auto &&func ) const;

    // 
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
