#pragma once

#include "utility/TridiagonalSymmetricMatrix.h"
#include "Interval.h"

namespace usdot {

/**
*/
template<class TF>
class PowerDiagram {
public:
    using              VC                       = std::vector<Interval<TF>>;
    using              VS                       = std::vector<VC>;
    using              VF                       = std::vector<TF>;
    using              VI                       = std::vector<PI>;

    /**/               PowerDiagram             ( const VF &seed_coords, TF sep_equ_coords = 1e-6 );

    // 
    void               set_original_weights     ( const VF &weights ) const; ///<
    VF                 get_original_weights     () const; ///<

    // visitors
    T_T void           for_each_interval        ( T &&func ) const;

    // basic information
    PI                 nb_original_seeds        () const { return sorted_seed_num_values.size(); }
    PI                 nb_sorted_seeds          () const { return sorted_seed_coords.size(); }
    T_T void           display                  ( T &ds ) const;
     
    // computations
    T_T void           get_sorted_newton_system ( TridiagonalSymmetricMatrix<TF> &M, VF &V, PI &nb_arcs, const T &density ) const;
    T_T VF             sorted_barycenters       ( const T &density ) const;
    T_T VF             sorted_masses            ( const T &density ) const;
    VS                 sorted_cells             () const; ///< left and right boundary for each cell
    
    T_T VF             original_barycenters     ( const T &density ) const;
    T_T VF             original_masses          ( const T &density ) const;
    VS                 original_cells           () const; ///< left and right boundary for each cell
    
    // approximated versions (for debug purpose)
    T_T void           get_newton_system_ap     ( TridiagonalSymmetricMatrix<TF> &M, VF &V, PI &nb_arcs, const T &density, TF eps = 1e-6 );
    T_T VF             barycenters_ap           ( const T &density, bool sorted_nums = true, PI ni = 10000 ) const;
    T_T VF             masses_ap                ( const T &density, bool sorted_nums = true, PI ni = 10000 ) const;
       
    // attributes
    VI                 sorted_seed_num_offsets; ///< offsets in sorted_dirac_num_values. size = nb_sorted_diracs + 1
    VI                 sorted_seed_num_values;  ///<
    VF                 sorted_seed_weights;     ///<
    VF                 sorted_seed_masses;      ///<
    VF                 sorted_seed_coords;      ///<
    bool               allow_ball_cut           = 1;
    TF                 epsilon                  = 0; ///<
    
private: 
    void               __for_each_sub_interval  ( DensityIterator<TF> &density_iterator, const Interval<TF> &cell_interval, PI cell_num, auto &&func ) const; ///< after handling of epsilon
    void               _for_each_sub_interval   ( DensityIterator<TF> &density_iterator, const Interval<TF> &cell_interval, PI cell_num, auto &&func ) const; ///< after x0, x1 swapping (if necessary)
};

    
} // namespace usdot

#include "PowerDiagram.cxx" // IWYU pragma: keep
