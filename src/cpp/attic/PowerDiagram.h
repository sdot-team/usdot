#pragma once

#include "SymmetricBandMatrix.h"
#include "Density/Density.h"
#include "Interval.h"

namespace usdot {


// #include <boost/multiprecision/gmp.hpp>
// T_T void display( Displayer &ds, const boost::multiprecision::number<T> &n ) { ds << n.str(); }

/**
*/
template<class TF>
class PowerDiagram : public WithRefCount {
public:
    /**/               PowerDiagram             ( const Vec<TF> &seed_coords, const Vec<TF> &seed_weights, TF sep_equ_coords = 1e-6 );

    // 
    void               set_weights              ( const Vec<TF> &weights, bool sorted_nums = true ) const;
    Vec<TF>            get_weights              ( bool sorted_nums = true ) const;

    // visitors
    void               for_each_sub_interval    ( DensityIterator<TF> &density_iterator, const Interval<TF> &cell_interval, PI cell_num, auto &&func ) const;
    void               for_each_cell            ( auto &&func ) const; ///< func( Interval( ... ) )

    // basic information
    PI                 nb_cells                 () const { return sorted_seed_coords.size(); }
    void               display                  ( Displayer &ds ) const;
     
    // computations
    void               get_newton_system        ( SymmetricBandMatrix<TF> &M, Vec<TF> &V, PI &nb_arcs, const Density<TF> &density, TF coeff = 1 ) const;
    Vec<Vec<TF,2>>     cell_boundaries          ( bool sorted_nums = true ) const; ///< left and right boundary for each cell
    Vec<TF>            barycenters              ( const Density<TF> &density, bool sorted_nums = true ) const;
    Vec<TF>            masses                   ( const Density<TF> &density, bool sorted_nums = true ) const;
    
    // approximated versions (for debug purpose)
    void               get_newton_system_ap     ( SymmetricBandMatrix<TF> &M, Vec<TF> &V, PI &nb_arcs, const Density<TF> &density, TF eps = 1e-6 );
    Vec<TF>            barycenters_ap           ( const Density<TF> &density, bool sorted_nums = true, PI ni = 10000 ) const;
    Vec<TF>            masses_ap                ( const Density<TF> &density, bool sorted_nums = true, PI ni = 10000 ) const;
       
    // attributes
    Vec<TF>            sorted_seed_weights;     ///<
    Vec<TF>            sorted_seed_coords;      ///<
    Vec<PI>            sorted_seed_nums;        ///<
    bool               allow_ball_cut           = 1;
    TF                 epsilon                  = 0; ///<
    
private: 
    void               __for_each_sub_interval  ( DensityIterator<TF> &density_iterator, const Interval<TF> &cell_interval, PI cell_num, auto &&func ) const; ///< after handling of epsilon
    void               _for_each_sub_interval   ( DensityIterator<TF> &density_iterator, const Interval<TF> &cell_interval, PI cell_num, auto &&func ) const; ///< after x0, x1 swapping (if necessary)
};

    
} // namespace usdot

#include "PowerDiagram.cxx" // IWYU pragma: keep
