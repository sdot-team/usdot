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
    /**/               PowerDiagram             ( const Vec<TF> &seed_coords, const Vec<TF> &seed_weights );

    // 
    void               set_weights              ( const Vec<TF> &weights ) const;
    Vec<TF>            get_weights              () const;

    // visitors
    void               for_each_sub_interval    ( DensityIterator<TF> &density_iterator, const Interval<TF> &cell_interval, auto &&func ) const;
    void               for_each_cell            ( auto &&func ) const; ///< func( Interval( ... ) )

    // basic information
    PI                 nb_cells                 () const { return sorted_seed_coords.size(); }
    void               display                  ( Displayer &ds ) const;
     
    // computations
    void               get_newton_system        ( SymmetricBandMatrix<TF> &M, Vec<TF> &V, PI &nb_arcs, const Density<TF> &density, TF coeff = 1 ) const;
    Vec<TF>            cell_boundaries          () const; ///< left and right boundary for each cell
    Vec<TF>            barycenters              ( const Density<TF> &density ) const;
    TF                 integral                 ( const Density<TF> &density, TF x0, TF x1 ) const;
    Vec<TF>            masses                   ( const Density<TF> &density ) const;
    
    // approximated versions (for debug purpose)
    void               get_newton_system_ap     ( SymmetricBandMatrix<TF> &M, Vec<TF> &V, PI &nb_arcs, const Density<TF> &density, TF eps = 1e-6 ) const;
    Vec<TF>            barycenters_ap           ( const Density<TF> &density, PI ni = 10000 ) const;
    TF                 integral_ap              ( const Density<TF> &density, TF x0, TF x1, PI ni = 10000 ) const;
    Vec<TF>            masses_ap                ( const Density<TF> &density, PI ni = 10000 ) const;
       
    // attributes
    Vec<TF>            sorted_seed_weights;     ///<
    Vec<TF>            sorted_seed_coords;      ///<
    Vec<PI>            sorted_seed_nums;        ///<
    
private: 
    void               _for_each_sub_interval   ( DensityIterator<TF> &density_iterator, const Interval<TF> &cell_interval, auto &&func ) const;
};

    
} // namespace usdot

#include "PowerDiagram.cxx" // IWYU pragma: keep
