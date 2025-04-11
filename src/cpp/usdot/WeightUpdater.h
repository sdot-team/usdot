#pragma once

#include "System.h"

namespace usdot {

/**
*/
template<class TF,class Density>
class WeightUpdater {
public:
    struct   Poly            { TF c0, c1, c2, ca; }; ///< poly to get optimal weight wrt radius of the first cell and blend ratio for target masses
    using    Sys             = System<TF,Density>;
    using    VC              = std::vector<std::array<PI,2>>;
    using    VP              = std::vector<Poly>;
    using    VF              = std::vector<TF>;
    using    PI              = std::size_t;
   
    /**/     WeightUpdater   ( Sys &sys );
   
    void     newton_update   ( PI max_iter = 1000 );
    void     run             ();
    
    // output
    PI       nb_iterations; ///<

private:
    TF       density_integral( TF x0, TF x1 ) const { return sys.density->integral( x0, x1 ); }
    int      get_weights_for ( VF &new_dirac_weights, const VC &connected_cells, const VP &polys, TF a );
    TF       density_value   ( TF a ) const { using namespace std; return max( 1e-3, sys.density->value( a ) );}
    TF       best_r_for_ib   ( const VP &polys, PI n, TF a ) const;
    void     get_system      ( VC &connected_cells, VP &polys, TF &max_a, TF &cell_error, int &has_bad_cell, const VF &sorted_dirac_weight ) const;
    VF       ce_errors       ( const VP &polys, PI nb, PI ne, TF a, TF r ) const;
    TF       bi_error        ( const VP &polys, PI n, TF a, TF r ) const;
    TF       ii_error        ( const VP &polys, PI n, TF a, TF r ) const;
    T_T void on_cell         ( const T &func ) const;
    
    T_T void on_newton_item  ( const T &func ) const;
    auto     newton_dir      () const -> std::pair<VF,TF>;
  
    Sys&     sys;            ///<
};

} // namespace usdot

#include "WeightUpdater.cxx"
