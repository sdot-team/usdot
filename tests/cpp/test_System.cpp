#include <usdot/utility/glot.h>
#include <usdot/System.h>
#include "catch_main.h"
#include <iostream>

// #include <usdot/utility/gmp.h>
// using namespace boost::multiprecision;
// using TF = number<backends::cpp_bin_float<40>>;

using namespace usdot;
using namespace std;
using TF = double;

TEST_CASE( "System", "" ) {
    GridDensity<TF> gd( { 1, 0, 1 } );

    System<TF,GridDensity<TF>> si;
    si.stream = &std::cout;
    si.verbosity = 2;

    // si.set_dirac_positions( Vec<TF>::cellspace( 1.6, 2.0, 15 ) );
    si.set_dirac_positions( cellspace<TF>( 0.0, 2.0, 5 ) );
    si.set_global_mass_ratio( 1 );
    si.set_density( &gd );

    si.initialize_weights();
    //si.update_weights();
    //si.solve();
    si.plot();
    P( si.l2_mass_error() );

    // P( si.cell_boundaries() );
    // P( si.cell_masses() );
    // P( si.l2_mass_error() );

    // si.dirac_positions = Vec<TF>::cellspace( 0, 1, 200 );
    // si.starting_filter_value = 0.0; //99;
    // si.target_filter_value = 0.0;
    // si.global_mass_ratio = 2.0 / 3.0;

    // si.density_values = { 1, 1, 1e-3, 1e-3, 1e-3, 1, 1 };
    // si.beg_x_density = -1;
    // si.end_x_density = 2;
 
    // ConvGridSolver<TF> solver( std::move( si ) );
    // // P( solver.normalized_cell_boundaries() );
    // solver.solve(); // 85 -> 74
    // solver.plot();

    // glot( Vec<TF>::linspace( -1, 2, 1000 ), 
    //     [&]( TF x ) { return solver.density_value( x ); }
    // );
}
