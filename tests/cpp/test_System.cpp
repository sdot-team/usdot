#include <partial1D/System.h>
#include <partial1D/glot.h>
#include "catch_main.h"

// #include <partial1D/gmp.h>
// using namespace boost::multiprecision;
// using TF = number<backends::cpp_bin_float<40>>;

using namespace usdot;
using namespace std;
using TF = double;

TEST_CASE( "System", "" ) {
    GridDensity<TF> gd( { 1, 0, 1 } );

    System<TF,GridDensity<TF>> si;
    si.set_dirac_positions( Vec<TF>::cellspace( 0.25, 1.75, 15 ) );
    si.set_global_mass_ratio( 0.25 );
    si.set_density( &gd );

    si.initialize_weights();
    si.plot();

    // P( si.cell_boundaries() );
    // P( si.cell_masses() );
    P( si.l2_mass_error() );

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
