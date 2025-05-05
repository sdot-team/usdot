#include <tl/support/P.h>
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

// void check_ders() {
//     GridDensity<TF,2> gd( { 1, 0.01, 1 } );

//     System<TF,GridDensity<TF,2>> si;
//     si.stream = &std::cout;
//     si.verbosity = 2;

//     // si.set_dirac_positions( Vec<TF>::cellspace( 1.6, 2.0, 15 ) );
//     si.set_dirac_positions( cellspace<TF>( 0.0, 2.0, 5 ) );
//     // si.set_dirac_positions( { -10, -11, -12, 4, 0, 1, 2 } );
//     si.set_global_mass_ratio( 1.0 );
//     si.set_density( &gd );

//     si.solve();

//     auto w0 = si.sorted_dirac_weights;

//     const TF eps = 1e-5;
//     gd.set_lag_ratio( 1 - eps );
//     si.newton_iterations();

//     auto w1 = si.sorted_dirac_weights;
//     for( PI i = 0; i < w0.size(); ++i )
//         P( ( w0[ i ] - w1[ i ] ) / eps );

//     si.plot();
//     P( si.l2_mass_error() );
// }

TEST_CASE( "System", "" ) {
    GridDensity<TF,2> gd( { 1, 0.01, 1 } );

    System<TF,GridDensity<TF,2>> si;
    si.stream = &std::cout;
    si.verbosity = 2;

    // si.set_dirac_positions( Vec<TF>::cellspace( 1.6, 2.0, 15 ) );
    si.set_dirac_positions( cellspace<TF>( 0.0, 2.0, 5 ) );
    // si.set_dirac_positions( { -10, -11, -12, 4, 0, 1, 2 } );
    si.set_global_mass_ratio( 1.0 );
    si.set_density( &gd );

    si.solve();
}
