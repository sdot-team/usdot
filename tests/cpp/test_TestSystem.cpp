#include <tl/support/string/to_string.h>
#include <usdot/utility/glot.h>
#include <tl/support/P.h>

#include <usdot/TestSystem.h>
#include "catch_main.h"
#include <iostream>

#include <usdot/utility/gmp.h>
using namespace boost::multiprecision;
using TF = number<backends::cpp_bin_float<64>>;
// using TF = usdot::FP64;

using namespace usdot;
using namespace std;

#include "problemos.h"

TEST_CASE( "System", "" ) {
    // DiffusionDensity<TF> gd( { 1, 0.1, 0.1, 1 } );
    DiffusionDensity<TF> gd( dp, dv );
    // glot_stream( [&]( std::ostream &fs ) {
    //     for ( TF r : linspace<TF>( 0, 1, 10 ) ) {
    //         gd.set_flattening_ratio( r );
    //         gd.plot( fs, r == 0 ? "-" : "--", r == 0 || r == 1 ? 2 : 1 );
    //     }
    // } );

    TestSystem<TF> si;
    // si.set_dirac_positions( cellspace<TF>( 1.0, 2.0, 3 ) );
    si.set_global_mass_ratio( .984 );
    si.set_dirac_positions( pd, 1e-2 );
    // si.set_dirac_positions( { -1, 0, 0, 0, 0, 1 }, 1e-3 );
    si.set_density( &gd );
    si._update_system();

    // si.initialize_with_flat_density();

    // TF mi = 1;
    // TF ma = 0;
    // for( PI i = 1; i < si.sorted_dirac_masses.size(); ++i ) {
    //     mi = min( mi, si.sorted_dirac_positions[ i ] - si.sorted_dirac_positions[ i - 1 ] );
    //     ma = max( ma, si.sorted_dirac_positions[ i ] - si.sorted_dirac_positions[ i - 1 ] );
    // }
    // P( mi );
    // P( ma );
        
    si.stream = &std::cout;
    si.verbosity = 2;

    si.solve();
    si.plot();
}

// TEST_CASE( "System", "" ) {
//  // check_ders();

//     // glot_vec_ys( dv );
        
//     DiffusionDensity<TF> gd( dp, dv );

//     // glot_stream( [&]( std::ostream &fs ) {
//     // } );

//     System<TF> si;
//     si.set_global_mass_ratio( 0.984 );
//     si.set_dirac_positions( pd );
//     si.set_density( &gd );

//     si.stream = &std::cout;
//     si.verbosity = 2;

//     si.solve();
//     // si.plot();
// }
