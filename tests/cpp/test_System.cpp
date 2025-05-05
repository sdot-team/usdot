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

// void check_ders( const TF eps = 1e-5 ) {
//     GridDensity<TF,2> gd( { 1, 0.01, 1 } );

//     System<TF,GridDensity<TF,2>> si;
//     si.stream = &std::cout;
//     si.verbosity = 2;

//     si.set_dirac_positions( cellspace<TF>( 0.0, 2.0, 5 ) );
//     si.set_global_mass_ratio( 1.0 );
//     si.set_density( &gd );
    
//     // w0
//     si.initialize_with_flat_density();
//     si.newton_iterations();
//     auto w0 = si.sorted_dirac_weights;
//     P( si.der_weights_wrt_lap_ratio( 2 ) );

//     // w1
//     gd.set_lag_ratio( 1 - eps );
//     si.newton_iterations();
//     auto w1 = si.sorted_dirac_weights;

//     // w2
//     gd.set_lag_ratio( 1 - 2 * eps );
//     si.newton_iterations();
//     auto w2 = si.sorted_dirac_weights;

//     //  
//     for( PI i = 0; i < w0.size(); ++i )
//         P( ( w0[ i ] - w1[ i ] ) / eps );
//     for( PI i = 0; i < w0.size(); ++i )
//         P( ( w0[ i ] - 2 * w1[ i ] + w2[ i ] ) / eps );
// }

TEST_CASE( "System", "" ) {
    // check_ders();
    std::vector<TF> dv( 100 );
    for( PI i = 0; i < dv.size(); ++i )
        dv[ i ] = ( i > dv.size() / 2 ? 0.1 : 1 );
        
    GridDensity<TF,2> gd( dv );
    // GridDensity<TF,2> ge = gd;
    // GridDensity<TF,2> gf = gd;
    // gd.set_lag_ratio( 1. );
    // ge.set_lag_ratio( .999 );
    // gf.set_lag_ratio( .998 );

    // glot( linspace<TF>( -1.0, dv.size(), 100 ), 
    //     [&]( TF x ) { return gd.value( x ); },
    //     [&]( TF x ) { return ge.value( x ); },
    //     [&]( TF x ) { return gf.value( x ); }
    // );

    System<TF,GridDensity<TF,2>> si;
    si.stream = &std::cout;
    si.verbosity = 2;

    // si.set_dirac_positions( Vec<TF>::cellspace( 1.6, 2.0, 15 ) );
    si.set_dirac_positions( cellspace<TF>( 0.0, dv.size() / 2.0, 10 ) );
    // si.set_dirac_positions( { -10, -11, -12, 4, 0, 1, 2 } );
    si.set_global_mass_ratio( 0.8 );
    si.set_density( &gd );

    si.solve();
}
