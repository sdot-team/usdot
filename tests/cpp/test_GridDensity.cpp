#include <tl/support/P.h>

#include <usdot/GridDensity.h>
#include <usdot/utility/glot.h>
#include <vector>
#include "catch_main.h"
#include "usdot/utility/linspace.h"

using namespace usdot;
using namespace std;
using TF = double;

// TEST_CASE( "GridDensity", "" ) {
//     GridDensity<TF,0> gd( { 0.5, 1.7, 2.3, 3 }, { 1, 2, 1, 2 } );
//     // GridDensity<TF,1> gd( 1, 3, { 1, 2, 1, 2 } );
//     // GridDensity<TF,2> gd( { 1, 2, 1, 2 } );

//     // glot( linspace<TF>( -1, 4, 1000 ), 
//     //     [&]( TF x ) { return ( gd.x_primitive( x + 1e-6 ) - gd.x_primitive( x ) ) / 1e-6; },
//     //     [&]( TF x ) { return x * gd.value( x ) + .01; }
//     // );

//     // P( gd.cdf( 10 ) );
//     // P( gd.cdf( 0 ) );
// }

// TEST_CASE( "GridDensity cdf", "" ) {
//     GridDensity<TF,0> gd( { 0.5, 1.7, 2.3, 3, 4 }, { 1, 2, 1, 2, 2 } );
//     // GridDensity<TF,1> gd( 1, 3, { 1, 2, 1, 2 } );
//     // GridDensity<TF,2> gd( { 1, 2, 1, 2 } );
//     // gd.set_lag_ratio( 0.5 );

//     GridDensity<TF,0> g0 = gd;
//     GridDensity<TF,0> g1 = gd;
//     GridDensity<TF,0> ga = gd;
//     GridDensity<TF,0> gb = gd;
//     GridDensity<TF,0> gc = gd;
//     g0.set_lag_ratio( 0.00 );
//     g1.set_lag_ratio( 1.00 );
//     ga.set_lag_ratio( 0.25 );
//     gb.set_lag_ratio( 0.50 );
//     gc.set_lag_ratio( 0.75 );
//     P( g0.integral( -1, 20 ) );
//     P( g1.integral( -1, 20 ) );
//     P( ga.integral( -1, 20 ) );
//     P( gb.integral( -1, 20 ) );
//     P( gc.integral( -1, 20 ) );
//     // P( ga.values );

//     TF ed = 1e-6;
//     glot( linspace<TF>( 0.0, 5.0, 100 ), 
//         // [&]( TF x ) { return ( gd.inv_cdf( x + 1e-6, ed ) - gd.inv_cdf( x, ed ) ) / 1e-6; },
//         [&]( TF x ) { return g0.value( x ); },
//         [&]( TF x ) { return g1.value( x ); },
//         [&]( TF x ) { return ga.value( x ); },
//         [&]( TF x ) { return gb.value( x ); },
//         [&]( TF x ) { return gc.value( x ); }
//     );
// }
TEST_CASE( "GridDensity der", "" ) {
    std::vector<double> values( 10, 0 );
    for( PI i = 0; i < values.size(); ++i )
        values[ i ] = rand() / double( RAND_MAX );
 
    TF eps = 1e-6;
    TF lr = 1;
    GridDensity<TF,2> gd( values );
    gd.set_lag_ratio( lr );

    GridDensity<TF,2> ge = gd;
    ge.set_lag_ratio( lr - eps );

    glot( linspace<TF>( -1.0, 11.0, 100 ), 
        [&]( TF x ) { return gd.value( x, 3 ); },
        [&]( TF x ) { return ( gd.value( x, 2 ) - ge.value( x, 2 ) ) / eps; }
    );
}
