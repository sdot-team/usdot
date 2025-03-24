#include <tl/support/containers/Vec.h>
#include <partial1D/GridDensity.h>
#include <partial1D/glot.h>
#include "catch_main.h"

using namespace usdot;
using namespace std;
using TF = double;

// TEST_CASE( "Grid density", "" ) {
//     PI mul_x = 10;
//     GridDensity<TF> gd( { 1, 0, 0, 1 }, 0.8, mul_x );
//     GridDensity<TF> ge( { 1, 0, 0, 1 }, 0.8 + 1e-6, mul_x );

//     glot( Vec<TF>::linspace( -3, 7, 1000 ), 
//         [&]( TF x ) { return gd.der_value( x * mul_x ); },
//         [&]( TF x ) { return ( ge.value( x * mul_x ) - gd.value( x * mul_x ) ) / 1e-6 + .1; }
//     );
// }

// TEST_CASE( "Grid density no filter", "" ) {
//     GridDensity<TF> gd( { 1, 0, 0, 1 }, 0.0 );

//     glot( Vec<TF>::linspace( -2, 6, 1000 ), 
//         [&]( TF x ) { return gd.value( x ); }
//     );
// }

TEST_CASE( "Grid primitive", "" ) {
    PI mul_x = 10;
    GridDensity<TF> gd( { 1, 0, 0, 1, 1 }, 0.9, mul_x );

    // glot( Vec<TF>::linspace( -3, 7, 1000 ), 
    //     [&]( TF x ) { return gd.value( x * mul_x ); },
    //     [&]( TF x ) { return ( gd.primitive( x * mul_x + 1e-6 ) - gd.primitive( x * mul_x ) ) / 1e-6 + .01; }
    // );
    // glot( Vec<TF>::linspace( -3, 7, 1000 ), 
    //     [&]( TF x ) { return x * mul_x * gd.value( x * mul_x ); },
    //     [&]( TF x ) { return ( gd.x_primitive( x * mul_x + 1e-6 ) - gd.x_primitive( x * mul_x ) ) / 1e-6 + .01; }
    // );
    glot( Vec<TF>::linspace( -7, 12, 1000 ), 
        [&]( TF x ) { return gd.der_value( x * mul_x ); },
        [&]( TF x ) { return ( gd.der_primitive( x * mul_x + 1e-6 ) - gd.der_primitive( x * mul_x ) ) / 1e-6 + .1; }
    );
}
