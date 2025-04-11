#include <tl/support/P.h>

#include <usdot/GridDensity.h>
#include <usdot/utility/glot.h>
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

TEST_CASE( "GridDensity cdf", "" ) {
    GridDensity<TF,0> gd( { 0.5, 1.7, 2.3, 3, 4 }, { 1, 2, 1, 2, 2 } );
    // GridDensity<TF,1> gd( 1, 3, { 1, 2, 1, 2 } );
    // GridDensity<TF,2> gd( { 1, 2, 1, 2 } );

    TF ed = 1e-6;
    glot( linspace<TF>( 0.01, 0.99, 100 ), 
        [&]( TF x ) { return ( gd.inv_cdf( x + 1e-6, ed ) - gd.inv_cdf( x, ed ) ) / 1e-6; },
        [&]( TF x ) { return gd.der_inv_cdf( x, ed ) + 0.01; }
    );
}
