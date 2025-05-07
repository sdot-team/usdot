#include <usdot/utility/glot.h>
#include <tl/support/P.h>

#include <usdot/DiffusionDensity.h>
#include "catch_main.h"

using namespace usdot;
using namespace std;
using TF = double;

TEST_CASE( "DiffusionDensity der", "" ) {
    std::vector<double> values( 50, 0 );
    for( PI i = 0; i < values.size(); ++i )
        values[ i ] = i > values.size() / 2 ? 1 : 0.1; //rand() / double( RAND_MAX );

    DiffusionDensity<TF> gd( 0, 1, values );
    gd.set_flattening_ratio( 0.015 );

    // glot_stream( [&]( std::ostream &fs ) {
    //     for( TF r = 0; r <= 1; r += 1. / 16 ) {
    //         gd.set_flattening_ratio( r );
    //         gd.plot( fs );
    //         // P( gd.integral( 0, 2 ) );
    //     }
    // } );

    TF eps = 1e-6;

    // DiffusionDensity<TF> ge = gd;
    // ge.set_flattening_ratio( 0.5 - eps );

    glot( linspace<TF>( -0.1, 1.1, 100 ), 
        [&]( TF x ) { return gd.value( x ) + 0.01; },
        [&]( TF x ) { return ( gd.primitive( x + eps ) - gd.primitive( x ) ) / eps; }
    );
}
