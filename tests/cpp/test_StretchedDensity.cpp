#include <tl/support/P.h>

// #include "usdot/utility/linspace.h"
#include <usdot/StretchedDensity.h>
#include <usdot/utility/glot.h>
#include "catch_main.h"

using namespace usdot;
using namespace std;
using TF = double;

TEST_CASE( "StretchedDensity", "" ) {
    const TF eps = 1e-4;

    std::vector<double> values( 100, 0 );
    for( PI i = 0; i < values.size(); ++i )
        values[ i ] = ( i > values.size() * 1 / 3 && i < values.size() * 2 / 3 ? 1 : 0.1 );
    P( values );

    StretchedDensity<TF> gd( values );
    // gd.set_flattening_ratio( 0.5 );

    // glot( linspace<TF>( -0.1, 1.1, 1000 ), 
    //     [&]( TF x ) { return gd.value( x ) + 0.01; },
    //     [&]( TF x ) { return ( gd.primitive( x + eps ) - gd.primitive( x ) ) / eps; }
    // );

    glot_stream( [&]( std::ostream &fs ) {
        for( TF r = 0; r <= 1; r += 1. / 16 ) {
            gd.set_flattening_ratio( r );
            gd.plot( fs );
            P( gd.integral( 0, 2 ) );
        }
    } );
    
    // StretchedDensity<TF> ge( values );
    // ge.set_flattening_ratio( 0.5 - eps );

    // glot( linspace<TF>( -0.1, 1.1, 1000 ), 
    //     // [&]( TF x ) { return gd.value( x ); }
    //     [&]( TF x ) { return ( gd.value( x ) - ge.value( x ) ) / eps; }
    // );
}
