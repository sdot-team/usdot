#include <usdot/utility/glot.h>
#include <tl/support/P.h>

#include <usdot/DiffusionDensity.h>
#include "catch_main.h"

#include <usdot/utility/gmp.h>
using namespace boost::multiprecision;
using TF = number<backends::cpp_bin_float<64>>;

using namespace usdot;
using namespace std;
// using TF = double;

TEST_CASE( "DiffusionDensity der", "" ) {
    std::vector<TF> positions( 50, 0 );
    for( PI i = 0; i < positions.size(); ++i )
        positions[ i ] = pow( TF( i ) / ( positions.size() - 1 ), 1 ); 

    std::vector<TF> values( positions.size(), 0 );
    for( PI i = 0; i < values.size(); ++i )
        values[ i ] = i > values.size() / 2 ? 1 : 0.1; //rand() / double( RAND_MAX );

    DiffusionDensity<TF> gd( positions, values );
    gd.set_flattening_ratio( 0.5 - 1e-6 );

    // glot_stream( [&]( std::ostream &fs ) {
    //     for ( PI n : { 100 } ) {
    //         std::vector<TF> positions( n, 0 );
    //         for( PI i = 0; i < positions.size(); ++i )
    //             positions[ i ] = pow( TF( i ) / ( positions.size() - 1 ), 2 ); 
        
    //         std::vector<TF> values( n, 0 );
    //         for( PI i = 0; i < values.size(); ++i )
    //             values[ i ] = positions[ i ] > 0.5 ? 1 : 0.1; //rand() / double( RAND_MAX );
        
    //         DiffusionDensity<TF> gd( positions, values );
    //         gd.set_flattening_ratio( 0.5 );
    //         for( TF r = 0; r <= 1; r += 1. / 16 ) {
    //         // for( TF r : { 0.5, 0.9 } ) {
    //             gd.set_flattening_ratio( r );
    //             gd.plot( fs );
    //             // P( gd.integral( 0, 2 ) );
    //         }
    //     }
    // } );
            
    TF eps = 1e-16;
    DiffusionDensity<TF> ge = gd;
    DiffusionDensity<TF> gf = gd;

    ge.set_flattening_ratio( gd.current_flattening_ratio - 1 * eps );
    gf.set_flattening_ratio( gd.current_flattening_ratio - 2 * eps );

    gd.compute_derivatives( 2 );
    ge.compute_derivatives( 2 );

    for( TF x = 0.1; x < 1; x += 0.1 )
        P( gd.primitive( x, 1 ), ( gd.primitive( x, 0 ) - ge.primitive( x, 0 ) ) / eps );

    glot( linspace<TF>( -0.1, 1.1, 100 ), 
        [&]( TF x ) { return gd.primitive( x, 2 ) + 0.0; },
        [&]( TF x ) { return ( gd.primitive( x, 1 ) - ge.primitive( x, 1 ) ) / eps; }
    );
}
