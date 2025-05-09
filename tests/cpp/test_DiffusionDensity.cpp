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

// TEST_CASE( "DiffusionDensity ratio", "" ) {
//     glot_stream( [&]( std::ostream &fs ) {
//         for ( PI n : { 100, 1000 } ) {
//             const TF b = 0;
//             for ( TF e : { 10 } ) {
//                 // const TF e = 10;
//                 const TF m = ( b + e ) / 2;

//                 std::vector<TF> positions( n );
//                 for( PI i = 0; i < positions.size(); ++i )
//                     positions[ i ] = b + pow( TF( i ) / ( positions.size() - 1 ), 3 ) * ( e - b );
            
//                 std::vector<TF> values( n );
//                 for( PI i = 0; i < values.size(); ++i )
//                     values[ i ] = 0.5 + max( TF( -0.5 ), min( TF( 0.5 ), ( positions[ i ] - m ) / ( e - b ) * 10 ) );
            
//                 DiffusionDensity<TF> gd( positions, values );
//                 for( TF r = 0; r <= 1; r += 1. / 4 ) {
//                     gd.set_flattening_ratio( r );
//                     gd.plot( fs );
//                 }
//             }
//         }
//     } );
// }

TEST_CASE( "DiffusionDensity der", "" ) {
    std::vector<TF> positions( 50, 0 );
    for( PI i = 0; i < positions.size(); ++i )
        positions[ i ] = 2 * pow( TF( i ) / ( positions.size() - 1 ), 1 ); 

    std::vector<TF> values( positions.size(), 0 );
    for( PI i = 0; i < values.size(); ++i )
        values[ i ] = i > values.size() / 2 ? 1 : 0.1; //rand() / double( RAND_MAX );

    DiffusionDensity<TF> gd( positions, values );
    gd.set_flattening_ratio( 0.1 - 1e-6 );

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
            // for( TF r = 0; r <= 1; r += 1. / 16 ) {
            // // for( TF r : { 0.5, 0.9 } ) {
            //     gd.set_flattening_ratio( r );
            //     gd.plot( fs );
            //     // P( gd.integral( 0, 2 ) );
            // }
    //     }
    // } );
            
    TF eps = 1e-20;
    DiffusionDensity<TF> ge = gd;
    DiffusionDensity<TF> gf = gd;

    ge.set_flattening_ratio( gd.current_flattening_ratio - 1 * eps );
    gf.set_flattening_ratio( gd.current_flattening_ratio - 2 * eps );

    gd.compute_derivatives( 2 );
    ge.compute_derivatives( 2 );

    PI d = 2;
    for( TF x = 0.1; x < 1; x += 0.1 )
        P( gd.value( x, d ) - ( gd.value( x, d - 1 ) - ge.value( x, d - 1 ) ) / eps );

    glot( linspace<TF>( gd.min_x() - 0.1, gd.max_x() + 0.1, 1000 ), 
        [&]( TF x ) { return gd.value( x, d ) + 0.0; },
        [&]( TF x ) { return ( gd.value( x, d - 1 ) - ge.value( x, d - 1 ) ) / eps; }
    );
}
