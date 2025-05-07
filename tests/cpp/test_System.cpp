#include <usdot/utility/glot.h>
#include <tl/support/P.h>

#include <usdot/System.h>
#include "catch_main.h"
#include <iostream>

#include <usdot/utility/gmp.h>
using namespace boost::multiprecision;
using TF = number<backends::cpp_bin_float<40>>;

using namespace usdot;
using namespace std;
// using TF = FP80;

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
        dv[ i ] = ( i > dv.size() / 2 ? 0.5 : 1 );
        
    DiffusionDensity<TF> gd( dv );

    System<TF> si;
    si.set_dirac_positions( cellspace<TF>( 0.0, dv.size() / 2.0, 1000 ) );
    si.set_global_mass_ratio( 0.9 );
    si.set_density( &gd );

    si.stream = &std::cout;
    si.verbosity = 2;

    // si.solve();

    gd.set_flattening_ratio( 1 );
    si.initialize_with_flat_density();
    // si.newton_iterations();
    // auto w0 = si.sorted_dirac_weights;
    
    // const TF eps = 1e-2;
    // gd.set_flattening_ratio( 1 - 1 * eps );
    // si.newton_iterations();
    // auto w1 = si.sorted_dirac_weights;

    // gd.set_flattening_ratio( 1 - 2 * eps );
    // si.newton_iterations();
    // auto w2 = si.sorted_dirac_weights;

    // gd.set_flattening_ratio( 1 - 3 * eps );
    // si.newton_iterations();
    // auto w3 = si.sorted_dirac_weights;

    // std::vector<TF> arrs;
    // std::vector<TF> brrs;
    // std::vector<TF> crrs;
    // std::vector<TF> drrs;
    // for( TF a = 0.0; a < 0.5; a += 1e-3 ) {
    //     gd.set_flattening_ratio( 1 - a );

    //     si.sorted_dirac_weights = w0;
    //     if ( si.l2_mass_error() < 1e10 )
    //         arrs.push_back( si.l2_mass_error() );

    //     for( PI i = 0; i < w0.size(); ++i )
    //         si.sorted_dirac_weights[ i ] = w0[ i ] + ( w1[ i ] - w0[ i ] ) / eps * a;
    //     if ( si.l2_mass_error() < 1e10 )
    //         brrs.push_back( si.l2_mass_error() );

    //     for( PI i = 0; i < w0.size(); ++i )
    //         si.sorted_dirac_weights[ i ] = w0[ i ] + ( w1[ i ] - w0[ i ] ) / eps * a + ( w0[ i ] + w2[ i ] - 2 * w1[ i ] ) / pow( eps, 2 ) * pow( a, 2 ) / 2;
    //     if ( si.l2_mass_error() < 1e10 )
    //         crrs.push_back( si.l2_mass_error() );

    //     for( PI i = 0; i < w0.size(); ++i )
    //         si.sorted_dirac_weights[ i ] = w0[ i ] + ( w1[ i ] - w0[ i ] ) / eps * a
    //                                      + ( w0[ i ] + w2[ i ] - 2 * w1[ i ] ) / pow( eps, 2 ) * pow( a, 2 ) / 2
    //                                      + ( w3[ i ] - 3 * w2[ i ] + 3 * w1[ i ] - w0[ i ] ) / pow( eps, 3 ) * pow( a, 3 ) / 6;
    //     if ( si.l2_mass_error() < 1e10 )
    //         drrs.push_back( si.l2_mass_error() );
    // }

    // glot_vec_ys( arrs, brrs, crrs, drrs );
}
