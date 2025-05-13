#include <tl/support/string/to_string.h>
#include <usdot/utility/glot.h>
#include <tl/support/P.h>

#include <usdot/System.h>
#include "catch_main.h"
#include <iostream>

// #include <usdot/utility/gmp.h>
// using namespace boost::multiprecision;
// using TF = number<backends::cpp_bin_float<64>>;

using namespace usdot;
using namespace std;
using TF = FP64;

// TEST_CASE( "System", "" ) {
//     // check_ders();
//     std::vector<TF> dv( 1000 );
//     for( PI i = 0; i < dv.size(); ++i ) {
//         TF f = -1 + 2 * TF( i ) / ( dv.size() - 1 );
//         dv[ i ] = max( 0.1, min( 1.0, 10 * ( abs( f ) - .3 ) ) );
//     }

//     // glot_vec_ys( dv );
        
//     DiffusionDensity<TF> gd( dv );

//     // glot_stream( [&]( std::ostream &fs ) {
//     //     for ( TF r : linspace<TF>( 0, 1, 10 ) ) {
//     //         gd.set_flattening_ratio( r );
//     //         gd.plot( fs, r == 0 ? "-" : "--", r == 0 || r == 1 ? 2 : 1 );
//     //     }
//     // } );

//     System<TF> si;
//     si.set_dirac_positions( cellspace<TF>( 0.0, 1.0, 10 ) );
//     si.set_global_mass_ratio( 1.0 );
//     si.set_density( &gd );

//     si.stream = &std::cout;
//     si.verbosity = 2;

//     si.solve();
//     // si.plot();
// }
TEST_CASE( "illustration diff", "" ) {
    // check_ders();
    std::vector<TF> dv( 1000 );
    for( PI i = 0; i < dv.size(); ++i ) {
        TF f = -1 + 2 * TF( i ) / ( dv.size() - 1 );
        dv[ i ] = max( 0.15, min( 1.0, 10 * ( abs( f ) - .3 ) ) );
        dv[ i ] += max( -0.1, min( 0.1, sin( 14 * f ) ) );
    }

    // glot_vec_ys( dv );
        
    DiffusionDensity<TF> gd( dv );

    // glot_stream( [&]( std::ostream &fs ) {
    //     for ( TF r : linspace<TF>( 0, 1, 14 ) ) {
    //         gd.set_flattening_ratio( r );
    //         gd.plot( fs, r == 0 ? "-" : "--", r == 0 || r == 1 ? 2 : 1 );
    //     }
    // } );

    System<TF> si;
    si.set_dirac_positions( cellspace<TF>( 0.0, 1.0, 10 ) );
    si.set_global_mass_ratio( 1.0 );
    si.set_density( &gd );
    si.stream = &std::cout;
    si.verbosity = 2;

    si.initialize_with_flat_density();
    gd.set_flattening_ratio( 1 - 1e-6 );
    si.newton_iterations( 1e-3 );
    std::vector<std::vector<std::array<TF,2>>> bnds;
    std::vector<TF> ys;
    PI cpt_y = 0;
    for( TF d : { 0.5, 0.25, 0.125, 0.0 } ) {
        auto w0 = si.sorted_dirac_weights;
        auto ders = si.der_weights_wrt_flat_ratio( 2 );
        const TF a = d - gd.current_flattening_ratio;
 
        gd.set_flattening_ratio( d );

        for( TF la : linspace<TF>( 0, a, 10 ) ) {
            for( PI i = 0; i < si.nb_sorted_diracs(); ++i )
                si.sorted_dirac_weights[ i ] = w0[ i ] + ders[ 0 ][ i ] * la + ders[ 1 ][ i ] * pow( la, 2 ) / 2;
            bnds.push_back( si.cell_boundaries() );
            ys.push_back( cpt_y++ * 4 / 36.0 );
        }
        cpt_y--;

        si.newton_iterations( 1 );

    }
    // bnds.push_back( si.cell_boundaries() );
    // ys.push_back( cpt_y++ );
    // si.plot_bnds_evolution( bnds );
    P( ys );

    std::ofstream fs( "glot.py" );
    fs << "from matplotlib import pyplot\n";
    for( PI i = 0; i < bnds[ 0 ].size(); ++i ) {
        for( PI j = 0; j < bnds[ 0 ][ 0 ].size(); ++j ) {
            std::vector<TF> xs( bnds.size() );
            for( PI n = 0; n < bnds.size(); ++n )
                xs[ n ] = bnds[ n ][ i ][ j ];
            fs << "pyplot.plot( " << to_string( xs ) << ", " << to_string( ys ) << " )\n";
        }
    }
    fs << "pyplot.ylabel( 'step number (modification of t)' )\n";
    fs << "pyplot.xlabel( 'cell boundaries' )\n";
    fs << "pyplot.legend()\n";
    fs << "pyplot.show()\n";
    

    // si.solve();
    // si.plot();
}
