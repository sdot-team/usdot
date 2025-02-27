#include <tl/support/operators/norm_2.h>
#include <tl/support/operators/sum.h>
#include <tl/support/ASSERT.h>

#include <partial1D/InitialSolution.h>

#include <matplotlibcpp.h>
#include <vector>

#include "catch_main.h"

using TF = PowerDiagram::TF;

T_T std::vector<T> as_std_vec( const Vec<T> &t ) {
    return { t.begin(), t.end() };
}

// void check_initial_solution( Vec<TF> density_coords, Vec<TF> density_values, Vec<TF> dirac_coords, TF mass_ratio, Vec<TF> exp_centroids ) {
//     InitialSolution is( dirac_coords, density_coords, density_values, mass_ratio );
//     P( is.cell_ys() );
// }

// TEST_CASE( "Init sol", "" ) {
//     // check_initial_solution( { 0, 1 }, { 1, 1 }, { 0.5      }, 0.5, {} );
//     // check_initial_solution( { 0, 1 }, { 1, 1 }, { 0.2, 0.8 }, 0.4, {} );
//     // check_initial_solution( { 0, 1 }, { 1, 1 }, { 0.1, 0.9 }, 0.6, {} );
//     // check_initial_solution( { 0, 1 }, { 1, 1 }, { 0.1, 0.9 }, { 0.5 } );
//     // check_initial_solution( { 0, 1 }, { 1, 1 }, { 0.6, 0.9 }, 0.6, {} );
//     // check_initial_solution( { 0, 1 }, { 1, 1 }, { 0.1 }, 0.6, {} );
//     // check_initial_solution( { 0, 1 }, { 1, 1 }, { 0.1, 0.3 }, 0.6, {} );
//     // check_initial_solution( { 0, 1 }, { 1, 1 }, { 0.1, 0.4 }, 0.6, {} );
//     check_initial_solution( { 0, 1 }, { 1, 1 }, { 0.8, 0.9 }, 0.6, {} );


//     // matplotlibcpp::plot( as_std_vec( is.dirac_coords ), as_std_vec( is.dirac_y_offsets ), "+" );
//     // matplotlibcpp::show();
// }

TEST_CASE( "Init sol gaussian", "" ) {
    // Vec<TF> density_coords;
    // Vec<TF> density_values;
    // for( TF x = -5; x < 5; x += 0.01 ) {
    //     density_values << exp( - x * x );
    //     density_coords << x;
    // }
    Vec<TF> density_coords{ 0, 0.5, 1 };
    Vec<TF> density_values{ 0, 0.5, 1 };

    Vec<TF> dirac_coords;
    for( TF i = 0; i < 8; ++i )
        dirac_coords << ( TF( rand() ) / RAND_MAX - 0.5 ) * 2;
    // std::sort( dirac_coords.begin(), dirac_coords.end() );

    InitialSolution is( dirac_coords, density_coords, density_values, 1 );

    // std::vector<TF> xs, ys;
    // for( TF x = 0; x < 1; x += 0.01 ) {
    //     ys.push_back( is.density_primitive( x ) );
    //     xs.push_back( x );
    // }
    // matplotlibcpp::plot( xs, ys, "-" );
    // matplotlibcpp::show();
    // P( is.cell_ys() );
    // // P( is );

    std::vector<TF> xs, ys;
    for( auto p : is.barycenters() ) {
        xs.push_back( p );
        ys.push_back( 1 );
    }
    P( xs );
    matplotlibcpp::plot( xs, ys, "." );
    matplotlibcpp::show();

}

// TEST_CASE( "Init sol", "" ) {
//     Vec<TF> density_coords;
//     Vec<TF> density_values;
//     for( TF x = -5; x < 5; x += 0.1 ) {
//         density_values << exp( - x * x );
//         density_coords << x;
//     }

//     Vec<TF> dirac_coords;
//     for( TF i = 0; i < 150; ++i )
//         dirac_coords << ( TF( rand() ) / RAND_MAX - 0.5 ) * 2;

//     InitialSolution is( dirac_coords, density_coords, density_values, 0.5 );
//     is.solve();

//     // matplotlibcpp::plot( as_std_vec( is.dirac_coords ), as_std_vec( is.dirac_y_offsets ), "+" );
//     // matplotlibcpp::show();
// }
