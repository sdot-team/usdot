#include <tl/support/operators/norm_2.h>
#include <tl/support/operators/sum.h>
#include <tl/support/ASSERT.h>

#include <partial1D/InitialSolution.h>

#include <matplotlibcpp.h>
#include <vector>

#include "catch_main.h"

using TF = PowerDiagram::TF;

// T_T std::vector<T> as_std_vec( const Vec<T> &t ) {
//     return { t.begin(), t.end() };
// }

// void check_initial_solution( Vec<TF> density_coords, Vec<TF> density_values, Vec<TF> dirac_coords, TF mass_ratio, Vec<TF> exp_centroids ) {
//     InitialSolution is( dirac_coords, density_coords, density_values, mass_ratio );
//     P( is.cell_ps() );
// }

// TEST_CASE( "Init sol", "" ) {
//     // check_initial_solution( { 0, 1 }, { 1 }, { 0.5      }, 0.5, {} );
//     // check_initial_solution( { 0, 1 }, { 1 }, { 0.2, 0.8 }, 0.4, {} );
//     // check_initial_solution( { 0, 1 }, { 1 }, { 0.1, 0.9 }, 0.6, {} );
//     // check_initial_solution( { 0, 1 }, { 1 }, { 0.1, 0.9 }, { 0.5 } );
//     // check_initial_solution( { 0, 1 }, { 1 }, { 0.6, 0.9 }, 0.6, {} );
//     // check_initial_solution( { 0, 1 }, { 1 }, { 0.1 }, 0.6, {} );
//     // check_initial_solution( { 0, 1 }, { 1 }, { 0.1, 0.3 }, 0.6, {} );
//     // check_initial_solution( { 0, 1 }, { 1 }, { 0.1, 0.4 }, 0.6, {} );
//     // check_initial_solution( { 0, 1 }, { 1 }, { 0.8, 0.9 }, 0.6, {} );
//     // check_initial_solution( { 0, 1 }, { 1 }, { 0.45, 0.55 }, 0.6, {} );
//     check_initial_solution( { 0, 0.4, 0.6, 1 }, { 1, 1e-5, 1 }, { 0.45, 0.55 }, 0.6, {} );
// }

// TEST_CASE( "Init sol gaussian", "" ) {
//     Vec<TF> density_coords;
//     Vec<TF> density_values;
//     for( TF x = -5; x < 5; x += 0.01 ) {
//         density_values << exp( - 6 * x * x );
//         density_coords << x;
//     }
//     density_values.pop_back();
//     // Vec<TF> density_coords{ -5, 5 };
//     // Vec<TF> density_values{ 1 };

//     Vec<TF> dirac_coords;
//     for( TF i = 0; i < 100; ++i )
//         dirac_coords << ( TF( rand() ) / RAND_MAX - 0.5 ) * 5;
//     std::sort( dirac_coords.begin(), dirac_coords.end() );

//     InitialSolution is( dirac_coords, density_coords, density_values, 1.0 );

//     std::vector<TF> xs, ys;
//     for( auto p : is.barycenters() ) {
//         xs.push_back( p );
//         ys.push_back( 1 );
//     }
//     matplotlibcpp::plot( xs, ys, "." );
//     matplotlibcpp::show();
// }

struct IsoResult {
    Vec<TF> beg_p;
    Vec<TF> end_p;
    Vec<TF> beg_x;
    Vec<TF> end_x;
    Vec<TF> bar_x;
};

void test_is( Vec<TF> density_coords, Vec<TF> density_values, Vec<TF> dirac_coords, TF ratio, IsoResult exp ) {
    InitialSolution inso = make_initial_solution( dirac_coords, density_coords, density_values, ratio );

    IsoResult res;
    inso.for_each_cell( [&]( const InitialSolution::Cell &cell ) {
        res.beg_p << cell.beg_p;
        res.end_p << cell.end_p;
        res.beg_x << cell.beg_x;
        res.end_x << cell.end_x;
        res.bar_x << inso.barycenter( cell );
    } );

    CHECK_PROX_IF_NOT_EMPTY( res.beg_p, exp.beg_p, 1e-6 );
    CHECK_PROX_IF_NOT_EMPTY( res.end_p, exp.end_p, 1e-6 );
    CHECK_PROX_IF_NOT_EMPTY( res.beg_x, exp.beg_x, 1e-6 );
    CHECK_PROX_IF_NOT_EMPTY( res.end_x, exp.end_x, 1e-6 );
    CHECK_PROX_IF_NOT_EMPTY( res.bar_x, exp.bar_x, 1e-6 );
}

// TEST_CASE( "Init sol basic", "" ) {
//     // single cell, hole in the density
//     test_is( { 0, 1, 2, 3 }, { 1, 0, 1 }, { 0.00 }, 0.2, { .beg_x = { 0.0 }, .end_x = { 0.4 }, .bar_x = { 0.2 } } );
//     test_is( { 0, 1, 2, 3 }, { 1, 0, 1 }, { 0.50 }, 0.2, { .beg_x = { 0.3 }, .end_x = { 0.7 }, .bar_x = { 0.5 } } );
//     test_is( { 0, 1, 2, 3 }, { 1, 0, 1 }, { 1.00 }, 0.2, { .bar_x = { 0.8 } } );
//     test_is( { 0, 1, 2, 3 }, { 1, 0, 1 }, { 1.50 }, 0.2, { .bar_x = { 1.5 } } );
//     test_is( { 0, 1, 2, 3 }, { 1, 0, 1 }, { 2.00 }, 0.2, { .bar_x = { 2.2 } } );
//     test_is( { 0, 1, 2, 3 }, { 1, 0, 1 }, { 2.10 }, 0.2, { .bar_x = { 2.2 } } );
//     test_is( { 0, 1, 2, 3 }, { 1, 0, 1 }, { 2.50 }, 0.2, { .bar_x = { 2.5 } } );

//     // several cells
//     test_is( { 0, 1 }, { 1 }, { 0.5, 0.5 }, 0.2, { .bar_x = { 0.4, 0.6 } } );
// }

TEST_CASE( "Init sol gaussian", "" ) {
    Vec<TF> density_coords;
    Vec<TF> density_values;
    for( TF s = 10.0 / 900, x = -2; x < 2 + s / 2; x += s ) {
        density_values << exp( - 1.5 * pow( x + s / 2, 2 ) );
        density_coords << x;
    }
    density_values.pop_back();

    Vec<TF> dirac_coords;
    for( TF i = 0; i < 100; ++i )
        dirac_coords << ( TF( rand() ) / RAND_MAX - 0.5 ) * 3;
    // std::sort( dirac_coords.begin(), dirac_coords.end() );
    // for( TF i = 0; i < 2; ++i )
    //     dirac_coords << 0;

    std::vector<std::vector<TF>> xs( 2 * dirac_coords.size() );
    std::vector<std::vector<TF>> ys( 2 * dirac_coords.size() );
    for( TF ratio = 0.01; ratio <= 1; ratio += 0.01 ) {
        InitialSolution inso = make_initial_solution( dirac_coords, density_coords, density_values, ratio );
        PI cpt = 0;
        inso.for_each_cell( [&]( CellInInitialSolution &cell ) {
            xs[ cpt + 0 ].push_back( cell.beg_x );
            xs[ cpt + 1 ].push_back( cell.end_x );
            ys[ cpt + 0 ].push_back( ratio );
            ys[ cpt + 1 ].push_back( ratio );
            cpt += 2;
        } );
    }
    
    for( PI i = 0; i < xs.size(); ++i )
        matplotlibcpp::plot( xs[ i ], ys[ i ], "-" );
    matplotlibcpp::show();
}
