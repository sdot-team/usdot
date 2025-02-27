#include <tl/support/operators/norm_2.h>
#include <tl/support/operators/sum.h>
#include <tl/support/ASSERT.h>
#include <partial1D/Solver.h>
#include "catch_main.h"

using TF = PowerDiagram::TF;

// TEST_CASE( "Solver linear", "" ) {
//     Vec<TF> density_coords = { -10, 10 };
//     Vec<TF> density_values = { 1, 1 };
//     Vec<TF> dirac_coords = { -1, 0, 1 };

//     Solver sr( dirac_coords, density_coords, density_values, 1 );
//     sr.solve();

//     VtkOutput vo;
//     sr.power_diagram().write_vtk( vo, true );
//     vo.save( "out.vtk" );
// }

TF my_norm_2( const Vec<TF> &a, const Vec<TF> &b ) {
    TF res = ( a[ 0 ] - b[ 0 ] ) * ( a[ 0 ] - b[ 0 ] );
    for( std::size_t i = 1; i < a.size(); ++i )
        res += ( a[ i ] - b[ i ] ) * ( a[ i ] - b[ i ] );
    return res;
}

TEST_CASE( "Solver ", "" ) {
    Vec<TF> density_coords;
    Vec<TF> density_values;
    for( TF x = -5; x < 5; x += 0.1 ) {
        density_values << exp( - x * x );
        density_coords << x;
    }

    Vec<TF> dirac_coords;
    for( TF i = 0; i < 15; ++i )
        dirac_coords << TF( rand() ) / RAND_MAX;


    Solver sr( dirac_coords, density_coords, density_values, 0.5 );
    sr.max_error_ratio_target = 1e-17;
    sr.min_newton_iterations = 3;
    sr.vtk_output = "out_$0.vtk";
    sr.verbosity = 0;
    sr.set_mass_ratio( 0.5 );

    sr.set_epsilon( 1e-6 );
    sr.solve();

    auto ref_weights = sr.dirac_weights;

    for ( TF epsilon = 2e-6; epsilon < 0.5; epsilon *= 2 ) {
        sr.set_epsilon( epsilon );
        sr.solve();
    
        // P( epsilon * sr.power_diagram().mass_of_density() / sum( sr.power_diagram().areas() ) );
        P( epsilon, my_norm_2( sr.dirac_weights, ref_weights ) );
    }
}
