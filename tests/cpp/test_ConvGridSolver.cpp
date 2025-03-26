#include <partial1D/ConvGridSolver.h>
#include <partial1D/glot.h>
#include "catch_main.h"

using namespace usdot;
using namespace std;
using TF = double;

TEST_CASE( "Conv Grid solver", "" ) {
    ConvGridSolverInput<TF> si;
    si.dirac_positions = Vec<TF>::cellspace( 0, 0.75, 500 );
    si.starting_filter_value = 0.5; //99;
    si.target_filter_value = 0.00001;

    si.density_values = { 0, 1, 1, 0, 0, 1, 1, 0 };
    si.beg_x_density = 0;
    si.end_x_density = 1;
 
    ConvGridSolver<TF> solver( std::move( si ) );
    solver.solve();
    // solver.plot();

    // glot( Vec<TF>::linspace( -7, 12, 1000 ), 
    //     [&]( TF x ) { return solver.density_value( x ); }
    // );
}
// TF cp( TF mi ) {
//     ConvGridSolverInput<TF> si;
//     si.dirac_positions = { 0, 0.5, 1 };
//     si.starting_filter_value = 1e-2;
//     si.target_filter_value = 0.0;

//     si.density_values = { mi, mi, 1, 1 };
//     si.beg_x_density = 0;
//     si.end_x_density = 1;
 
//     ConvGridSolver<TF> solver( std::move( si ) );
//     auto d = solver.newton_dir().first;
//     for( PI i = 0; i < solver.nb_diracs(); ++i )
//         solver.sorted_dirac_weights[ i ] += d[ i ];

//     auto bnds = solver.normalized_cell_boundaries();
//     P( d, bnds );

//     return bnds[ 0 ][ 1 ];
// }

// TEST_CASE( "Conv Grid solver", "" ) {
//     for( TF mi = 1; mi > 1e-3; mi /= 2 )
//         P( mi, cp( mi ) );
// }
