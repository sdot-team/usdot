#include <partial1D/GradientSolver.h>
#include <partial1D/glot.h>
#include "catch_main.h"

using namespace usdot;
using namespace std;
using TF = double;

TEST_CASE( "Conv Grid solver", "" ) {
    GradientSolverInput<TF> si;
    si.dirac_positions = Vec<TF>::cellspace( 0, 1, 50 );
    si.density_values = Vec<TF>::fill( 10, 1e-4 );
    si.density_values[ 4 ] = 1;
    si.density_values[ 5 ] = 1;
 
    GradientSolver<TF> solver( std::move( si ) );
    // solver.plot();
    solver.solve();

    // P( solver.cell_masses() );
    // P( solver.normalized_cell_boundaries() );

    // glot( Vec<TF>::linspace( -1, 2, 10000 ), 
    //     [&]( TF x ) { return solver.density_value( x, 0 ); },
    //     [&]( TF x ) { return solver.density_value( x, 1 ); },
    //     [&]( TF x ) { return solver.density_value( x, 2 ); }
    // );
}
