#include <partial1D/ConvGridSolver.h>
#include <partial1D/glot.h>
#include "catch_main.h"

using namespace usdot;
using namespace std;
using TF = double;

TEST_CASE( "Conv Grid solver", "" ) {
    ConvGridSolverInput<TF> si;
    si.dirac_positions = Vec<TF>::cellspace( 0, 1, 50 );
    si.starting_filter_value = 0.95; //99;
    si.target_filter_value = 0.01;

    si.density_values = { 1, 1, 1 };
    si.beg_x_density = 0;
    si.end_x_density = 1;
 
    ConvGridSolver<TF> solver( std::move( si ) );
    solver.solve();
    // solver.glot();

    glot( Vec<TF>::linspace( -7, 12, 1000 ), 
        [&]( TF x ) { return solver.density_value( x ); }
    );

}
