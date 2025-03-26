#include <partial1D/GradientSolver.h>
#include <partial1D/glot.h>
#include "catch_main.h"

using namespace usdot;
using namespace std;
using TF = double;

TEST_CASE( "Conv Grid solver", "" ) {
    GradientSolverInput<TF> si;
    si.dirac_positions = Vec<TF>::cellspace( 0, 1, 200 );
    si.density_values = { 1, 0, 1 };
 
    GradientSolver<TF> solver( std::move( si ) );
    solver.solve();

    // glot( Vec<TF>::linspace( -7, 12, 1000 ), 
    //     [&]( TF x ) { return solver.density_value( x ); }
    // );
}
