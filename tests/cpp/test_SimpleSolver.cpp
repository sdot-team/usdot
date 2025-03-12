#include <partial1D/SimpleSolver.h>
#include "catch_main.h"
#include <fstream>

using namespace usdot;
using namespace std;
using TF = double;

void glot( Vec<TF> xs, auto &&...funcs ) {
    std::ofstream fs( "glot.py" );
    fs << "from matplotlib import pyplot\n";
    auto pf = [&]( auto &&func ) {
        Vec<FP64> ys;
        for( auto x : xs )
            ys << func( x );
        fs << "pyplot.plot( " << to_string( xs ) << ", " << to_string( ys ) << " )\n";
    };
    ( pf( funcs ), ... );
    fs << "pyplot.show()\n";
}
 
TEST_CASE( "Simple solver prog", "" ) {
    SimpleSolverInput<TF> si;
    si.dirac_positions = Vec<TF>::linspace( -0.1, 1.1, 10 );
    si.density_values = { 1, 0, 1 };
    si.starting_contrast_ratio = .1;
    si.beg_x_density = 0;
    si.end_x_density = 1;
 
    SimpleSolver<TF> solver( si );

    glot( Vec<TF>::linspace( -1, 5, 500 ),
        [&]( TF x ) { return solver.normalized_density_value( x ); },
        [&]( TF x ) { return solver.normalized_dirac_convolution( x ); }
    );
}
