#include <partial1D/SimpleSolver.h>
#include "catch_main.h"
#include <fstream>
#include <chrono>

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
    si.dirac_positions = Vec<TF>::linspace( 0.01, 0.99, 1000 );
    si.starting_contrast_ratio = 1e-1;

    si.density_values = { 1, 0.95, 1 };
    si.beg_x_density = 0;
    si.end_x_density = 1;
 
    SimpleSolver<TF> solver( si );

    // {
    //     std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    //     solver.newton_dir();
    //     std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    //     std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;
    // }

    // {
    //     std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    //     solver.newton_system();
    //     std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    //     std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;
    // }

    // {
    //     std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    //     solver.normalized_error();
    //     std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    //     std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;
    // }

    // solver.update_weights();
    // P( solver.normalized_integrals() );
    solver.solve();

    // glot( Vec<TF>::linspace( -1, 4, 1500 ),
    //     [&]( TF x ) { return solver.normalized_density_value( x ); },
    //     [&]( TF x ) { return solver.normalized_dirac_convolution( x, 5e-2 ); }
    // );
}
