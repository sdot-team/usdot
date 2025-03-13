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
    si.dirac_positions = Vec<TF>::linspace( 0.01, 0.99, 100 );
    si.starting_contrast_ratio = 1e-1;

    si.density_values = { 1, 0.5, 1 };
    si.beg_x_density = 0;
    si.end_x_density = 1;
 
    SimpleSolver<TF> solver( si );

    glot( Vec<TF>::linspace( -3, 7, 1500 ),
        [&]( TF x ) { return solver.normalized_density_integral( -10, x ); }
        // [&]( TF x ) { return solver.normalized_density_integral( x, x + 1e-6 ) / 1e-6 + 1; },
        // [&]( TF x ) { return solver.normalized_density_value( x ); }
        // [&]( TF x ) { return solver.normalized_dirac_convolution( x, 5e-2 ); }
    );

    TF s = 0;
    for( TF x : Vec<TF>::linspace( -10, 20, 100000 ) )
        s += solver.normalized_density_value( x );
    P( s * 30 / 100000, solver.normalized_density_integral( -10, 20 ) );
    P( solver.normalized_density_integral( 0, 3 ) );

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
    // solver.solve();

}
