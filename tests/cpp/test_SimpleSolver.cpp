#include <partial1D/SimpleSolver.h>
#include "catch_main.h"
#include "partial1D/ThreadPool.h"
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
    // thread_pool.init( 1 ); 
    /// 
    // Time difference = 11183[µs]
    // Time difference = 2320[µs]
    // Time difference = 686[µs]

    // Time difference = 16159[µs]
    // Time difference = 5993[µs]
    // Time difference = 3001[µs]

    // Time difference = 2712[µs]

    SimpleSolverInput<TF> si;
    si.dirac_positions = Vec<TF>::linspace( 0.01, 0.99, 800000 );
    si.starting_contrast_ratio = 1e-6;
    si.target_contrast_ratio = 1e-6;
 
    // for( auto &p : si.dirac_positions )
    //     p += 0.05 * TF( rand() ) / RAND_MAX;

    si.density_values = { 1, 0.01, 1 };
    si.beg_x_density = 0;
    si.end_x_density = 1;
 
    SimpleSolver<TF> solver( si );
    solver.multithread = false;

    // solver.sorted_dirac_weights = solver.sorted_dirac_weights * 0 + 1.5;
    // solver.sorted_dirac_weights[ 5 ] = 4;
    // P( solver.sorted_dirac_weights );

    // P( solver.newton_system_ap() );
    // P( solver.newton_system() );

    // solver.solve();
    // solver.plot();

    // P( solver.normalized_density_x_integral( -1.1, 6.1 ) );

    // TF s = 0;
    // for( TF x : Vec<TF>::linspace( -1.1, 6.1, 1000000 ) )
    //     s += x * solver.normalized_density_value( x );
    // P( s * 7.2 / 1000000 );

    // glot( Vec<TF>::linspace( -3, 7, 1500 ),
    //     [&]( TF x ) { return solver.normalized_density_x_primitive( x ); }
    //     // [&]( TF x ) { return solver.normalized_density_integral( x, x + 1e-6 ) / 1e-6 + 1; },
    //     // [&]( TF x ) { return solver.normalized_density_value( x ); }
    //     // [&]( TF x ) { return solver.normalized_dirac_convolution( x, 5e-2 ); }
    // );

    // TF s = 0;
    // for( TF x : Vec<TF>::linspace( -10, 20, 100000 ) )
    //     s += solver.normalized_density_value( x );
    // P( s * 30 / 100000, solver.normalized_density_integral( -10, 20 ) );
    // P( solver.normalized_density_integral( 0, 3 ) );
    // A x + B^T y = a
    // B x + C   y = b
    // x = A^-1 * ( a - B^T y ) 
    // B * A^-1 * ( a - B y ) + C y = b
    // B * A^-1 * a - B * A^-1 * B y + C y = b
    //   ( C - B * A^-1 * B ) y = b - B * A^-1 * a
    {
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        solver.newton_dir();
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;
    }

    {
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        solver.newton_system();
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;
    }

    {
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        solver.errors();
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;
    }
}
