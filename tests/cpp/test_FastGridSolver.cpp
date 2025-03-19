#include <partial1D/FastGridSolver.h>
#include <partial1D/SimpleSolver.h>
#include "catch_main.h"
#include "glot.h"

#include <chrono>

using namespace usdot;
using namespace std;
using TF = double;

// TEST_CASE( "Simple solver", "" ) {
//     // solver.newton_dir(); -> 6783[µs] x 2.6 vs 16159[µs] avant
//     // solver.l2_error();   -> 2767[µs] x 1   vs 3001 [µs] avant

//     SimpleSolverInput<TF> si;
//     si.dirac_positions = Vec<TF>::cellspace( 0, 0.5, 8 );
//     si.starting_contrast_ratio = 1e-5;
//     si.target_contrast_ratio = 1e-5;

//     si.density_values = { 1, 1, 1 };
//     si.beg_x_density = 0;
//     si.end_x_density = 1;
 
//     SimpleSolver<TF> simple_solver( std::move( si ) );
//     // P( simple_solver.newton_dir() );
//     simple_solver.solve();
//     // solver.plot();

//     // {
//     //     std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
//     //     solver.newton_dir();
//     //     std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
//     //     std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;
//     // }

//     // {
//     //     std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
//     //     solver.l2_error();
//     //     std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
//     //     std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;
//     // }

//     // P( solver.newton_dir() );


//     // glot( Vec<TF>::linspace( -6, 10, 1500 ),
//     //     // [&]( TF x ) { return solver.normalized_density_integral( x, x + 1e-6 ) / 1e-6; },
//     //     // [&]( TF x ) { return solver.normalized_density_value( x ); }
//     //     // [&]( TF x ) { return solver.normalized_dirac_convolution( x, 5e-2 ); }
//     //     [&]( TF x ) { return solver.normalized_density_x_integral_ap( 1, x ); },
//     //     [&]( TF x ) { return solver.normalized_density_x_integral( 1, x ) + 0.1; }
//     // );
// }

TEST_CASE( "Fast grid solver", "" ) {
    // solver.newton_dir(); -> 6783[µs] x 2.6 vs 16159[µs] avant
    // solver.l2_error();   -> 2767[µs] x 1   vs 3001 [µs] avant

    FastGridSolverInput<TF> si;
    si.dirac_positions = Vec<TF>::cellspace( 0, 0.75, 8 );
    si.starting_contrast_ratio = 0;
    si.target_contrast_ratio = 0;

    si.density_values = { 1 };
    si.beg_x_density = 0;
    si.end_x_density = 1;
 
    FastGridSolver<TF> solver( std::move( si ) );
    // P( solver.newton_dir() );
    solver.solve();
    // solver.plot();

    // {
    //     std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    //     solver.newton_dir();
    //     std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    //     std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;
    // }

    // {
    //     std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    //     solver.l2_error();
    //     std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    //     std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;
    // }

    // P( solver.newton_dir() );


    // glot( Vec<TF>::linspace( -6, 10, 1500 ),
    //     // [&]( TF x ) { return solver.normalized_density_integral( x, x + 1e-6 ) / 1e-6; },
    //     // [&]( TF x ) { return solver.normalized_density_value( x ); }
    //     // [&]( TF x ) { return solver.normalized_dirac_convolution( x, 5e-2 ); }
    //     [&]( TF x ) { return solver.normalized_density_x_integral_ap( 1, x ); },
    //     [&]( TF x ) { return solver.normalized_density_x_integral( 1, x ) + 0.1; }
    // );

}
