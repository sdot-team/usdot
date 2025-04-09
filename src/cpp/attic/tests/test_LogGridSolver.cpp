#include <partial1D/LogGridSolver.h>
#include "catch_main.h"

using namespace usdot;
using namespace std;
using TF = double;

TEST_CASE( "Log grid solver", "" ) {
    // solver.newton_dir(); -> 6783[µs] x 2.6 vs 16159[µs] avant
    // solver.l2_error();   -> 2767[µs] x 1   vs 3001 [µs] avant

    LogGridSolverInput<TF> si;
    si.dirac_positions = Vec<TF>::cellspace( 0, 0.75, 10000 ); // 30
    si.starting_contrast_ratio = 0.5;
    si.target_contrast_ratio = 1e-6;

    si.density_values = { 1, 0.2, 1 };
    si.beg_x_density = 0;
    si.end_x_density = 1;
 
    LogGridSolver<TF> solver( std::move( si ) );
    solver.solve();

    solver.plot();

    // try {
    //     solver.update_weights( 1 );
    // } catch ( std::runtime_error ) {}
    // auto o_w = solver.sorted_dirac_weights;
    // auto dir = solver.newton_dir();
    
    // Vec<Vec<TF>> traj( FromSize(), 2 * solver.nb_diracs() );
    // auto as = Vec<TF>::linspace( 0, 1, 30 );
    // for( TF a : as ) {
    //     solver.sorted_dirac_weights = o_w + a * dir;

    //     auto bnds = solver.normalized_cell_boundaries();
    //     for( PI i = 0; i < solver.nb_diracs(); ++i ) {
    //         traj[ 2 * i + 0 ] << bnds[ i ][ 0 ];
    //         traj[ 2 * i + 1 ] << bnds[ i ][ 1 ];
    //     }            
    // }
    // std::ofstream fs( "glot.py" );
    // fs << "from matplotlib import pyplot\n";
    // for( auto t : traj )
    //     fs << "pyplot.plot( " << to_string( t ) << ", " << to_string( as ) << " )\n";
    // fs << "pyplot.show()\n";


    // try {
    //     P( solver.newton_dir() );
    //     solver.solve();
    // } catch ( std::runtime_error ) {

    // }

    // P( solver.cell_barycenters() );
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
