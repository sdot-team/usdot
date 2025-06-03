#include <tl/support/string/to_string.h>
#include <usdot/utility/linspace.h>
#include <usdot/utility/glot.h>
#include <tl/support/P.h>
#include "catch_main.h"

#include <usdot/UsdotSystem.h>

// #include <usdot/utility/gmp.h>
// using namespace boost::multiprecision;
// using TF = number<backends::cpp_bin_float<64>>;
using TF = usdot::FP64;

using namespace usdot;
using namespace std;

#include "problemos.h"

TEST_CASE( "SdotSystem", "" ) {
    std::vector<TF> my_dv = problemos_dv;
    // for( PI i = 0; i + 1 < problemos_dv.size(); ++i )
    //     for( PI j = 0, n = 500; j < n; ++j )
    //         my_dv.push_back( problemos_dv[ i ] + ( problemos_dv[ i + 1 ] - problemos_dv[ i ] ) * j / n );
    // my_dv.push_back( problemos_dv.back() );

    UsdotDensity<TF> density(  -33.8823, 42.7147, my_dv );

    UsdotSystem<TF> sd( &density, problemos_pd, 0.984 ); // , 0.6, 1.0
    // UsdotDensity<TF> density( 0, 1, { 0.1, 1.0, 0.1 } );
    // UsdotSystem<TF> sd( &density, cellspace<TF>( -1, 2, 15 ), 0.84 ); // , 0.6, 1.0
    // sd.verbosity = 2;
    auto t0 = std::chrono::high_resolution_clock::now();

    sd.solve();

    auto t1 = std::chrono::high_resolution_clock::now();
    PE( std::chrono::duration<double>{ t1 - t0 } );
    PE( sd.nb_newton_iterations );

    // P( sd.original_cell_barycenters() );
    sd.plot();


    // P( sd.sorted_cell_masses() );
}
// check_newton_system: ok
//   relax: 0.00195312 nb_newton_iterations: 0 max_error_ratio: 1.29731
//   relax: 0.00195312 nb_newton_iterations: 1 max_error_ratio: 1.29478
//   relax: 0.0078125 nb_newton_iterations: 2 max_error_ratio: 1.28467
//   relax: 0.0078125 nb_newton_iterations: 3 max_error_ratio: 1.27464
//   relax: 0.015625 nb_newton_iterations: 4 max_error_ratio: 1.25474
//   relax: 0.015625 nb_newton_iterations: 5 max_error_ratio: 1.23515
//   relax: 0.015625 nb_newton_iterations: 6 max_error_ratio: 1.21587
//   relax: 0.015625 nb_newton_iterations: 7 max_error_ratio: 1.19688
//   relax: 0.015625 nb_newton_iterations: 8 max_error_ratio: 1.18012
//   relax: 0.015625 nb_newton_iterations: 9 max_error_ratio: 1.16172
//   relax: 0.015625 nb_newton_iterations: 10 max_error_ratio: 1.1436
//   relax: 0.015625 nb_newton_iterations: 11 max_error_ratio: 1.12576
//   relax: 0.015625 nb_newton_iterations: 12 max_error_ratio: 1.1082
//   relax: 0.015625 nb_newton_iterations: 13 max_error_ratio: 1.09092
//   relax: 0.015625 nb_newton_iterations: 14 max_error_ratio: 1.0739
//   relax: 0.015625 nb_newton_iterations: 15 max_error_ratio: 1.04683
//   relax: 0.015625 nb_newton_iterations: 16 max_error_ratio: 1.02984
//   relax: 0.015625 nb_newton_iterations: 17 max_error_ratio: 1.01307
//   relax: 0.015625 nb_newton_iterations: 18 max_error_ratio: 0.996511
//   relax: 0.015625 nb_newton_iterations: 19 max_error_ratio: 0.980153
//   relax: 0.015625 nb_newton_iterations: 20 max_error_ratio: 0.963987
//   relax: 0.015625 nb_newton_iterations: 21 max_error_ratio: 0.960093
//   relax: 0.015625 nb_newton_iterations: 22 max_error_ratio: 0.945417
//   relax: 0.015625 nb_newton_iterations: 23 max_error_ratio: 0.930603
//   relax: 0.015625 nb_newton_iterations: 24 max_error_ratio: 0.916018
//   relax: 0.015625 nb_newton_iterations: 25 max_error_ratio: 0.90166
//   relax: 0.03125 nb_newton_iterations: 26 max_error_ratio: 0.948855
//   relax: 0.03125 nb_newton_iterations: 27 max_error_ratio: 0.962473
//   relax: 0.03125 nb_newton_iterations: 28 max_error_ratio: 0.927127
//   relax: 0.03125 nb_newton_iterations: 29 max_error_ratio: 0.895344
//   relax: 0.03125 nb_newton_iterations: 30 max_error_ratio: 0.865518
//   relax: 0.03125 nb_newton_iterations: 31 max_error_ratio: 0.83712
//   relax: 0.03125 nb_newton_iterations: 32 max_error_ratio: 0.809916
//   relax: 0.03125 nb_newton_iterations: 33 max_error_ratio: 0.783766
//   relax: 0.03125 nb_newton_iterations: 34 max_error_ratio: 0.772071
//   relax: 0.03125 nb_newton_iterations: 35 max_error_ratio: 0.748099
//   relax: 0.03125 nb_newton_iterations: 36 max_error_ratio: 0.724869
//   relax: 0.03125 nb_newton_iterations: 37 max_error_ratio: 22.7431
//   relax: 0.0078125 nb_newton_iterations: 38 max_error_ratio: 22.5217
//   relax: 0.0078125 nb_newton_iterations: 39 max_error_ratio: 22.3109
//   relax: 1.52588e-05 nb_newton_iterations: 40 max_error_ratio: 22.3106