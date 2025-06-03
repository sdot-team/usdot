#include <tl/support/string/to_string.h>
#include <usdot/utility/glot.h>
#include <tl/support/P.h>
#include "catch_main.h"
#include "usdot/utility/linspace.h"

#include <usdot/SdotSystem.h>

// #include <usdot/utility/gmp.h>
// using namespace boost::multiprecision;
// using TF = number<backends::cpp_bin_float<64>>;
using TF = usdot::FP64;

using namespace usdot;
using namespace std;

#include "problemos.h"

// TEST_CASE( "SdotSystem", "" ) {
//     DiffusionDensity<TF> density( { 0, 1, 0 } );
//     SdotSystem<TF> sd( &density, { 0.4, 0.5, 0.6 } ); // , 0.6, 1.0
//     sd.sorted_dirac_weights[ 0 ] = 0.1;
//     sd.epsilon = 0.29;
    
//     TridiagonalSymmetricMatrix<TF> Ma;
//     std::vector<TF> Va;
//     PI nb_arcs;
//     sd.get_sorted_newton_system_ap( Ma, Va, nb_arcs, 1e-6 );
//     P( Ma, Va );
    
//     TridiagonalSymmetricMatrix<TF> M;
//     std::vector<TF> V;
//     sd.get_sorted_newton_system( M, V, nb_arcs );
//     P( M, V );

//     sd.plot();

//     // P( sd.sorted_cell_masses() );
// }

TEST_CASE( "SdotSystem", "" ) {
    DiffusionDensity<TF> density( dp, dv );
    SdotSystem<TF> sd( &density, pd, 0.984 ); // , 0.6, 1.0
    // DiffusionDensity<TF> density( { 0, 1, 0 } );
    // SdotSystem<TF> sd( &density, cellspace<TF>( -1, 1, 3 ), 0.984 ); // , 0.6, 1.0
    sd.verbosity = 2;
    for( TF epsilon : { 10. } ) {
        sd.epsilon = epsilon;
        sd.initialize_weights();
        sd.plot();
        sd.solve();
    }


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
// ^[[A^[[A  relax: 0.015625 nb_newton_iterations: 22 max_error_ratio: 0.945417
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