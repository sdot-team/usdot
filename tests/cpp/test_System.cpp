#include <usdot/utility/glot.h>
#include <tl/support/P.h>

#include <usdot/System.h>
#include "catch_main.h"
#include <iostream>

// #include <usdot/utility/gmp.h>
// using namespace boost::multiprecision;
// using TF = number<backends::cpp_bin_float<64>>;

using namespace usdot;
using namespace std;
using TF = FP64;

TEST_CASE( "System", "" ) {
    // check_ders();
    std::vector<TF> dv( 100 );
    for( PI i = 0; i < dv.size(); ++i )
        dv[ i ] = max( 0.01, min( 1.0, 10 * ( i - dv.size() / 2.0 ) / dv.size() ) );

    glot_vec_ys( dv );
        
    DiffusionDensity<TF> gd( dv );

    System<TF> si;
    si.set_dirac_positions( cellspace<TF>( 0.0, 0.5, 10 ) );
    si.set_global_mass_ratio( 0.8 );
    si.set_density( &gd );

    si.stream = &std::cout;
    si.verbosity = 2;

    si.solve();
    si.plot();
}
