#include <tl/support/P.h>
#include <usdot/utility/glot.h>

#include <usdot/System.h>
#include "catch_main.h"
#include <iostream>

// #include <usdot/utility/gmp.h>
// using namespace boost::multiprecision;
// using TF = number<backends::cpp_bin_float<40>>;

using namespace usdot;
using namespace std;
using TF = double;

TEST_CASE( "System", "" ) {
    GridDensity<TF,2> gd( { 1, 0.01, 1 } );

    System<TF,GridDensity<TF,2>> si;
    si.stream = &std::cout;
    si.verbosity = 2;

    // si.set_dirac_positions( Vec<TF>::cellspace( 1.6, 2.0, 15 ) );
    // si.set_dirac_positions( cellspace<TF>( 0.0, 2.0, 5 ) );
    si.set_dirac_positions( { -10, -11, -12, 4, 0, 1, 2 } );
    si.set_global_mass_ratio( 0.5 );
    si.set_density( &gd );

    si.solve();

    si.plot();
    P( si.l2_mass_error() );
}
