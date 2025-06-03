#include <tl/support/string/to_string.h>
#include <usdot/utility/glot.h>
#include <tl/support/P.h>

#include <usdot/PowerDiagram.h>
#include "catch_main.h"

using namespace usdot;
using namespace std;

// #include <usdot/utility/gmp.h>
// using namespace boost::multiprecision;
// using TF = number<backends::cpp_bin_float<128>>;
using TF = FP64;

TEST_CASE( "PowerDiagram", "" ) {
    // DiffusionDensity<TF> gd( linspace<TF>( 10, 11, 3 ), std::vector<TF>{ 1, 1, 1 } );
    PowerDiagram<TF> gd( dp, dv );
    // auto diracs = cellspace<TF>( 0.0, 1.0, 10 );
    // auto diracs = cellspace<TF>( 0.0, 1.0, 10 );

    System<TF> si;
    si.set_dirac_positions( pd, 1e-2 );
    si.set_global_mass_ratio( 0.984 );
    si.set_density( &gd );

    // si.stream = &std::cout;
    // si.verbosity = 2;

    auto t0 = std::chrono::high_resolution_clock::now();
    si.solve();
    auto t1 = std::chrono::high_resolution_clock::now();
    PE( std::chrono::duration<double>{ t1 - t0 } );

    si.plot();

    // std::vector<TF> f_grad( si.nb_original_diracs() );
    // for( PI i = 0; i < f_grad.size(); ++i )
    //     f_grad[ i ] = 2 * si.global_mass_ratio / si.nb_original_diracs() * ( si.dirac_positions()[ i ] - si.cell_barycenters()[ i ] );
       

    // const TF eps = 1e-20;
    // const TF cost_ref = pow( si.cost(), 2 );
    // for( PI i = 0; i < si.nb_original_diracs(); ++i ) {
    //     auto ov = diracs[ i ];
    //     diracs[ i ] += eps;

    //     System<TF> si;
    //     si.target_newton_error = 1e-40;
    //     si.set_dirac_positions( diracs );
    //     si.set_global_mass_ratio( 0.1 );
    //     si.set_density( &gd );
    //     si.solve();

    //     P( ( pow( si.cost(), 2 ) - cost_ref ) / eps - f_grad[ i ] );

    //     diracs[ i ] = ov;
    // }
        

    // si.plot();
}
