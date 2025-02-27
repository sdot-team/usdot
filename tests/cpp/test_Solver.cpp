#include <partial1D/Density/PiecewiseConstantDensity.h>
#include <partial1D/Convolution/GaussianConvolution.h>
#include <partial1D/Solver.h>
#include "catch_main.h"
#include <fstream>

using namespace usdot;
using namespace std;
using TF = double;

void glot_bnds_evolution( const Vec<Vec<Vec<TF>>> &vbnds, TF mi, TF ma ) {
    std::ofstream fs( "glot.py" );
    fs << "from matplotlib import pyplot\n";
    for( const Vec<Vec<TF>> &bnds : vbnds ) {
        Vec<FP64> ys;
        for( PI nw = 0; nw < bnds.size(); ++nw )
            ys << nw;
        for( PI nd = 0; nd < bnds[ 0 ].size(); ++nd ) {
            Vec<FP64> xs;
            for( PI nw = 0; nw < bnds.size(); ++nw )
                xs << min( ma, max( mi, bnds[ nw ][ nd ] ) );
            fs << "pyplot.plot( " << to_string( xs ) << ", " << to_string( ys ) << " )\n";
        }
    }
    fs << "pyplot.show()\n";
}

// void disp_conv_hull() {
//     PI nb_diracs = 31;

//     auto cb = RcPtr<PiecewiseConstantDensity<TF>>::New( 
//         Vec<TF>::linspace( 0.0, 1.0, 4 ),
//         Vec<TF>{ 1, 0, 1 }
//     );

//     auto pd = RcPtr<PowerDiagram<TF>>::New(
//         Vec<TF>::cellspace( 0, 1, nb_diracs ),
//         Vec<TF>::fill( nb_diracs, 1 )
//     );

//     Solver<TF> solver( pd, cb );

//     // Rq: si le domaine est borné on peut limiter les différences de potentiels
//     // De façon générale, on pourrait chercher à limiter les contrastes de densité
//     // par exemple, on pourrait dire que si le système conduit à augmenter 

//     // 
//     Vec<Vec<TF>> bnds;
//     solver.convex_hull_density_ratio = 1;
//     solver.update_weights();
//     auto s = solver.get_state();

//     for( TF r = 1; r >= 5e-2; r *= 0.9 ) {
//         solver.set_state( s );
//         solver.update_convex_hull_density_ratio( { .target_value = r, .polynomial_order = 2 } );

//         bnds << pd->cell_boundaries();
//     }

//     glot_bnds_evolution( { bnds }, 0, 1 );
// }

TEST_CASE( "Solver prog", "" ) {
    PI nb_diracs = 31;

    auto cb = RcPtr<PiecewiseConstantDensity<TF>>::New( 
        Vec<TF>::linspace( 0.0, 1.0, 4 ),
        Vec<TF>{ 1, 0, 1 }
    );

    auto pd = RcPtr<PowerDiagram<TF>>::New(
        Vec<TF>::linspace( 0, 1, nb_diracs ),
        Vec<TF>::fill( nb_diracs, 1 )
    );

    Solver<TF> solver( pd, cb );
    solver.convolution_factory = []( TF width ) -> RcPtr<Convolution<TF>> { return new GaussianConvolution<TF>( width ); };
    // solver.convolution_factory = []( TF width ) -> RcPtr<Convolution<TF>> { return new InvX2Convolution<TF>( width ); };

    // 
    Vec<Vec<TF>> bnds;
    for( TF w : Vec<TF>::linspace( 0.2, 1e-2, 100 ) ) {
        solver.convolution_width = w;
        solver.update_weights();
        bnds << pd->cell_boundaries();
    }

    // auto s = solver.get_state();

    // for( TF r = 1; r >= 5e-2; r *= 0.9 ) {
    //     solver.set_state( s );
    //     solver.update_convex_hull_density_ratio( { .target_value = r, .polynomial_order = 2 } );

    //     bnds << pd->cell_boundaries();
    // }

    glot_bnds_evolution( { bnds }, 0, 1 );
}
