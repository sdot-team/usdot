#include <tl/support/operators/norm_2.h>
#include <matplotlibcpp.h>
#include "catch_main.h"
#include "Proj.h"

using TF = PowerDiagram::TF;

TEST_CASE( "Flower ", "" ) {
    Proj proj;

    TF dev        = 1e-2;
    PI nb_samples = 10000;
    PI nb_lobes   = 5;
    for( PI n = 0; n < nb_samples; ++n ) {
        TF a = 2 * std::numbers::pi * n / nb_samples;
        TF r = 2 + sin( nb_lobes * a );
        proj.points << Proj::Pt{
            r * cos( a ) + ( TF( rand() ) / RAND_MAX - 0.5 ) * dev,
            r * sin( a ) + ( TF( rand() ) / RAND_MAX - 0.5 ) * dev
        };
    }


    // dirac_points;
    for( PI n = 0; n < 100; ++n ) {
        proj.dirac_points << Proj::Pt{
            TF( rand() ) / RAND_MAX + .5,
            TF( rand() ) / RAND_MAX
        };
    }


    proj.plot_points( "results/points.vtk" );
    proj.plot_diracs( "results/diracs_0.vtk" );
    for( PI num_step = 1; num_step <= 130; ++num_step ) {
        P( num_step );
        proj.quant_step();
        proj.plot_diracs( std::format( "results/diracs_{}.vtk", num_step ) );
    }
    // auto d = f.make_density( { 1, 0 }, 300 );
    // matplotlibcpp::plot( as_std_vec( d.first ), as_std_vec( d.second ), {} );
    // matplotlibcpp::show();
}
