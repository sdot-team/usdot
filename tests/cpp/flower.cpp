// #include <tl/support/operators/norm_2.h>
#include <cstdlib>
#include <tl/support/operators/abs.h>
#include <tl/support/operators/sp.h>
#include <eigen3/Eigen/LU>
#include <partial1D/Solver.h>
#include <matplotlibcpp.h>
#include "catch_main.h"
#include "eigen3/Eigen/src/Core/Matrix.h"
#include "eigen3/Eigen/src/LU/FullPivLU.h"
#include "partial1D/InitialSolution.h"
#include <numbers>
#include <limits>
#include <format>

using TF = PowerDiagram::TF;

struct Flower {
    using Pt = Vec<TF,2>;

    void init() {
        // points
        points.clear();
        for( PI n = 0; n < nb_samples; ++n ) {
            TF a = 2 * std::numbers::pi * n / nb_samples;
            TF r = 2 + sin( nb_lobes * a );
            points << Pt{
                r * cos( a ) + ( TF( rand() ) / RAND_MAX - 0.5 ) * dev,
                r * sin( a ) + ( TF( rand() ) / RAND_MAX - 0.5 ) * dev
            };
        }


        // dirac_points;
        dirac_points.clear();
        for( PI n = 0; n < 1000; ++n ) {
            dirac_points << Pt{
                TF( rand() ) / RAND_MAX + .5,
                TF( rand() ) / RAND_MAX
            };
        }
    }

    static VtkOutput::Pt to_vtk_point( const Pt &pt ) {
        return { pt[ 0 ], pt[ 1 ], 0.0 };
    }

    void plot_points( Str filename ) {
        VtkOutput vo;
        for( const auto &p : points )
            vo.add_point( to_vtk_point( p ) );
        vo.save( filename );
    }

    void plot_diracs( Str filename ) {
        VtkOutput vo;
        for( const auto &p : dirac_points )
            vo.add_point( to_vtk_point( p ) );
        vo.save( filename );
    }

    std::pair<Vec<TF>,Vec<TF>> make_density( Pt proj_dir ) {
        TF max_pos = std::numeric_limits<TF>::lowest();
        TF min_pos = std::numeric_limits<TF>::max();

        for( const Pt &p : targets ) {
            TF v = sp( p, proj_dir );
            max_pos = max( max_pos, v );
            min_pos = min( min_pos, v );
        }

        Vec<TF> x{ FromSizeAndFunctionOnIndex(), nb_bins, [&]( PI ind ) { return min_pos + ( max_pos - min_pos ) * ind / ( nb_bins - 1 ); } };
        Vec<TF> y{ FromSizeAndItemValue(), nb_bins, 0 };
        for( const Pt &p : points ) {
            TF v = sp( p, proj_dir );
            
            TF pos = ( v - min_pos ) * ( nb_bins - 1 ) / ( max_pos - min_pos );
            PI ind = min( PI( pos ), nb_bins - 1 );
            TF fra = pos - ind;

            y[ ind + 0 ] += ( 1 - fra );
            y[ ind + 1 ] += fra;
        }
        y.pop_back();

        return { x, y };
    }

    Vec<TF> delta_for_dir( Pt proj_dir ) {
        Vec<TF> projected_dirac_points( FromReservationSize(), dirac_points.size() );
        for( const Pt &p : dirac_points )
            projected_dirac_points << sp( p, proj_dir );

        std::pair<Vec<TF>,Vec<TF>> projected_density = make_density( proj_dir );

        InitialSolution solver = make_initial_solution( projected_dirac_points, projected_density.first, projected_density.second, /*mass_ratio =*/ 1.0 );

        return solver.barycenters() - projected_dirac_points;
    }

    void quant_step() {
        Vec<Vec<Pt,2>> M( FromSizeAndItemValue(), dirac_points.size(), Vec<Pt,2>{ Pt{ 0, 0 }, Pt{ 0, 0 } } );
        Vec<Pt> V( FromSizeAndItemValue(), dirac_points.size(), Pt{ 0, 0 } );
        for( PI num_dir = 0; num_dir < nb_dirs; ++num_dir ) {
            // TF a = num_dir * std::numbers::pi / nb_dirs;
            TF a = rand() * std::numbers::pi / RAND_MAX;
            Pt proj_dir{ cos( a ), sin( a ) };
            
            Vec<TF> de = delta_for_dir( proj_dir );
            for( PI num_dirac = 0; num_dirac < dirac_points.size(); ++num_dirac ) {
                M[ num_dirac ][ 0 ][ 0 ] += proj_dir[ 0 ] * proj_dir[ 0 ];
                M[ num_dirac ][ 0 ][ 1 ] += proj_dir[ 0 ] * proj_dir[ 1 ];
                M[ num_dirac ][ 1 ][ 0 ] += proj_dir[ 1 ] * proj_dir[ 0 ];
                M[ num_dirac ][ 1 ][ 1 ] += proj_dir[ 1 ] * proj_dir[ 1 ];
                V[ num_dirac ][ 0 ] += proj_dir[ 0 ] * de[ num_dirac ];
                V[ num_dirac ][ 1 ] += proj_dir[ 1 ] * de[ num_dirac ];
            }
        }

        for( PI num_dirac = 0; num_dirac < dirac_points.size(); ++num_dirac ) {
            using TM = Eigen::Matrix<TF,2,2>;
            using TV = Eigen::Matrix<TF,2,1>;
            TM eM{ 
                { M[ num_dirac ][ 0 ][ 0 ], M[ num_dirac ][ 0 ][ 1 ] },
                { M[ num_dirac ][ 1 ][ 0 ], M[ num_dirac ][ 1 ][ 1 ] },
            };
            TV eV = { V[ num_dirac ][ 0 ], V[ num_dirac ][ 1 ] };
            Eigen::FullPivLU<TM> lu( eM );
            Pt new_dir = lu.solve( eV );

            dirac_points[ num_dirac ] += new_dir;
        }
    }

    PI      nb_samples = 1000000;
    PI      nb_lobes   = 5;
    PI      nb_bins    = 200;
    PI      nb_dirs    = 2000;
    TF      dev        = 1e-2;

    Vec<Pt> projection_dirs;
    Vec<Pt> dirac_points;
    Vec<Pt> points;
};

TEST_CASE( "Flower ", "" ) {
    Flower f;
    f.init();

    f.plot_diracs( "diracs_0.vtk" );
    for( PI num_step = 1; num_step <= 300; ++num_step ) {
        P( num_step );
        f.quant_step();
        f.plot_diracs( std::format( "results/diracs_{}.vtk", num_step ) );
    }
    // auto d = f.make_density( { 1, 0 }, 300 );
    // matplotlibcpp::plot( as_std_vec( d.first ), as_std_vec( d.second ), {} );
    // matplotlibcpp::show();
}
