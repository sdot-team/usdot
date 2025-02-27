#include <tl/support/operators/norm_2.h>
#include <tl/support/operators/abs.h>
#include <tl/support/operators/sp.h>
#include <eigen3/Eigen/LU>
#include <partial1D/Solver.h>
#include <matplotlibcpp.h>
#include "catch_main.h"
#include "eigen3/Eigen/src/LU/FullPivLU.h"
#include "partial1D/InitialSolution.h"
#include <numbers>
#include <limits>
#include <format>

using TF = PowerDiagram::TF;
using Pt = Vec<TF,2>;

struct Discs {
    static VtkOutput::Pt to_vtk_point( const Pt &pt ) {
        return { pt[ 0 ], pt[ 1 ], 0.0 };
    }

    void plot_diracs( Str filename ) {
        VtkOutput vo;
        for( const auto &p : dirac_points )
            vo.add_point( to_vtk_point( p ) );
        vo.save( filename );
    }

    std::pair<Vec<TF>,Vec<TF>> make_density( Pt proj_dir ) {
        using std::sqrt;
        using std::pow;

        TF min_pos = std::numeric_limits<TF>::max();
        TF max_pos = std::numeric_limits<TF>::lowest();

        for( const std::pair<Pt,TF> &disc : discs ) {
            min_pos = min( min_pos, sp( disc.first, proj_dir ) - disc.second );
            max_pos = max( max_pos, sp( disc.first, proj_dir ) + disc.second );
        }

        Vec<TF> x{ FromSizeAndFunctionOnIndex(), nb_bins, [&]( PI ind ) { return min_pos + ( max_pos - min_pos ) * ind / ( nb_bins - 1 ); } };
        Vec<TF> y{ FromSizeAndItemValue(), nb_bins, 0 };
        for( PI i = 0; i < nb_bins; ++i ) {
            for( const std::pair<Pt,TF> &disc : discs ) {
                TF p = sp( disc.first, proj_dir );
                TF r = pow( disc.second, 2 );
                TF d = pow( x[ i ] - p, 2 );
                if ( d < r )
                    y[ i ] += sqrt( r - d );
            }
        }

        y.pop_back();

        return { x, y };
    }

    Vec<TF> delta_for_dir( Pt proj_dir ) {
        Vec<TF> projected_dirac_points( FromReservationSize(), dirac_points.size() );
        for( const Pt &p : dirac_points )
            projected_dirac_points << sp( p, proj_dir );

        std::pair<Vec<TF>,Vec<TF>> projected_density = make_density( proj_dir );

        InitialSolution solver = make_initial_solution( projected_dirac_points, projected_density.first, projected_density.second, /*mass_ratio =*/ 5.0 / 5.0 );

        return solver.barycenters() - projected_dirac_points;
    }

    void quant_step() {
        Vec<Vec<Pt,2>> M( FromSizeAndItemValue(), dirac_points.size(), Vec<Pt,2>{ Pt{ 0, 0 }, Pt{ 0, 0 } } );
        Vec<Pt> V( FromSizeAndItemValue(), dirac_points.size(), Pt{ 0, 0 } );
        TF o = rand() * std::numbers::pi / RAND_MAX;
        for( PI num_dir = 0; num_dir < nb_dirs; ++num_dir ) {
            TF a = o + num_dir * std::numbers::pi / nb_dirs;
            // TF a = rand() * std::numbers::pi / RAND_MAX;
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

    PI      nb_bins    = 1000;
    PI      nb_dirs    = 3;

    Vec<Pt> dirac_points;
    Vec<std::pair<Pt,TF>> discs;
};

TEST_CASE( "Discs ", "" ) {
    Discs f;

    f.discs << std::pair<Pt,TF>{ { 0, 0 }, 1 };
    f.discs << std::pair<Pt,TF>{ { 1.5, 0.0 }, 0.5 };

    // for( PI n = 0; n < 2000; ++n ) {
    //     Pt p{
    //         ( 2.0 * rand() / RAND_MAX - 1.0 ),
    //         ( 2.0 * rand() / RAND_MAX - 1.0 ),
    //     };
    //     if ( norm_2_p2( p ) <= 1 )
    //         f.dirac_points << p;
    // }
    for( TF x = -1; x <= 1; x += 3e-2 ) {
        for( TF y = -1; y <= 1; y += 3e-2 ) {
            Pt p{ x, y };
            if ( norm_2_p2( p ) <= 1 )
                f.dirac_points << p;
        }
    }

    f.plot_diracs( std::format( "results/diracs_{}_{}.vtk", f.nb_dirs, 0 ) );

    for( PI num_step = 1; num_step <= 400; ++num_step ) {
        P( num_step );
        f.quant_step();
        f.plot_diracs( std::format( "results/diracs_{}_{}.vtk", f.nb_dirs, num_step ) );
    }

    // auto d = f.make_density( { 1, 0 }, 300 );
    // matplotlibcpp::plot( as_std_vec( d.first ), as_std_vec( d.second ), {} );
    // matplotlibcpp::show();
}
