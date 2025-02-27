// #include <tl/support/operators/norm_2.h>
#include <partial1D/InitialSolution.h>
#include <tl/support/operators/abs.h>
#include <tl/support/operators/sp.h>
#include <matplotlibcpp.h>
#include <numbers>

/**
 */
struct Proj {
    using TF = PowerDiagram::TF;
    using Pt = Vec<TF,2>;

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

        for( const Pt &p : points ) {
            TF v = sp( p, proj_dir );
            max_pos = max( max_pos, v );
            min_pos = min( min_pos, v );
        }

        // Vec<TF> x{ FromSizeAndFunctionOnIndex(), nb_bins, [&]( PI ind ) { return min_pos + ( max_pos - min_pos ) * ind / ( nb_bins - 1 ); } };
        // Vec<TF> y{ FromSizeAndItemValue(), nb_bins, 0 };
        // for( const Pt &p : points ) {
        //     TF v = sp( p, proj_dir );
            
        //     TF pos = ( v - min_pos ) * ( nb_bins - 1 ) / ( max_pos - min_pos );
        //     PI ind = min( PI( pos ), nb_bins - 1 );
        //     TF fra = pos - ind;

        //     y[ ind + 0 ] += ( 1 - fra );
        //     y[ ind + 1 ] += fra;
        // }
        Vec<TF> x{ FromSizeAndFunctionOnIndex(), nb_bins + 1, [&]( PI ind ) { return min_pos + ( max_pos - min_pos ) * ind / nb_bins; } };
        Vec<TF> y{ FromSizeAndItemValue(), nb_bins, 0 };
        for( const Pt &p : points ) {
            TF v = sp( p, proj_dir );
            
            TF pos = ( v - min_pos ) * ( nb_bins - 1 ) / ( max_pos - min_pos );
            PI ind = min( PI( pos ), nb_bins - 1 );
            y[ ind ] += 1;
        }

        return { x, y };
    }

    Vec<Pt> delta_for_dir( Pt proj_dir ) {
        Vec<TF> projected_dirac_points( FromReservationSize(), dirac_points.size() );
        for( const Pt &p : dirac_points )
            projected_dirac_points << sp( p, proj_dir );

        std::pair<Vec<TF>,Vec<TF>> projected_density = make_density( proj_dir );

        InitialSolution solver( projected_dirac_points, projected_density.first, projected_density.second, mass_ratio );

        Vec<TF> projected_delta = solver.barycenters() - projected_dirac_points;

        Vec<Pt> delta( FromReservationSize(), projected_delta.size() );
        for( const TF &p : projected_delta )
            delta << p * proj_dir;
        return delta;
    }

    void quant_step() {
        // Vec<Pt> new_dir( FromSizeAndItemValue(), dirac_points.size(), Pt{ 0, 0 } );
        for( PI num_dir = 0; num_dir < nb_dirs; ++num_dir ) {
            TF a = num_dir * std::numbers::pi / nb_dirs;
            Pt proj_dir{ cos( a ), sin( a ) };
            
            auto d = delta_for_dir( proj_dir );
            for( PI n = 0; n < d.size(); ++n )
                dirac_points[ n ] += 0.85 * d[ n ] / nb_dirs;
        }
    }

    TF      mass_ratio = 0.25;
    PI      nb_samples = 1000000;
    PI      nb_bins    = 200;
    PI      nb_dirs    = 100;

    PI      num_dir    = 0;

    Vec<Pt> projection_dirs;
    Vec<Pt> dirac_points;-
    Vec<Pt> points;
};
