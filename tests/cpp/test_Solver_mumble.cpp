#include <tl/support/operators/mean.h>
#include <tl/support/operators/all.h>
#include <tl/support/operators/sp.h>
#include <partial1D/InitialSolution.h>
#include <partial1D/Solver.h>
#include "catch_main.h"
#include <fstream>

// #include <matplotlibcpp.h>

static constexpr int dim = 3;
using TF = PowerDiagram::TF;
using Pt = Vec<TF,dim>;

T_T std::vector<T> as_std_vec( const Vec<T> &t ) {
    return { t.begin(), t.end() };
}

TF my_norm_2( const Vec<TF> &a, const Vec<TF> &b ) {
    TF res = ( a[ 0 ] - b[ 0 ] ) * ( a[ 0 ] - b[ 0 ] );
    for( std::size_t i = 1; i < a.size(); ++i )
        res += ( a[ i ] - b[ i ] ) * ( a[ i ] - b[ i ] );
    using std::sqrt;
    return sqrt( res );
}

void load_points( Vec<Pt> &res, std::string filename ) {
    std::ifstream is( filename );
    while ( true ) {
        TF x, y, z;
        char t;
        is >> t >> x >> y >> z;
        if ( ! is )
            break;
        res << Pt{ x, y, z };
    }
}

std::pair<Vec<TF>,Vec<TF>> make_density( Pt proj_dir, const Vec<Pt> &target, PI nb_bins = ((((16+1)*2+1)*2+1)*2+1) ) {
    TF max_pos = std::numeric_limits<TF>::lowest();
    TF min_pos = std::numeric_limits<TF>::max();

    for( const Pt &p : target ) {
        TF v = sp( p, proj_dir );
        max_pos = max( max_pos, v );
        min_pos = min( min_pos, v );
    }
    TF delta = max_pos - min_pos;
    min_pos -= delta * 1e-3;
    max_pos += delta * 1e-3;

    Vec<TF> x{ FromSizeAndFunctionOnIndex(), nb_bins, [&]( PI ind ) { return min_pos + ( max_pos - min_pos ) * ind / ( nb_bins - 1 ); } };
    Vec<TF> y{ FromSizeAndItemValue(), nb_bins, 0 };
    for( const Pt &p : target ) {
        TF v = sp( p, proj_dir );
        
        TF pos = ( v - min_pos ) * ( nb_bins - 1 ) / ( max_pos - min_pos );
        PI ind = min( PI( pos ), nb_bins - 2 );
        TF fra = pos - ind;

        if ( ind + 0 < y.size() ) y[ ind + 0 ] += ( 1 - fra );
        if ( ind + 1 < y.size() ) y[ ind + 1 ] += fra;
    }
    y.pop_back();

    for( TF &v : y )
        if ( v == 0 )
            v = 1e-3;

    return { x, y };
}


double solve_for_dir( Pt proj_dir, Vec<Pt> diracs, const Vec<Pt> &target, TF ratio, Vec<PI> &nb_systems, Vec<PI> &nb_errors, bool use_init = true ) {
    // projected density
    std::pair<Vec<TF>,Vec<TF>> de = make_density( proj_dir, target );

    // projected_diracs
    Vec<TF> projected_diracs( FromReservationSize(), diracs.size() );
    for( PI n = 0; n < diracs.size(); ++n )
        projected_diracs << sp( diracs[ n ], proj_dir );

    auto t0 = std::chrono::high_resolution_clock::now();

    Solver is( projected_diracs, de.first, de.second, 1e-3, ratio );
    is.solve();

    auto t1 = std::chrono::high_resolution_clock::now();
    
    return std::chrono::duration<double>{ t1 - t0 }.count();
}

TEST_CASE( "Solver 80000", "" ) {
    // load points
    Vec<Pt> target;
    load_points( target, "../spot/Datasets/Pointsets/3D/mumble_sitting_100000.pts" );

    Vec<Pt> ciracs;
    load_points( ciracs, "../spot/Datasets/Pointsets/3D/mumble_sitting_80000.pts" );

    Vec<Pt> diracs;
    for( PI n = 0; n < ciracs.size(); ++n ) {
        Pt dir{ 0.42546746810417513, 0.35717031177892833, -0.8315087503861676 };
        Pt ori{ -26.23, 1.35, 50.83 };
        if ( sp( dir, ciracs[ n ] ) < sp( dir, ori ) )
            diracs << ciracs[ n ];
    }

    //
    Vec<double> timings;
    Vec<PI> nb_systems;
    Vec<PI> nb_errors;
    for( PI i = 0; i < 1000; ++i ) {
        TF u = rand() * 1.0 / RAND_MAX;
        TF v = rand() * 1.0 / RAND_MAX;
        TF theta = 2 * M_PI * u;
        TF phi = acos( 2 * v - 1 );
        TF x = sin( phi ) * cos( theta );
        TF y = sin( phi ) * sin( theta );
        TF z = cos( phi );

        Pt proj_dir{ x, y, z };

        timings << solve_for_dir( proj_dir, diracs, target, TF( diracs.size() ) / ciracs.size(), nb_systems, nb_errors );
    }
    P( mean( nb_systems ) );
    P( mean( nb_errors ) );
    P( mean( timings ) );
}
