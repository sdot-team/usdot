#pragma once

#include <tl/support/operators/sp.h>
#include "ShapeRegistration.h"
#include "ConvGridSolver.h"
#include "partial1D/ConvGridSolverInput.h"
#include <fstream>

namespace usdot {
    
#define DTP template<class TF,int dim>
#define UTP ShapeRegistration<TF,dim>

DTP void UTP::load_triangles( Vec<Triangle> &triangles, Str filename ) {
    std::ifstream is( filename );
    std::string line, word;
    Triangle t;
    PI n = 0;
    while ( std::getline( is, line ) ) {
        std::istringstream iss( line );
        iss >> word;
        if ( ! is )
            break;

        if ( word == "vertex" ) {
            iss >> t.points[ n ][ 0 ] >> t.points[ n ][ 1 ] >> t.points[ n ][ 2 ];    
            if ( n == 2 ) {
                triangles << t;
                n = 0;
            } else
                n += 1;
        }
    }
}

DTP void UTP::load_points( Vec<Pt> &points, Str filename ) {
    std::ifstream is( filename );
    while ( true ) {
        char t;
        Pt p;
        is >> t;
        for( PI i = 0; i < dim; ++i )
            is >> p[ i ];    
        if ( ! is )
            break;
        points << p;
    }
}

DTP void UTP::load_shape( Str filename ) {
    TODO;
}

DTP void UTP::glot( Vec<TF> xs, auto &&...funcs ) {
    std::ofstream fs( "glot.py" );
    fs << "from matplotlib import pyplot\n";
    auto pf = [&]( auto &&func ) {
        Vec<FP64> ys;
        for( auto x : xs )
            ys << func( x );
        fs << "pyplot.plot( " << to_string( xs ) << ", " << to_string( ys ) << " )\n";
    };
    ( pf( funcs ), ... );
    fs << "pyplot.show()\n";
}

DTP void UTP::display_barycenters_vtk( VtkOutput &vo ) {
    for( PI n = 0; n < diracs.size(); ++n )  {
        Pt di = transformation.apply( diracs[ n ] );
        vo.add_line( Vec<VtkOutput::Pt>{ di, new_diracs[ n ] } );
    }
}

DTP void UTP::display_shape_vtk( VtkOutput &vo ) {
    for( const Pt &pt : shape_points )
        vo.add_point( VtkOutput::Pt( pt ) );
    for( const Triangle &t : shape_triangles )
        vo.add_triangle( { VtkOutput::Pt( t.points[ 0 ] ), VtkOutput::Pt( t.points[ 1 ] ), VtkOutput::Pt( t.points[ 2 ] ) } );
}

DTP void UTP::display_diracs_vtk( VtkOutput &vo ) {
    for( const Pt &pt : diracs )
        vo.add_point( VtkOutput::Pt( transformation.apply( pt ) ) );
}

DTP void UTP::display_vtk( Str filename, PI num ) {
    VtkOutput v1;
    display_diracs_vtk( v1 );
    v1.save( filename + "_diracs_" + std::to_string( num ) + ".vtk" );
    
    VtkOutput v2;
    display_shape_vtk( v2 );
    v2.save( filename + "_shape_" + std::to_string( num ) + ".vtk" );
    
    VtkOutput v3;
    display_barycenters_vtk( v3 );
    v3.save( filename + "_barycenters_" + std::to_string( num ) + ".vtk" );    
}

DTP Vec<typename UTP::Pt> UTP::projection_dirs() {
    // TF a = num_dir * std::numbers::pi / nb_dirs;
    Vec<Pt> res;

    if ( dim == 2 ) {
        TF o = rand() * M_PI / RAND_MAX;
        for( PI n = 0; n < nb_projection_dirs; ++n ) {
            TF a = o + n * M_PI / nb_projection_dirs;
            res << Pt{ cos( a ), sin( a ) };
        }
        return res;
    } else if ( dim == 3 ) {
        for( PI n = 0; n < nb_projection_dirs; ++n ) {
            TF u = rand() * 1.0 / RAND_MAX;
            TF v = rand() * 1.0 / RAND_MAX;
            TF theta = 2 * M_PI * u;
            TF phi = acos( 2 * v - 1 );
            TF x = sin( phi ) * cos( theta );
            TF y = sin( phi ) * sin( theta );
            TF z = cos( phi );

            res << Pt{ x, y, z };
        }
    } else {
        TODO;
    }

    return res;
}

DTP std::tuple<TF,TF,Vec<TF>> UTP::projected_density( Pt proj_dir ) {
    if ( shape_triangles.size() ) {
        ASSERT( shape_points.empty() );

        struct Event { TF pos, slope; };
        Vec<Event> events;
        for( const Triangle &tr : shape_triangles ) {
            TF x0 = sp( tr.points[ 0 ], proj_dir );
            TF x1 = sp( tr.points[ 1 ], proj_dir );
            TF x2 = sp( tr.points[ 2 ], proj_dir );
        }

        // return RcPtr<PiecewiseAffConstantDensity<TF>>::New( x, y );
    }

    // histogram
    TF max_pos = std::numeric_limits<TF>::lowest();
    TF min_pos = std::numeric_limits<TF>::max();
    for( const Pt &p : shape_points ) {
        TF v = sp( p, proj_dir );
        max_pos = max( max_pos, v );
        min_pos = min( min_pos, v );
    }

    Vec<TF> y{ FromSizeAndItemValue(), nb_bins, 0 };
    for( const Pt &p : shape_points ) {
        TF pos = ( sp( p, proj_dir ) - min_pos ) * nb_bins / ( max_pos - min_pos );
        PI ind = min( PI( pos ), nb_bins - 1 );
        y[ ind ] += 1;
    }

    return { min_pos, max_pos, y };
}

DTP Vec<TF> UTP::delta_for_dir( Pt proj_dir ) {
    ConvGridSolverInput<TF> solver_input;
    solver_input.global_mass_ratio = mass_ratio;
    solver_input.starting_filter_value = 0.5;
    solver_input.target_mass_error = 1e-3;
    solver_input.min_dirac_separation = 1e-4;

    // diracs
    solver_input.dirac_positions.resize( new_diracs.size() );
    for( PI n = 0; n < new_diracs.size(); ++n )
        solver_input.dirac_positions[ n ] = sp( new_diracs[ n ], proj_dir );

    // density
    std::tuple<TF,TF,Vec<TF>> de = projected_density( proj_dir );
    solver_input.density_values = std::get<2>( de );
    solver_input.beg_x_density = std::get<0>( de );
    solver_input.end_x_density = std::get<1>( de );

    // solver
    ConvGridSolver<TF> solver( std::move( solver_input ) );
    solver.sorted_dirac_weights *= 40;
    try {
        solver.solve();
    } catch ( std::runtime_error e ) {
        // P( solver.dirac_positions() );
        // P( mass_ratio );
        // P( de );
        solver.plot();
        P( e.what() );
        ASSERT( 0 );
    }
    // ASSERT( 0 );

    // delta
    return solver.cell_barycenters() - solver.dirac_positions();
}

DTP void UTP::compute_new_diracs( PI nb_iter ) {
    new_diracs.resize( diracs.size() );
    for( PI n = 0; n < diracs.size(); ++n )
        new_diracs[ n ] = transformation.apply( diracs[ n ] );

    for( PI num_iter_cnd = 0; num_iter_cnd < nb_iter; ++num_iter_cnd ) {
        P( num_iter_cnd );
        Vec<Pt> proj_dirs = projection_dirs();

        Vec<Pt,dim> M;
        for( PI r = 0; r < dim; ++r )
            for( PI c = 0; c < dim; ++c )
                M[ r ][ c ] = 0;

        Vec<Pt> V( FromSize(), diracs.size() );
        for( auto &v : V )
            for( PI r = 0; r < dim; ++r )
                v[ r ] = 0;

        for( Pt proj_dir : proj_dirs ) {
            // P( proj_dir );

            Vec<TF> de = delta_for_dir( proj_dir );

            for( PI r = 0; r < dim; ++r )
                for( PI c = 0; c < dim; ++c )
                    M[ r ][ c ] += proj_dir[ r ] * proj_dir[ c ];

            for( PI num_dirac = 0; num_dirac < diracs.size(); ++num_dirac )
                for( PI r = 0; r < dim; ++r )
                    V[ num_dirac ][ r ] += proj_dir[ r ] * de[ num_dirac ];
        }
            
        // make matrix        
        using TM = Eigen::Matrix<TF,dim,dim>;
        TM eM;
        for( PI r = 0; r < dim; ++r )
            for( PI c = 0; c < dim; ++c )
                eM.coeffRef( r, c ) = M[ r ][ c ];
        Eigen::FullPivLU<TM> lu( eM );

        // solve
        using TV = Eigen::Matrix<TF,dim,1>;
        for( PI num_dirac = 0; num_dirac < diracs.size(); ++num_dirac ) {
            TV eV;
            for( PI r = 0; r < dim; ++r )
                eV[ r ] = V[ num_dirac ][ r ];

            new_diracs[ num_dirac ] += 0.95 * Pt( lu.solve( eV ) );
        }
    }
}

#undef DTP
#undef UTP

// TF iteration_SVD() {
//     // centers
//     Pt center_old( FromItemValue(), 0 );
//     Pt center_new( FromItemValue(), 0 );
//     for( PI num_dirac = 0; num_dirac < diracs.size(); ++num_dirac ) {
//         Pt old_pos = dirac_trans.apply( diracs[ num_dirac ] );
//         Pt new_pos = new_diracs[ num_dirac ];
//         center_old += old_pos;
//         center_new += new_pos;
//     }
//     center_old /= diracs.size();
//     center_new /= diracs.size();

//     // covariance
//     TM cov = TM::Zero();
//     for( PI num_dirac = 0; num_dirac < diracs.size(); ++num_dirac ) {
//         Pt old_pos = dirac_trans.apply( diracs[ num_dirac ] );
//         Pt new_pos = new_diracs[ num_dirac ];
//         for( PI r = 0; r < dim; ++r )
//             for( PI c = 0; c < dim; ++c )
//                 cov.coeffRef( r, c ) += ( new_pos[ r ] - center_new[ r ] ) * ( old_pos[ c ] - center_old[ c ] );
//     }
//     Eigen::JacobiSVD<TM> svd( cov, Eigen::ComputeFullU | Eigen::ComputeFullV );
//     // std::cout << svd.singularValues() << std::endl;
//     // std::cout << svd.matrixU() << std::endl;
//     // std::cout << svd.matrixV() << std::endl;

//     TM orth = svd.matrixU() * svd.matrixV().transpose();
//     // std::cout << orth << std::endl;
//     // std::cout << "det " << orth.determinant() << std::endl;

//     TM rot = svd.matrixU() * svd.matrixV().transpose();
//     // PE( ( rot - TM::Identity() ).norm() );
//     double err = 0;
//     err += pow( 1 - rot( 0 ), 2 );
//     err += pow(     rot( 1 ), 2 );
//     err += pow(     rot( 2 ), 2 );
//     err += pow(     rot( 3 ), 2 );
//     err += pow( 1 - rot( 4 ), 2 );
//     err += pow(     rot( 5 ), 2 );
//     err += pow(     rot( 6 ), 2 );
//     err += pow(     rot( 7 ), 2 );
//     err += pow( 1 - rot( 8 ), 2 );
//     std::cout << std::sqrt( err ) << std::endl;

//     dirac_trans *= Tr::translation( - center_old ) * Tr::linear( rot ) * Tr::translation( center_new );

//     return std::sqrt( err );
// }

} // namespace usdot
