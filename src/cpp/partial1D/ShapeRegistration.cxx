#pragma once

#include "Density/PiecewiseConstantDensity.h"
#include <string>
#include <tl/support/operators/sp.h>
#include "ShapeRegistration.h"
#include "Solver.h"
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

DTP RcPtr<Density<TF>> UTP::projected_density( Pt proj_dir ) {
    TF max_pos = std::numeric_limits<TF>::lowest();
    TF min_pos = std::numeric_limits<TF>::max();

    for( const Pt &p : shape_points ) {
        TF v = sp( p, proj_dir );
        max_pos = max( max_pos, v );
        min_pos = min( min_pos, v );
    }

    Vec<TF> x = Vec<TF>::linspace( min_pos, max_pos, nb_bins + 1 );
    Vec<TF> y{ FromSizeAndItemValue(), nb_bins, 0 };
    for( const Pt &p : shape_points ) {
        TF pos = ( sp( p, proj_dir ) - min_pos ) * nb_bins / ( max_pos - min_pos );
        PI ind = min( PI( pos ), nb_bins - 1 );
        y[ ind ] += 1;
    }

    return RcPtr<PiecewiseConstantDensity<TF>>::New( x, y );
}

DTP Vec<TF> UTP::delta_for_dir( Pt proj_dir ) {
    // projected_density
    RcPtr<Density<TF>> de = projected_density( proj_dir );

    // projected_diracs
    Vec<TF> di( FromReservationSize(), new_diracs.size() );
    for( PI n = 0; n < diracs.size(); ++n )
        di << sp( new_diracs[ n ], proj_dir );

    // glot( Vec<TF>::linspace( min( pde.first ), max( pde.first ), 10000, false ), [&]( TF x ) { return de->value( x ); } );
    // TODO;
    
    // power diagram
    auto pd = RcPtr<PowerDiagram<TF>>::New( di, Vec<TF>::fill( diracs.size(), 1 ) );

    // // auto t0 = std::chrono::high_resolution_clock::now();
    // // auto t1 = std::chrono::high_resolution_clock::now();
    // // PE( std::chrono::duration<double>{ t1 - t0 } );

    // solver
    Solver<TF> solver( pd, de, mass_ratio );
    solver.initialize_weights();
    // solver.update_weights();

    // delta
    return pd->barycenters( *de, false ) - di;
}

DTP void UTP::compute_new_diracs( PI nb_iter ) {
    new_diracs.resize( diracs.size() );
    for( PI n = 0; n < diracs.size(); ++n )
        new_diracs[ n ] = transformation.apply( diracs[ n ] );

    for( PI num_iter = 0; num_iter < nb_iter; ++num_iter ) {
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
            for( PI r = 0; r < dim; ++r )
                for( PI c = 0; c < dim; ++c )
                    M[ r ][ c ] += proj_dir[ r ] * proj_dir[ c ];

            Vec<TF> de = delta_for_dir( proj_dir );
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

            new_diracs[ num_dirac ] += 0.5 * Pt( lu.solve( eV ) );
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
