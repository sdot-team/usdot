#pragma once

#include "Density/PiecewiseConstantDensity.h"
#include <tl/support/operators/sp.h>
#include "ShapeRegistration.h"
#include "Solver.h"
#include <fstream>

namespace usdot {
    
#define DTP template<class TF,int dim>
#define UTP ShapeRegistration<TF,dim>

DTP void UTP::load_triangles( Vec<Triangle> &points, Str filename ) {
    TODO;
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

DTP void UTP::display_displacements_vtk( VtkOutput &vo ) {
    for( PI n = 0; n < diracs.size(); ++n )  {
        Pt di = transformation.apply( diracs[ n ] );
        vo.add_line( Vec<VtkOutput::Pt>{ di, di + displacements[ n ] } );
    }
}

DTP void UTP::display_shape_vtk( VtkOutput &vo ) {
    for( const Pt &pt : shape_points )
        vo.add_point( VtkOutput::Pt( pt ) );
}

DTP void UTP::display_diracs_vtk( VtkOutput &vo ) {
    for( const Pt &pt : diracs )
        vo.add_point( VtkOutput::Pt( transformation.apply( pt ) ) );
}

DTP Vec<typename UTP::Pt> UTP::get_proj_dirs() {
    // TF a = num_dir * std::numbers::pi / nb_dirs;
    Vec<Pt> res;

    if ( dim == 2 ) {
        TF o = rand() * M_PI / RAND_MAX;
        for( PI n = 0; n < nb_proj_dirs; ++n ) {
            TF a = o + n * M_PI / nb_proj_dirs;
            res << Pt{ cos( a ), sin( a ) };
        }
        return res;
    } else if ( dim == 3 ) {
        for( PI n = 0; n < nb_proj_dirs; ++n ) {
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

DTP std::pair<Vec<TF>,Vec<TF>> UTP::projected_density( Pt proj_dir ) {
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

    return { x, y };
}

DTP Vec<TF> UTP::delta_for_dir( Pt proj_dir ) {
    // projected_density
    std::pair<Vec<TF>,Vec<TF>> pde = projected_density( proj_dir );
    auto de = RcPtr<PiecewiseConstantDensity<TF>>::New( pde.first, pde.second );

    // projected_diracs
    Vec<TF> pdi( FromReservationSize(), diracs.size() );
    for( PI n = 0; n < diracs.size(); ++n )
        pdi << sp( transformation.apply( diracs[ n ] ), proj_dir );

    // glot( Vec<TF>::linspace( min( pde.first ), max( pde.first ), 10000, false ), [&]( TF x ) { return de->value( x ); } );
    // TODO;
    
    // power diagram
    auto pd = RcPtr<PowerDiagram<TF>>::New( pdi, Vec<TF>::fill( diracs.size(), 1 ) );

    // // auto t0 = std::chrono::high_resolution_clock::now();
    // // auto t1 = std::chrono::high_resolution_clock::now();
    // // PE( std::chrono::duration<double>{ t1 - t0 } );
    // solver
    Solver<TF> solver( pd, de, mass_ratio );
    solver.initialize_weights();
    solver.update_weights();

    // delta
    Vec<TF> bar = pd->barycenters( *de );
    for( PI n = 0; n < diracs.size(); ++n )
        pdi[ pd->sorted_seed_nums[ n ] ] = bar[ n ] - pdi[ pd->sorted_seed_nums[ n ] ];
    return pdi;
}

DTP void UTP::update_displacements() {
    Vec<Pt> proj_dirs = get_proj_dirs();

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
    displacements.resize( diracs.size() );
    for( PI num_dirac = 0; num_dirac < diracs.size(); ++num_dirac ) {
        TV eV;
        for( PI r = 0; r < dim; ++r )
            eV[ r ] = V[ num_dirac ][ r ];

        displacements[ num_dirac ] = Pt( lu.solve( eV ) );
    }
}

#undef DTP
#undef UTP


} // namespace usdot
