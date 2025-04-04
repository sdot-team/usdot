#pragma once

#include <tl/support/ASSERT.h>
#include <tl/support/P.h>
#include <numeric>
#include <fstream>

#include "WeightInitializer.h"
#include "System.h"

namespace usdot {
    
#define DTP template<class TF,class Density>
#define UTP System<TF,Density>

DTP UTP::System() {
    // using namespace std;
    global_mass_ratio = 1;
    total_dirac_mass = -1;
    density = nullptr;
}

DTP void UTP::initialize_weights() {
    _update_system( true );
    
    WeightInitializer<TF,Density> wi( *this );
    wi.run();
}

DTP PI UTP::nb_diracs() const {
    return sorted_dirac_positions.size();
}

DTP TF UTP::density_value( TF pos ) const {
    return density->value( pos );
}

DTP typename UTP::VF UTP::dirac_barycenters() const {
    VF res( FromSize(), nb_diracs() );
    for( TI i = 0; i < nb_diracs(); ++i )
        res[ sorted_dirac_nums[ i ] ] = sorted_dirac_positions[ i ];
    return res;
}

DTP void UTP::plot( Str filename ) const {
    std::ofstream fs( filename );

    // TF de = ( end_x_density - beg_x_density ) / 3;
    // TV xs = Vec<TF>::linspace( beg_x_density - de, end_x_density + de, 10000 );
    // TV ys;
    // for( TF x : xs )
    //     ys << density_value( x );

    // // for( PI i = 0; i < original_density_values.size(); ++i ) {
    // //     const TF x = beg_x_density + ( end_x_density - beg_x_density ) * i / ( original_density_values.size() - 1 );
    // //     ys.push_back( original_density_values[ i ] );
    // //     xs.push_back( x );
    // // }
    // // ys.push_back( original_density_values.back() );
    // // xs.push_back( end_x_density );



    // TV bx = cell_barycenters();
    // TV by( FromSizeAndItemValue(), bx.size(), -0.1 );

    // TV dx = dirac_positions();
    // TV dy( FromSizeAndItemValue(), dx.size(), -0.2 );

    // fs << "from matplotlib import pyplot\n";
    // fs << "pyplot.plot( " << to_string( xs ) << ", " << to_string( ys ) << " )\n";
    // fs << "pyplot.plot( " << to_string( bx ) << ", " << to_string( by ) << ", '+' )\n";
    // fs << "pyplot.plot( " << to_string( dx ) << ", " << to_string( dy ) << ", '+' )\n";
    // fs << "pyplot.show()\n";
}

template<class TF>
void plot_bnds_evolution( const Vec<Vec<Vec<TF,2>>> &bnds ) {
    auto ys = Vec<TF>::linspace( 0, 1, bnds.size() );
    std::ofstream fs( "glot.py" );
    fs << "from matplotlib import pyplot\n";
    for( PI i = 0; i < bnds[ 0 ].size(); ++i ) {
        for( PI j = 0; j < bnds[ 0 ][ 0 ].size(); ++j ) {
            Vec<TF> xs( FromSizeAndFunctionOnIndex(), bnds.size(), [&]( PI n ) { return bnds[ n ][ i ][ j ]; } );
            fs << "pyplot.plot( " << to_string( xs ) << ", " << to_string( ys ) << " )\n";
        }
    }
    fs << "pyplot.show()\n";
}

DTP void UTP::for_each_cell( auto &&func ) const {
    using namespace std;

    if ( nb_diracs() == 0 )
        return;

    TF i0 = numeric_limits<TF>::lowest();
    TF d0 = sorted_dirac_positions[ 0 ];
    TF w0 = sorted_dirac_weights[ 0 ];

    for( PI n = 1; n < nb_diracs(); ++n ) {
        const TF d1 = sorted_dirac_positions[ n ];
        const TF w1 = sorted_dirac_weights[ n ];
        const TF i1 = ( d0 + d1 + ( w0 - w1 ) / ( d1 - d0 ) ) / 2;
        const TF r0 = sqrt( max( w0, 0 ) );
        const TF b0 = d0 - r0;
        const TF e0 = d0 + r0;

        func( n - 1, max( i0, b0 ), min( i1, e0 ) );

        i0 = i1;
        d0 = d1;
        w0 = w1;
    }

    const TF r0 = sqrt( max( w0, 0 ) );
    const TF b0 = d0 - r0;
    const TF e0 = d0 + r0;

    func( nb_diracs() - 1, max( i0, b0 ), e0 );
}

DTP typename UTP::VF UTP::cell_masses() const {
    _update_system();
    VF res( FromSize(), nb_diracs() );
    for_each_cell( [&]( TI n, TF b, TF e ) {
        res[ sorted_dirac_nums[ n ] ] = density->integral( b, e );
    } );
    return res;
}

DTP TF UTP::l2_mass_error() const {
    using namespace std;
    _update_system();
    TF res = 0;
    for_each_cell( [&]( TI n, TF b, TF e ) {
        res += pow( sorted_dirac_masses[ n ] - density->integral( b, e ), 2 );
    } );
    return sqrt( res );
}

DTP void UTP::set_dirac_positions( const VF &dirac_positions ) {
    using namespace std;

    // sorted_dirac_nums
    sorted_dirac_nums.resize( dirac_positions.size() );
    iota( sorted_dirac_nums.begin(), sorted_dirac_nums.end(), PI( 0 ) );
    sort( sorted_dirac_nums.begin(), sorted_dirac_nums.end(), [&]( PI a, PI b ) {
        return dirac_positions[ a ] < dirac_positions[ b ];
    } );

    // sorted_dirac_positions
    // const TF sep = input.min_dirac_separation * mul_dirac;
    sorted_dirac_positions.resize( dirac_positions.size() );
    for( PI i = 0; i < dirac_positions.size(); ++i ) {
        sorted_dirac_positions[ i ] = dirac_positions[ sorted_dirac_nums[ i ] ];
        // if ( i && sorted_dirac_positions[ i ] < sorted_dirac_positions[ i - 1 ] + sep ) 
        //     sorted_dirac_positions[ i ] = sorted_dirac_positions[ i - 1 ] + sep;
    }
}

DTP void UTP::set_global_mass_ratio( TF mass_ratio ) {
    global_mass_ratio = mass_ratio;
}

DTP void UTP::set_density( const Density *density ) {
    this->density = density;
}

DTP void UTP::_update_system( bool need_weights ) const {
    if ( sorted_dirac_masses.empty() ) {
        total_dirac_mass = global_mass_ratio * density->mass();
        sorted_dirac_masses = Vec<TF>::fill( nb_diracs(), total_dirac_mass / nb_diracs() );
    }

    if ( sorted_dirac_weights.empty() ) {
        const TF b = sorted_dirac_positions.front();
        const TF e = sorted_dirac_positions.back();
        const TF w = pow( global_mass_ratio * ( e - b ) / nb_diracs(), 2 );
        sorted_dirac_weights = Vec<TF>::fill( nb_diracs(), w );
    }
}

#undef DTP
#undef UTP

} // namespace usdot
