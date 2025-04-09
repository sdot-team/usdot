#pragma once


#include "WeightInitializer.h"
#include "utility/linspace.h"
#include "WeightUpdater.h"
#include "System.h"
#include <limits>
#include <numeric>
#include <fstream>
#include <stdexcept>

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

    nb_iterations_init = wi.nb_iterations;
}

DTP void UTP::update_weights() {
    _update_system( true );

    WeightUpdater<TF,Density> wi( *this );
    wi.run();

    nb_iterations_update = wi.nb_iterations;
}

DTP void UTP::solve() {
    initialize_weights();
    plot();

    update_weights();

    if ( verbosity >= 2 && stream )
        *stream << "nb iteration init: " << nb_iterations_init << " update: " << nb_iterations_update << "\n";
}

DTP PI UTP::nb_original_diracs() const {
    return sorted_dirac_num_values.size();
}

DTP PI UTP::nb_sorted_diracs() const {
    return sorted_dirac_positions.size();
}

DTP TF UTP::density_value( TF pos ) const {
    return density->value( pos );
}

DTP typename UTP::VF UTP::dirac_positions() const {
    VF res( nb_original_diracs() );
    for( PI i = 0; i < nb_sorted_diracs(); ++i )
        for( PI o = sorted_dirac_num_offsets[ i + 0 ]; o < sorted_dirac_num_offsets[ i + 1 ]; ++o )
            res[ sorted_dirac_num_values[ o ] ] = sorted_dirac_positions[ i ];
    return res;
}

DTP typename UTP::VF UTP::dirac_weights() const {
    VF res( nb_original_diracs() );
    for( PI i = 0; i < nb_sorted_diracs(); ++i )
        for( PI o = sorted_dirac_num_offsets[ i + 0 ]; o < sorted_dirac_num_offsets[ i + 1 ]; ++o )
            res[ sorted_dirac_num_values[ o ] ] = sorted_dirac_weights[ i ];
    return res;
}

DTP TF UTP::x_tol() const {
    return target_max_mass_error * global_mass_ratio * density->width() / nb_sorted_diracs();
}

DTP void UTP::plot( Str filename ) const {
    std::ofstream fs( filename );

    fs << "from matplotlib import pyplot\n";
    density->plot( fs );

    TF y = 0;
    for( auto c : cell_boundaries() ) {
        const TF x0 = c[ 0 ];
        const TF x1 = c[ 1 ];
        const TF y0 = ( y++ ) / nb_sorted_diracs();
        const TF y1 = y0 + 1;
        fs << "pyplot.plot( [ " << x0 << ", " << x1 << ", " << x1 << ", " << x0 << ", " << x0 << " ], ";
        fs << "[ " << y0 << ", " << y0 << ", " << y1 << ", " << y1 << ", " << y0 << " ] )\n";
    }

    fs << "pyplot.show()\n";
}

DTP void UTP::plot_bnds_evolution( const std::vector<VB> &bnds ) {
    auto ys = linspace<TF>( 0, 1, bnds.size() );

    std::ofstream fs( "glot.py" );
    fs << "from matplotlib import pyplot\n";
    for( PI i = 0; i < bnds[ 0 ].size(); ++i ) {
        for( PI j = 0; j < bnds[ 0 ][ 0 ].size(); ++j ) {
            VF xs( bnds.size() );
            for( PI n = 0; n < bnds.size(); ++n )
                xs[ n ] = bnds[ n ][ i ][ j ];
            fs << "pyplot.plot( " << to_string( xs ) << ", " << to_string( ys ) << " )\n";
        }
    }
    fs << "pyplot.show()\n";
}

DTP T_T void UTP::for_each_cell( const T &func ) const {
    using namespace std;

    if ( nb_sorted_diracs() == 0 )
        return;

    TF i0 = numeric_limits<TF>::lowest();
    TF d0 = sorted_dirac_positions[ 0 ];
    TF w0 = sorted_dirac_weights[ 0 ];

    for( PI n = 1; n < nb_sorted_diracs(); ++n ) {
        const TF d1 = sorted_dirac_positions[ n ];
        const TF w1 = sorted_dirac_weights[ n ];
        const TF i1 = ( d0 + d1 + ( w0 - w1 ) / ( d1 - d0 ) ) / 2;
        const TF r0 = sqrt( max( w0, TF( 0 ) ) );
        const TF b0 = d0 - r0;
        const TF e0 = d0 + r0;

        func( n - 1, max( i0, b0 ), min( i1, e0 ) );

        i0 = i1;
        d0 = d1;
        w0 = w1;
    }

    const TF r0 = sqrt( max( w0, TF( 0 ) ) );
    const TF b0 = d0 - r0;
    const TF e0 = d0 + r0;

    func( nb_sorted_diracs() - 1, max( i0, b0 ), e0 );
}

DTP typename UTP::VF UTP::cell_barycenters() const {
    _update_system();

    VF res( nb_original_diracs() );
    for_each_cell( [&]( PI n, TF b, TF e ) {
        const TF it = density->integral( b, e );
        for( PI o = sorted_dirac_num_offsets[ n + 0 ]; o < sorted_dirac_num_offsets[ n + 1 ]; ++o )
            res[ sorted_dirac_num_values[ o ] ] = it ? density->x_integral( b, e ) / it : 0;
    } );
    return res;
}

DTP typename UTP::VB UTP::cell_boundaries() const {
    using namespace std;

    _update_system();

    VB res( nb_original_diracs() );
    TF mi = density->min_x();
    TF ma = density->max_x();
    for_each_cell( [&]( PI n, TF b, TF e ) {
        b = max( mi, min( ma, b ) );
        e = max( mi, min( ma, e ) );
        for( PI o = sorted_dirac_num_offsets[ n + 0 ]; o < sorted_dirac_num_offsets[ n + 1 ]; ++o )
            res[ sorted_dirac_num_values[ o ] ] = { b, e }; 
    } );
    return res;
}

DTP typename UTP::VF UTP::cell_masses() const {
    _update_system();
    VF res( nb_original_diracs() );
    for_each_cell( [&]( PI n, TF b, TF e ) {
        PI bo = sorted_dirac_num_offsets[ n + 0 ];
        PI eo = sorted_dirac_num_offsets[ n + 1 ];
        for( PI o = bo; o < eo; ++o )
            res[ sorted_dirac_num_values[ o ] ] = density->integral( b, e ) / ( eo - bo );
    } );
    return res;
}

DTP TF UTP::l2_mass_error() const {
    using namespace std;
    _update_system();
    TF res = 0;
    for_each_cell( [&]( PI n, TF b, TF e ) {
        res += pow( sorted_dirac_masses[ n ] - density->integral( b, e ), 2 );
    } );
    return sqrt( res );
}

DTP void UTP::set_dirac_positions( const VF &dirac_positions, TF min_dirac_separation ) {
    using namespace std;

    // sorted_dirac_nums
    VI sorted_dirac_nums( dirac_positions.size() );
    iota( sorted_dirac_nums.begin(), sorted_dirac_nums.end(), PI( 0 ) );
    sort( sorted_dirac_nums.begin(), sorted_dirac_nums.end(), [&]( PI a, PI b ) {
        return dirac_positions[ a ] < dirac_positions[ b ];
    } );

    // sorted_dirac_positions
    // const TF eps = numeric_limits<TF>::epsilon() * 100 * density->width();
    sorted_dirac_num_offsets.reserve( dirac_positions.size() );
    sorted_dirac_num_values.reserve( dirac_positions.size() );
    sorted_dirac_positions.reserve( dirac_positions.size() );
    for( PI i = 0; i < sorted_dirac_nums.size(); ++i ) {
        const PI s = sorted_dirac_num_values.size();
        const PI n = sorted_dirac_nums[ i ];
        const TF p = dirac_positions[ n ];
        if ( s && sorted_dirac_positions.back() + min_dirac_separation >= p ) {
            sorted_dirac_num_values.push_back( n );
            continue;
        }

        sorted_dirac_num_offsets.push_back( s );
        sorted_dirac_num_values.push_back( n );
        sorted_dirac_positions.push_back( p );
    }
    sorted_dirac_num_offsets.push_back( sorted_dirac_num_values.size() );
}

DTP void UTP::set_relative_dirac_masses( const VF &relative_mass_ratios ) {
    this->relative_mass_ratios = relative_mass_ratios;
}

DTP void UTP::set_global_mass_ratio( const TF &global_mass_ratio ) {
    this->global_mass_ratio = global_mass_ratio;
}

DTP void UTP::set_density( const Density *density ) {
    this->density = density;
}

DTP void UTP::_update_system( bool need_weights ) const {
    if ( ! density )
        throw std::runtime_error( "density is not defined" );
    
    if ( sorted_dirac_masses.empty() ) {
        total_dirac_mass = global_mass_ratio * density->mass();
        const TF bm = total_dirac_mass / nb_original_diracs();
        sorted_dirac_masses.resize( nb_sorted_diracs() );
        for( PI i = 0; i < nb_sorted_diracs(); ++i )
            sorted_dirac_masses[ i ] = bm * ( sorted_dirac_num_offsets[ i + 1 ] - sorted_dirac_num_offsets[ i + 0 ] );
    }

    if ( sorted_dirac_weights.empty() ) {
        const TF b = sorted_dirac_positions.front();
        const TF e = sorted_dirac_positions.back();
        const TF w = pow( global_mass_ratio * ( e - b ) / nb_sorted_diracs(), 2 );
        sorted_dirac_weights = fill<TF>( nb_sorted_diracs(), w );
    }
}

#undef DTP
#undef UTP

} // namespace usdot
