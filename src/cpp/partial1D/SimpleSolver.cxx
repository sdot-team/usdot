#pragma once

// #include <tl/support/operators/norm_inf.h>
// #include <tl/support/operators/sum.h>
// #include <tl/support/operators/pow.h>
// #include <tl/support/operators/abs.h>
// #include <tl/support/string/va_string.h>

#include <stdexcept>
#include <tl/support/operators/norm_2.h>
#include <tl/support/operators/max.h>
#include <tl/support/operators/sum.h>
#include <tl/support/ASSERT.h>
#include <tl/support/P.h>
#include <limits>

#include "SimpleSolver.h"
#include "ThreadPool.h"

namespace usdot {

    
#define DTP template<class TF>
#define UTP SimpleSolver<TF>

DTP UTP::SimpleSolver( const SimpleSolverInput<TF> &input ) {
    using namespace std;
    
    if ( input.dirac_positions.empty() || input.density_values.empty() )
        return;

    // sorted_dirac_nums
    sorted_dirac_nums = { FromSizeAndFunctionOnIndex(), input.dirac_positions.size(), []( auto i ) { return i; } };
    sort( sorted_dirac_nums.begin(), sorted_dirac_nums.end(), [&]( PI a, PI b ) {
        return input.dirac_positions[ a ] < input.dirac_positions[ b ];
    } );

    //
    const TF end_x_density = input.end_x_density ? *input.end_x_density : input.density_values.size();
    const TF beg_x_density = input.beg_x_density ? *input.beg_x_density : 0;
    const TF step_size = ( end_x_density - beg_x_density ) / input.density_values.size();

    const TF flt_left = input.dirac_positions[ sorted_dirac_nums.front() ] - beg_x_density;
    const SI int_left = floor( flt_left / step_size );
    
    const TF flt_right = input.dirac_positions[ sorted_dirac_nums.back() ] - beg_x_density;
    const SI int_right = floor( flt_right / step_size + 1 );

    add_right = max( input.density_values.size(), int_right ) - input.density_values.size();
    add_left = max( SI( 0 ), - int_left );
    
    end_x = end_x_density + add_right * step_size;
    beg_x = beg_x_density - add_left * step_size;

    // sorted_dirac_positions, sorted_dirac_masses
    const PI total_size = add_left + input.density_values.size() + add_right;
    const PI nb_diracs = input.dirac_positions.size();

    const TF sep = input.relative_dirac_separation * total_size / nb_diracs;
    const TF dirac_mul = total_size / ( end_x - beg_x );
    sorted_dirac_positions.resize( nb_diracs );
    sorted_dirac_masses.resize( nb_diracs );
    sum_of_dirac_masses = 0;
    for( PI i = 0; i < nb_diracs; ++i ) {
        sorted_dirac_positions[ i ] = ( input.dirac_positions[ sorted_dirac_nums[ i ] ] - beg_x ) * dirac_mul;
        if ( i && sorted_dirac_positions[ i ] < sorted_dirac_positions[ i - 1 ] + sep ) 
            sorted_dirac_positions[ i ] = sorted_dirac_positions[ i - 1 ] + sep;
        
        const TF dirac_mass = input.dirac_masses.size() ? input.dirac_masses[ sorted_dirac_nums[ i ] ] : 1;
        sorted_dirac_masses[ i ] = dirac_mass;
        sum_of_dirac_masses += dirac_mass;
    }

    // sorted_dirac_weights
    sorted_dirac_weights = Vec<TF>::fill( input.dirac_positions.size(), 1 );

    // original_density_values
    max_of_original_density_values = max( input.density_values );
    original_density_values = input.density_values;
    global_mass_ratio = input.global_mass_ratio;

    //
    set_density_contrast( input.starting_contrast_ratio );
}

DTP PI UTP::nb_diracs() const {
    return sorted_dirac_positions.size();
}

DTP void UTP::for_each_cell_mt( auto &&func ) const {
    using namespace std;

    PI nb_jobs = 1 * thread_pool.nb_threads();
    thread_pool.execute( nb_jobs, [&]( PI num_job, PI num_thread ) {
        const PI beg = ( num_job + 0 ) * nb_diracs() / nb_jobs;
        const PI end = ( num_job + 1 ) * nb_diracs() / nb_jobs;

        //
        TF i0 = numeric_limits<TF>::lowest(); // previous interface
        TF o0 = numeric_limits<TF>::max(); // previous distance
        TF d0 = sorted_dirac_positions[ beg ];
        TF w0 = sorted_dirac_weights[ beg ];
        if ( beg ) {
            const TF dp = sorted_dirac_positions[ beg - 1 ];
            const TF wp = sorted_dirac_weights[ beg - 1 ];
            i0 = ( dp + d0 + ( wp - w0 ) / ( d0 - dp ) ) / 2;
            o0 = d0 - dp;    
        }

        for( PI n = beg + 1; n < end; ++n ) {
            const TF d1 = sorted_dirac_positions[ n ];
            const TF w1 = sorted_dirac_weights[ n ];
    
            const TF i1 = ( d0 + d1 + ( w0 - w1 ) / ( d1 - d0 ) ) / 2;
            func( d0, w0, n - 1, o0, d1 - d0, i0, i1, num_thread );
    
            o0 = d1 - d0;
            i0 = i1;
            d0 = d1;
            w0 = w1;
        }

        TF i1 = numeric_limits<TF>::max();
        TF o1 = numeric_limits<TF>::max();
        if ( end < nb_diracs() ) {
            const TF dn = sorted_dirac_positions[ end ];
            const TF wn = sorted_dirac_weights[ end ];
            i1 = ( d0 + dn + ( w0 - wn ) / ( dn - d0 ) ) / 2;
            o1 = dn - d0;
        }

        func( d0, w0, end - 1, o0, o1, i0, i1, num_thread );
    } );
}

DTP void UTP::for_each_cell( auto &&func ) const {
    using namespace std;

    TF i0 = numeric_limits<TF>::lowest();
    TF d0 = sorted_dirac_positions[ 0 ];
    TF w0 = sorted_dirac_weights[ 0 ];
    TF od = numeric_limits<TF>::max();
    for( PI n = 1; n < nb_diracs(); ++n ) {
        const TF d1 = sorted_dirac_positions[ n ];
        const TF w1 = sorted_dirac_weights[ n ];

        const TF i1 = ( d0 + d1 + ( w0 - w1 ) / ( d1 - d0 ) ) / 2;
        func( d0, w0, n - 1, od, d1 - d0, i0, i1 );

        od = d1 - d0;
        i0 = i1;
        d0 = d1;
        w0 = w1;
    }

    func( d0, w0, nb_diracs() - 1, od, numeric_limits<TF>::max(), i0, numeric_limits<TF>::max() );
}

DTP TF UTP::density_primitive( TF x ) const {
    if ( x >= density_values.size() ) return density_integrals.back();
    if ( x < 0 ) return 0;

    PI ix( x );
    TF fx = x - ix;
    return density_integrals[ ix ] + fx * density_values[ ix ];
}

DTP TF UTP::density_integral( TF x0, TF x1 ) const {
    return density_primitive( x1 ) - density_primitive( x0 );
}

DTP TF UTP::density_value( TF x ) const {
    if ( x < 0 || x >= density_values.size() ) return 0;
    return density_values[ x ];
}

DTP std::tuple<SymmetricBandMatrix<TF>,Vec<TF>,TF,bool> UTP::newton_system_ap( TF eps ) {
    using std::max;
    using std::min;

    // M
    SymmetricBandMatrix<TF> M( FromSize(), nb_diracs() );    
    Vec<TF> R = normalized_integrals();
    for( PI n = 0; n < nb_diracs(); ++n ) {
        TF &ref_weight = sorted_dirac_weights[ n ];
        TF old_weight = ref_weight;
        ref_weight += eps;

        Vec<TF> A = normalized_integrals();
        for( PI m = max( n, 1ul ) - 1; m <= n; ++m )
            M( m, n ) += 2 * ( A[ m ] - R[ m ] ) / eps;

        ref_weight = old_weight;
    }

    // V
    Vec<TF> V = 2 * ( sorted_dirac_masses - R );

    return { M, V, norm_2_p2( V ), false };
}


DTP std::tuple<SymmetricBandMatrix<TF>,Vec<TF>,TF,bool> UTP::newton_system() {
    using namespace std;

    const PI mul_thread = max( 1, 32 / sizeof( TF ) );
    const PI nb_threads = thread_pool.nb_threads();

    Vec<TF> has_arc = Vec<TF>::fill( mul_thread * nb_threads, 0 );
    Vec<TF> R = Vec<TF>::fill( mul_thread * nb_threads, 0 );
    SymmetricBandMatrix<TF> M( FromSize(), nb_diracs() );
    Vec<TF> V( FromSize(), nb_diracs() );
    for_each_cell_mt( [&]( TF dirac_position, TF dirac_weight, PI num_dirac, TF ldist, TF rdist, TF c0, TF c1, PI num_thread ) {
        const TF rd = sqrt( dirac_weight );
        const TF b0 = dirac_position - rd;
        const TF b1 = dirac_position + rd;

        TF v = sorted_dirac_masses[ num_dirac ];
        TF mc, mr;
        if ( c0 > c1 ) {
            TODO;
        } else {
            if ( b0 > c0 ) { // ball cut on the left
                const TF sqrtw = sqrt( dirac_weight );
                if ( b1 < c1 ) { // ball cut on the right
                    mc = ( density_value( b0 ) + density_value( b1 ) ) / sqrtw;
                    mr = 0;

                    v -= density_integral( b0, b1 );
                } else { // interface on the right
                    const TF cr = density_value( c1 ) / rdist;
                    mc = density_value( b0 ) / sqrtw + cr;
                    mr = - cr;
                    
                    v -= density_integral( b0, c1 );
                }
                has_arc[ mul_thread * num_thread ] = true;
            } else { // interface on the left
                if ( b1 < c1 ) { // ball cut on the right
                    const TF cl = density_value( c0 ) / ldist;
                    const TF sqrtw = sqrt( dirac_weight );
                    mc = cl + density_value( b1 ) / sqrtw;
                    mr = 0;
                    
                    v -= density_integral( c0, b1 );
                    has_arc[ mul_thread * num_thread ] = true;
                } else { // interface on the right
                    const TF cl = density_value( c0 ) / ldist;
                    const TF cr = density_value( c1 ) / rdist;
                    mc = cl + cr;
                    mr = - cr;

                    v -= density_integral( c0, c1 );
                }
            }
        }        

        v *= 2;

        if ( num_dirac + 1 < nb_diracs() )
            M( num_dirac, num_dirac + 1 ) = mr;
        M( num_dirac, num_dirac ) = mc;
        V[ num_dirac ] = v;
        R[ mul_thread * num_thread ] += v * v;
    } );

    return { M, V, sum( R ), sum( has_arc ) };
}

DTP Vec<TF> UTP::newton_dir() {
    auto sys = newton_system();

    auto &M = std::get<0>( sys );
    auto &V = std::get<1>( sys );
    bool ha = std::get<3>( sys );

    // P( M, V );
    M( 0, 0 ) += 1;

    if ( ! ha ) {
        M( 0, 0 ) += 1;
        TODO;
    }

    // solve
    return M.solve( V );
}

DTP Vec<TF> UTP::normalized_integrals() {
    using namespace std;
    Vec<TF> res( FromSize(), nb_diracs() );
    for_each_cell( [&]( TF dirac_position, TF dirac_weight, PI num_dirac, TF, TF, TF c0, TF c1 ) {
        const TF rd = sqrt( dirac_weight );
        if ( c0 > c1 ) {
            res[ num_dirac ] = - density_integral(
                max( dirac_position - rd, c1 ),
                min( dirac_position + rd, c0 ) 
            );
        } else {
            res[ num_dirac ] = density_integral(
                max( dirac_position - rd, c0 ),
                min( dirac_position + rd, c1 ) 
            );
        }        
    } );
    return res;
}

DTP TF UTP::normalized_error() {
    using namespace std;

    const PI mul_thread = max( 1, 32 / sizeof( TF ) );
    const PI nb_threads = thread_pool.nb_threads();
    auto res = Vec<TF>::fill( mul_thread * nb_threads, 0 );
    for_each_cell_mt( [&]( TF dirac_position, TF dirac_weight, PI num_dirac, TF, TF, TF c0, TF c1, PI num_thread ) {
        const TF rd = sqrt( dirac_weight );
        TF loc;
        if ( c0 > c1 ) {
            loc = - density_integral(
                max( dirac_position - rd, c1 ),
                min( dirac_position + rd, c0 ) 
            );
        } else {
            loc = density_integral(
                max( dirac_position - rd, c0 ),
                min( dirac_position + rd, c1 ) 
            );
        }        
        res[ mul_thread * num_thread ] += pow( sorted_dirac_masses[ num_dirac ] - loc, 2 );
    } );
    return sum( res );
}

DTP bool UTP::update_weights() {
    using namespace std;

    return true;
}

DTP void UTP::solve() {
    //
    for( PI num_iter = 0; num_iter < 30; ++num_iter ) {        
        Vec<TF> olw = sorted_dirac_weights;

        P( num_iter );

        // newton system
        auto sys = newton_system();
        auto &M = std::get<0>( sys );
        auto &V = std::get<1>( sys );
        auto er = std::get<2>( sys );
        bool ha = std::get<3>( sys );
        M( 0, 0 ) += 1;
    
        // find a search dir
        Vec<TF> dir = M.solve( V );

        // find a first relaxation coeff
        for( TF a = 1; ; a /= 2 ) {
            if ( a == 1e-20 )
                throw std::runtime_error( "bad direction" );
            sorted_dirac_weights = olw + a * dir;
            TF nr = normalized_error();
            P( a, er, nr );
            if ( nr < er ) {
                break;
            }
        }
    }

    //set_density_contrast( 1e-6 );
}

DTP TF UTP::normalized_dirac_convolution( TF normalized_pos, TF convolution_width ) const {
    using namespace std;

    TF res = 0;
    for( PI i = 0; i < nb_diracs(); ++i )
        res += sorted_dirac_masses[ i ] * exp( - pow( ( normalized_pos - sorted_dirac_positions[ i ] ) / convolution_width, 2 ) );
    return res;
}

DTP TF UTP::normalized_density_value( TF normalized_pos ) const {
    return normalized_pos >= 0 && normalized_pos < density_values.size() ? density_values[ normalized_pos ] : 0;
}

DTP void UTP::set_density_contrast( TF max_ratio ) {
    // first take
    const PI total_size = add_left + original_density_values.size() + add_right;
    const TF min_value = max_of_original_density_values * max_ratio;
    need_lower_contrast_ratio = add_left || add_right;
    density_integrals.resize( total_size + 1 );
    density_values.resize( total_size );

    TF sum_of_density_values = 0;
    for( PI i = 0; i < add_left; ++i ) {
        density_integrals[ i ] = sum_of_density_values;
        sum_of_density_values += min_value;
        density_values[ i ] = min_value;
    }

    for( PI i = 0; i < original_density_values.size(); ++i ) {
        TF dv = original_density_values[ i ];
        if ( dv < min_value ) {
            need_lower_contrast_ratio = true;
            dv = min_value;
        }

        density_integrals[ add_left + i ] = sum_of_density_values;
        density_values[ add_left + i ] = dv;
        sum_of_density_values += dv;
    }

    for( PI i = add_left + original_density_values.size(); i < total_size; ++i ) {
        density_integrals[ i ] = sum_of_density_values;
        sum_of_density_values += min_value;
        density_values[ i ] = min_value;
    }

    density_integrals.back() = sum_of_density_values;

    // normalization
    const TF coeff = sum_of_dirac_masses / ( sum_of_density_values * global_mass_ratio );
    for( TF &v : density_integrals )
        v *= coeff;
    for( TF &v : density_values )
        v *= coeff;
}

#undef DTP
#undef UTP


} // namespace usdot
