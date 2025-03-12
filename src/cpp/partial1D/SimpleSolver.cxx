#pragma once

// #include <tl/support/operators/norm_inf.h>
// #include <tl/support/operators/norm_2.h>
// #include <tl/support/operators/sum.h>
// #include <tl/support/operators/pow.h>
// #include <tl/support/operators/abs.h>
// #include <tl/support/string/va_string.h>

#include <tl/support/operators/max.h>
#include <tl/support/ASSERT.h>
#include <tl/support/P.h>
#include "SimpleSolver.h"

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

    P( sorted_dirac_positions );

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

DTP void UTP::solve() {
    set_density_contrast( 1e-6 );
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
    const PI total_size = add_left + original_density_values.size() + add_right;

    // first take
    density_integrals.resize( total_size + 1, 0 );
    density_values.resize( total_size );

    const TF min_value = max_of_original_density_values * max_ratio;
    need_lower_contrast_ratio = false;
    
    for( PI i = 0; i < add_left; ++i )
        density_values[ i ] = min_value;
    
    TF sum_of_density_values = 0;
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

    for( PI i = add_left + original_density_values.size(); i <= total_size; ++i ) {
        density_integrals[ i ] = sum_of_density_values;
        density_values[ i ] = min_value;
    }

    // normalization
    const TF coeff = sum_of_dirac_masses / ( sum_of_density_values * global_mass_ratio );
    for( PI i = 0; i < density_values.size(); ++i ) {
        density_integrals[ i ] *= coeff;
        density_values[ i ] *= coeff;
    }
}

#undef DTP
#undef UTP


} // namespace usdot
