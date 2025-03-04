#pragma once

#include <tl/support/operators/norm_inf.h>
#include <tl/support/operators/norm_2.h>
#include <tl/support/operators/sum.h>
#include <tl/support/operators/pow.h>
#include <tl/support/operators/abs.h>
#include <tl/support/operators/max.h>

#include <tl/support/string/va_string.h>
#include <tl/support/P.h>

#include "Convolution/InvX2Convolution.h"
#include "Density/SumOfDensities.h"
#include "Density/BoundedDensity.h"
#include "Density/Lebesgue.h"
#include "Extrapolation.h"
#include "CdfSolver.h"
#include "Solver.h"

namespace usdot {

    
#define DTP template<class TF>
#define UTP Solver<TF>

DTP UTP::Solver( RcPtr<PowerDiagram<TF>> power_diagram, RcPtr<Density<TF>> density, TF global_mass_ratio ) : global_mass_ratio( global_mass_ratio ), power_diagram( power_diagram ), density( density ) {
    convolution_factory = []( TF width ) -> RcPtr<Convolution<TF>> { return new InvX2Convolution<TF>( width ); };
    relative_mass_ratios = { FromSizeAndItemValue(), power_diagram->nb_cells(), 1 };
}

DTP void UTP::find_approx_weights_using_cdf( InitializeWeightsPrm parms ) {
    Vec<TF> mass_ratios = global_mass_ratio / sum( relative_mass_ratios ) * relative_mass_ratios;
    CdfSolver cs( density->cdf_approximation( 1e-4 / power_diagram->nb_cells() ), power_diagram->sorted_seed_coords, mass_ratios );
    cs.solve( power_diagram->sorted_seed_weights );
}

DTP void UTP::find_exact_weights_using_cdf( InitializeWeightsPrm parms ) {
    if ( power_diagram->nb_cells() == 0 )
        return;
    
        // for nown this procedure only works with global_mass_ratio == 1
    if ( global_mass_ratio < 1 ) {
        TODO;
    }

    // variables used to find boundaries
    const auto mass_coeff = density->mass() / sum( relative_mass_ratios ) * global_mass_ratio;
    auto diter = density->iterator();
    TF c1 = diter->mass();
    TF y0 = 0;
    TF c0 = 0;
 
    // find boundaries
    Vec<TF> bnds;
    bnds << diter->x0;
    for( PI n = 0; n < power_diagram->nb_cells(); ++n ) {
        TF y1 = y0 + relative_mass_ratios[ power_diagram->sorted_seed_nums[ n ] ] * mass_coeff;
        while ( c1 < y1 ) {
            diter->move_forward();
            c0 = std::exchange( c1, c1 + diter->mass() );
        }

        bnds << diter->inv_cdf( y1 - c0 );

        y0 = y1;
    }

    // first cell
    power_diagram->sorted_seed_weights[ 0 ] = 0;
    TF w = 0;
    
    // dc = current seed coord
    TF dc = power_diagram->sorted_seed_coords[ 0 ];
    TF w_to_add = 2 * max( pow( bnds[ 1 ] - dc, 2 ), pow( dc - bnds[ 0 ], 2 ) );
    for( PI n = 1; n < power_diagram->nb_cells(); ++n ) {
        // update w (dp = previous seed coord)
        TF dp = std::exchange( dc, power_diagram->sorted_seed_coords[ n ] );
        w += ( dp + dc - 2 * bnds[ n ] ) * ( dc - dp );
        
        // check that the ball is not goind to cut the cell
        const TF min_w = 2 * max( pow( bnds[ n ] - dc, 2 ), pow( dc - bnds[ n - 1 ], 2 ) );
        w_to_add = max( w_to_add, min_w - w );

        // store the result
        power_diagram->sorted_seed_weights[ n ] = w;
    }

    if ( w_to_add > 0 )
        for( PI n = 0; n < power_diagram->nb_cells(); ++n )
            power_diagram->sorted_seed_weights[ n ] += w_to_add;
}

DTP bool UTP::converged() const {
    if ( max_mass_ratio_error_history.empty() )
        return false;
    return max_mass_ratio_error_history.back() <= max_mass_ratio_error_target;
}

DTP void UTP::update_convex_hull_density_ratio( UpdateConvexHullDensityRatioPrm parms ) {
    if ( parms.target_value && *parms.target_value == convex_hull_density_ratio )
        return;

    if ( parms.epsilon ) {
        TF base_convex_hull_density_ratio = convex_hull_density_ratio;
        Vec<Vec<TF>> weights;
        weights << power_diagram->sorted_seed_weights;
        for( PI i = 1; i <= parms.polynomial_order; ++i ) {
            convex_hull_density_ratio = base_convex_hull_density_ratio - i * parms.epsilon;
            update_weights_using_newton();

            weights << power_diagram->sorted_seed_weights;
        }

        ASSERT( parms.target_value );

        power_diagram->sorted_seed_weights = extrapolation( weights, ( base_convex_hull_density_ratio - *parms.target_value ) / parms.epsilon );
        convex_hull_density_ratio = *parms.target_value;
        return;
    }

    ASSERT( parms.target_value );
    convex_hull_density_ratio = *parms.target_value;
}

DTP void UTP::update_convolution_width( UpdateConvolutionWidthPrm parms ) {
    TODO;
}

DTP UTP::State UTP::get_state() const {
    return {
        .convex_hull_density_ratio = convex_hull_density_ratio, 
        .convolution_width = convolution_width,
        .weights = power_diagram->sorted_seed_weights
    };
}

DTP void UTP::set_state( const State &state ) {
    convex_hull_density_ratio = state.convex_hull_density_ratio;
    convolution_width = state.convolution_width;
    power_diagram->sorted_seed_weights = state.weights;
}

DTP RcPtr<Density<TF>> UTP::current_density() {
    // convoluted
    if ( _convoluted_density_width != convolution_width || ! _convoluted_density ) {
        if ( convolution_width ) {
            TF min_x = min( power_diagram->sorted_seed_coords.front(), density->min_x() );
            TF max_x = min( power_diagram->sorted_seed_coords.back(), density->max_x() );
            _convoluted_density = new BoundedDensity<TF>( 
                density->convoluted( convolution_factory( convolution_width ) ),
                min_x,
                max_x
            );
            _convoluted_density_width = convolution_width;
        } else
            _convoluted_density = density;
    }

    // 
    if ( convex_hull_density_ratio ) {
        TF min_x = min( power_diagram->sorted_seed_coords.front(), _convoluted_density->min_x() );
        TF max_x = min( power_diagram->sorted_seed_coords.back(), _convoluted_density->max_x() );
        return new SumOfDensities<TF>( { 
            { convex_hull_density_ratio, new Lebesgue<TF>( min_x, max_x ) },
            { 1 - convex_hull_density_ratio, _convoluted_density }
        } );            
    }

    return _convoluted_density;
}

DTP int UTP::update_weights_using_newton( UpdateWeightsPrm parms ) {
    using namespace std;

    // used density
    RcPtr<Density<TF>> de = current_density();
    const TF mass_de = de->mass();

    // common variables
    SymmetricBandMatrix<TF> M( FromSize(), power_diagram->nb_cells() );
    Vec<TF> V( FromSize(), power_diagram->nb_cells() );
    PI nb_arcs;

    // history
    norm_2_rhs_history.clear();

    // until convergence
    Vec<TF> old_sorted_seed_weights;
    for( PI nb_iter = 0; ; ++nb_iter ) {
        if ( nb_iter == parms.max_newton_iterations ) {
            if ( verbosity >= 2 )
                P( "max_newton_iterations reached" );
            return 1;
        }

        // the the newton system
        V = global_mass_ratio * relative_mass_ratios / sum( relative_mass_ratios );
        M.fill_with( 0 );
        nb_arcs = 0;

        power_diagram->get_newton_system( M, V, nb_arcs, *de, 1 / mass_de );

        //
        TF norm_2_rhs = norm_2( V );
        if ( verbosity >= 2 )
            P( norm_2_rhs );

        // check the system
        TF mid = M( 0, 0 );
        TF mad = mid;
        for( PI i = 0; i < power_diagram->nb_cells(); ++i ) {
            // mav = std::max( mav, V[ i ] - target_mass_ratios[ i ] );
            mid = min( mid, abs( M( i, i ) ) );
            mad = max( mad, abs( M( i, i ) ) );
            if ( isnan( M( i, i ) ) )
                mid = 0;
        }

        
        // stop if system is in bad shape
        if ( mid == 0 || ( norm_2_rhs_history.size() && norm_2_rhs > norm_2_rhs_history.back() ) ) { // || mid / mad <= std::numeric_limits<TF>::epsilon() * power_diagram->nb_cells() TODO: a more precise criterion
            if ( nb_iter == 0 ) {
                if ( verbosity >= 2 )
                    P( "bad initialization" );
                return 2;
            }

            if ( verbosity >= 2 ) 
                P( "backtracking", mid, mad );
            power_diagram->sorted_seed_weights = 0.5 * old_sorted_seed_weights + 0.5 * power_diagram->sorted_seed_weights;
            continue;
        }

        //
        if ( nb_arcs == 0 )
            M( 0, 0 ) += 1;

        // solve
        Vec<TF> R = M.solve( V );

        // update weights
        old_sorted_seed_weights = power_diagram->sorted_seed_weights;
        power_diagram->sorted_seed_weights += R;

        // history
        TF max_mass_ratio_error = norm_inf( V / relative_mass_ratios );
        TF norm_2_residual = norm_2( R );

        max_mass_ratio_error_history << max_mass_ratio_error;
        norm_2_residual_history << norm_2_residual;
        norm_2_rhs_history << norm_2_rhs;

        if ( isnan( norm_2_residual ) ) {
            if ( verbosity >= 2 )
                P( "nan in the residual" );
            return 3;

        }

        // log
        if ( iteration_callback )
            iteration_callback();

        //
        if ( nb_iter >= parms.min_newton_iterations && converged() )
            break;     
    }

    return 0;
}

DTP void UTP::solve() {
    // not partial -> solve using cdf and that's it
    if ( global_mass_ratio == 1 ) {
        find_exact_weights_using_cdf();
        return;
    }

    // else, try with approx cdf + newton steps
    find_approx_weights_using_cdf();
    int res = update_weights_using_newton();
    if ( res == 0 )
        return;

    // else, solve with the characteristic of a convex hull
    convex_hull_density_ratio = 1;
    find_approx_weights_using_cdf();
    if ( update_weights_using_newton() )
        throw std::runtime_error( "not expected to fail here" );

    // then make a blend toward the required density
    for( TF ratio = 1; ratio > 1e-20; ) {
        P( ratio );

        // try 0, then 0.5, 0.25, ...
        Vec<TF> old_weights = power_diagram->sorted_seed_weights;
        for( TF next_ratio = ratio / 1e3; ; next_ratio = ratio - ( ratio - next_ratio ) / 2 ) {
            if ( next_ratio == ratio )
                throw std::runtime_error( "stuck at ratio convex_hull" );

            update_convex_hull_density_ratio( { .target_value = next_ratio, .epsilon = 0 } );
            if ( update_weights_using_newton() == 0 ) {
                ratio = next_ratio;
                break;
            }

            power_diagram->sorted_seed_weights = old_weights;
        }
    }
}

#undef DTP
#undef UTP


} // namespace usdot
