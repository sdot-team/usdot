#include <tl/support/operators/argmin.h>
#include <tl/support/ASSERT.h>
#include <tl/support/ERROR.h>

#include "partial1D/GluedCells.h"
#include "InitialSolution.h"

// #include <matplotlibcpp.h>
// #include <utility>


InitialSolution::TF InitialSolution::p_to_x( PI density_ind, TF p ) const {
    TF de = density_ps[ density_ind + 1 ] - density_ps[ density_ind + 0 ];
    return de ?
        density_xs[ density_ind + 0 ] + ( density_xs[ density_ind + 1 ] - density_xs[ density_ind + 0 ] ) * ( p - density_ps[ density_ind + 0 ] ) / de :
        density_xs[ density_ind + 0 ];
}

InitialSolution::TF InitialSolution::x_to_p( PI density_ind, TF x ) const {
    TF de = density_xs[density_ind + 1] - density_xs[density_ind + 0];
    return de ?
        density_ps[ density_ind + 0 ] + ( density_ps[ density_ind + 1 ] - density_ps[ density_ind + 0 ] ) * ( x - density_xs[ density_ind + 0 ] ) / de :
        density_ps[ density_ind + 0 ];
}

InitialSolution::TF InitialSolution::barycenter( const Cell &cell ) const {
    if ( cell.first_density_ind == cell.last_density_ind )
        return ( cell.beg_x + cell.end_x ) / 2;

    TF sum_x = 0, sum_1 = 0;
    auto add_part = [&]( TF density_value, TF bx, TF ex ) {
        TF pond = density_value * ( ex - bx );
        sum_x += pond * ( ex + bx ) / 2;
        sum_1 += pond;
    };

    add_part( density_vs[ cell.first_density_ind ], cell.beg_x, density_xs[ cell.first_density_ind + 1 ] );
    for( PI density_ind = cell.first_density_ind + 1; density_ind < cell.last_density_ind; ++density_ind )
        add_part( density_vs[ density_ind ], density_xs[ density_ind + 0 ], density_xs[ density_ind + 1 ] );
    add_part( density_vs[ cell.last_density_ind ], density_xs[ cell.last_density_ind ], cell.end_x );

    return sum_x / sum_1;
}

Vec<InitialSolution::TF> InitialSolution::barycenters() const {
    Vec<TF> res( FromSize(), nb_diracs() );
    for_each_cell( [&]( const Cell &cell ) {
        res[ cell.dirac_num ] = barycenter( cell );
    } );
    return res;
}

void InitialSolution::sanity_check() const {
    PI nb_gc = 0, nb_cells = 0;
    for( GluedCells *gc = first_glued_cells; gc; gc = gc->next_gc ) {
        if ( gc->prev_gc )
            ASSERT( end_p( *gc->prev_gc ) < gc->beg_p );
        nb_cells += gc->nb_cells();
        ++nb_gc;

        // check min cost
        // TF minc = cost( *gc, 0, 0 );
        // for( TF o = -1e-2; o < 1e2; o += 1e-4 ) {
        //     if ( gc->beg_p + o >= 0 and end_p( *gc ) + o <= 1 and minc - 1e-6 > cost( *gc, 0, o ) ) {
        //         P( minc, cost( *gc, 0, o ), *gc, o );
        //         ERROR( "" );
        //     }
        // }

        // PI n = 0;
        // PI ind = gc.beg_ind, density_ind = beg_density_ind;
        // const TF base_beg_p = gc.beg_p + ( n + 0 ) * cell_size;
        // TF d0 = dirac_xs[ dirac_nums[ ind ] ] - p_to_x( density_ind, base_beg_p );

        // n = gc.nb_cells() - 1;
        // ind = gc.end_ind - 1;
        // const TF base_end_p = gc.beg_p + ( n + 1 ) * cell_size;
        // while ( density_ind + 2 < density_ps.size() && density_ps[ density_ind + 1 ] <= base_end_p )
        //     ++density_ind;
        // TF d1 = p_to_x( density_ind, base_end_p ) - dirac_xs[ dirac_nums[ ind ] ];

        // if ( ( d1 - d0 ) / ( d1 ) > 1e-6 ) {
        //     if ( gc.beg_p > 1e-3 && end_p( gc ) < 1 - 1e-3 ) {
        //         P( d0, d1, ( d1 - d0 ) / ( d1 ), gc.nb_cells(), gc.beg_p );
        //         TF minc = cost( gc, 0, 0 );
        //         for( TF o = -1e-2; o < 1e2; o += 1e-4 )
        //             if ( gc.beg_p + o >= 0 and end_p( gc ) + o <= 1 and minc - 1e-6 > cost( gc, 0, o ) ) {
        //                 P( minc, cost( gc, 0, o ), gc, o );
        //                 ERROR( "" );
        //             }
        //     }
        // }

    }
    ASSERT( nb_cells == dirac_xs.size() );
    // P( nb_gc );
}

void InitialSolution::push_touching_cells() {
    // const TF cell_size = mass_ratio / nb_diracs();
    // ASSERT( glued_cells.size() );
    // bool blocked = false;

    // while ( true ) {
    //     PI gs = glued_cells.size();
    //     if ( gs < 2 )
    //         break;

    //     GluedCells &last = glued_cells[ gs - 1 ];
    //     GluedCells &prev = glued_cells[ gs - 2 ];
    //     TF delta = prev.end_p - last.beg_p;
    //     if ( delta < 0 )
    //         break;

    //     // append the references in prev
    //     prev.end_ind = last.end_ind;

    //     // start with position assuming last is not moved
    //     prev.beg_p = last.end_p - prev.nb_cells() * cell_size;
    //     prev.end_p = last.end_p;

    //     //
    //     if ( ! blocked ) {
    //         optimize_position( prev );

    //         if ( prev.end_p > 1 ) {
    //             prev.beg_p = 1 - prev.nb_cells() * cell_size;
    //             prev.end_p = 1;
    //             blocked = true;
    //         }
    //     }

    //     //
    //     glued_cells.pop_back();
    // }

    // // < 0 ?
    // GluedCells &gc = glued_cells.back();
    // if ( gc.beg_p < 0 ) {
    //     gc.end_p = cell_size * gc.nb_cells();
    //     gc.beg_p = 0;
    // }
}

void InitialSolution::display( Displayer &ds ) const {
    ds.start_array();
    for( const GluedCells *gc = first_glued_cells; gc; gc = gc->next_gc )
        ds << *gc;
    ds.end_array();
}

void InitialSolution::add_off_x( Cell &cell, TF off_x ) const {
    if ( off_x == 0 )
        return;

    // update interval beginning
    cell.beg_x += off_x;

    if ( off_x > 0 )
        while ( cell.first_density_ind + 1 < density_xs.size() && density_xs[ cell.first_density_ind + 1 ] <= cell.beg_x )
            ++cell.first_density_ind;
    else
        while ( cell.first_density_ind && density_xs[ cell.first_density_ind ] > cell.beg_x )
            --cell.first_density_ind;
    cell.beg_p = x_to_p( cell.first_density_ind, cell.beg_x );

    // update interval ending
    cell.end_p = cell.beg_p + cell_size;

    if ( off_x > 0 )
        while ( cell.last_density_ind + 1 < density_ps.size() && density_ps[ cell.last_density_ind + 1 ] <= cell.end_p )
            ++cell.last_density_ind;
    else
        while ( cell.last_density_ind && density_ps[ cell.last_density_ind ] > cell.end_p )
            --cell.last_density_ind;
    cell.end_x = p_to_x( cell.last_density_ind, cell.end_p );
}

void InitialSolution::add_off_p( Cell &cell, TF off_p ) const {
    if ( off_p == 0 )
        return;

    // update p
    cell.beg_p += off_p;
    cell.end_p += off_p;
    
    // update x
    if ( off_p > 0 ) {
        while ( cell.first_density_ind + 1 < density_ps.size() && density_ps[ cell.first_density_ind + 1 ] < cell.beg_p )
            ++cell.first_density_ind;
        while ( cell.last_density_ind + 1 < density_ps.size() && density_ps[ cell.last_density_ind + 1 ] < cell.end_p )
            ++cell.last_density_ind;
    } else{
        while ( cell.first_density_ind && density_ps[ cell.first_density_ind ] > cell.beg_p )
            --cell.first_density_ind;
        while ( cell.last_density_ind && density_ps[ cell.last_density_ind ] > cell.end_p )
            --cell.last_density_ind;
    }

    cell.beg_x = p_to_x( cell.first_density_ind, cell.beg_p );
    cell.end_x = p_to_x( cell.last_density_ind, cell.end_p );
}

InitialSolution::TF InitialSolution::cost( const GluedCells &gc, TF off_x, TF off_p ) const {
    using std::pow;
    TF res = 0;
    // int( ( x - c )^2 * density_value )
    for_each_cell( gc, [&]( Cell cell ) {
        add_off_x( cell, off_x );
        add_off_p( cell, off_p );
    
        TF c = dirac_xs[ cell.dirac_num ];
        auto add_int = [&]( TF beg_x, TF end_x, PI density_ind ) {
            res += density_vs[ density_ind ] * ( pow( end_x - c, 3 ) - pow( beg_x - c, 3 ) ) / 3;
            // P( cell.beg_x, cell.end_x, density_vs[ density_ind ], beg_x, end_x );
        };
        if ( cell.first_density_ind == cell.last_density_ind ) {
            add_int( cell.beg_x, cell.end_x, cell.first_density_ind );
        } else {
            add_int( cell.beg_x, density_xs[ cell.first_density_ind + 1 ], cell.first_density_ind );
            for( PI ind = cell.first_density_ind + 1; ind < cell.last_density_ind; ++ind )
                add_int( density_xs[ ind + 0 ], density_xs[ ind + 1 ], ind );
            add_int( density_xs[ cell.last_density_ind ], cell.end_x, cell.last_density_ind );
        }
    } );
    return res;
}

void InitialSolution::optimize_position( GluedCells &gc, TF beg_off_p, TF end_off_p ) {
    using std::pow;

    // ensures that p stays between 0 and 1
    if ( gc.beg_p + beg_off_p < 0 )
        beg_off_p = - gc.beg_p;
    if ( end_p( gc ) + end_off_p > 1 )
        end_off_p = 1 - end_p( gc );

    // std::vector<TF> xs, ys;
    // for( TF x = beg_off_p; x <= end_off_p; x += 0.01 ) {
    //     TF p = cost( gc, 0, x + 1e-4 ) - cost( gc, 0, x );
    //     xs.push_back( x );
    //     ys.push_back( p );
    // }
    // matplotlibcpp::plot( xs, ys, "-" );
    // matplotlibcpp::show();
    // P( ys );
    // TODO;
    
    // helper function that returns 0 if the best offset is found (gc will be modified accordingly)
    PI beg_density_ind;
    auto test_off_p = [&]( TF off_p ) -> int {
        // beg_density_ind (gc.beg_density_ind after the off_p modification)
        beg_density_ind = gc.beg_density_ind;
        if ( off_p < 0 )
            while ( beg_density_ind && density_ps[ beg_density_ind ] > gc.beg_p + off_p )
                --beg_density_ind;
        if ( off_p > 0 )
            while ( beg_density_ind + 2 < density_ps.size() && density_ps[ beg_density_ind + 1 ] <= gc.beg_p + off_p )
                ++beg_density_ind;

        // get limits of "knowledge" (min and max of off_p for which we known the density values)
        TF min_off_p = beg_off_p;
        TF max_off_p = end_off_p;
        for( PI ind = gc.beg_ind, n = 0, density_ind = beg_density_ind; ind < gc.end_ind; ++ind, ++n ) {
            const TF base_beg_p = gc.beg_p + ( n + 0 ) * cell_size;
            const TF base_end_p = gc.beg_p + ( n + 1 ) * cell_size;

            // beg cell
            min_off_p = max( min_off_p, density_ps[ density_ind + 0 ] - base_beg_p );
            max_off_p = min( max_off_p, density_ps[ density_ind + 1 ] - base_beg_p );

            // end cell
            const TF end_p = off_p + base_end_p;
            while ( density_ind + 2 < density_ps.size() && density_ps[ density_ind + 1 ] <= end_p )
                ++density_ind;
            min_off_p = max( min_off_p, density_ps[ density_ind + 0 ] - base_end_p );
            max_off_p = min( max_off_p, density_ps[ density_ind + 1 ] - base_end_p );
        }

        // get cost for each bound of off_p
        TF cost_at_min_off_p = 0;
        TF cost_at_max_off_p = 0;
        for( PI ind = gc.beg_ind, n = 0, density_ind = beg_density_ind; ind < gc.end_ind; ++ind, ++n ) {
            const TF base_beg_p = gc.beg_p + ( n + 0 ) * cell_size;
            const TF base_end_p = gc.beg_p + ( n + 1 ) * cell_size;
            const TF dirac_x = dirac_xs[ dirac_nums[ ind ] ];

            // beg cell
            cost_at_min_off_p -= pow( p_to_x( density_ind, base_beg_p + min_off_p ) - dirac_x, 2 );
            cost_at_max_off_p -= pow( p_to_x( density_ind, base_beg_p + max_off_p ) - dirac_x, 2 );

            // end cell
            const TF end_p = base_end_p + ( min_off_p + max_off_p ) / 2;
            while ( density_ind + 2 < density_ps.size() && density_ps[ density_ind + 1 ] <= end_p )
                ++density_ind;
            cost_at_min_off_p += pow( p_to_x( density_ind, base_end_p + min_off_p ) - dirac_x, 2 );
            cost_at_max_off_p += pow( p_to_x( density_ind, base_end_p + max_off_p ) - dirac_x, 2 );
        }
 
        //
        if ( cost_at_min_off_p > 0 && cost_at_max_off_p > 0 ) {
            end_off_p = min_off_p;
            if ( end_off_p == beg_off_p ) {
                gc.beg_density_ind = beg_density_ind;        
                gc.beg_p += beg_off_p;        
                return 0;
            }
            return + 1;
        }

        //
        if ( cost_at_min_off_p < 0 && cost_at_max_off_p < 0 ) {
            beg_off_p = max_off_p;
            if ( end_off_p == beg_off_p ) {
                gc.beg_density_ind = beg_density_ind;        
                gc.beg_p += beg_off_p;        
                return 0;
            }
            return - 1;
        }

        const TF best_off_p = min_off_p + cost_at_min_off_p / ( cost_at_min_off_p - cost_at_max_off_p ) * ( max_off_p - min_off_p );
        gc.beg_density_ind = beg_density_ind;        
        gc.beg_p += best_off_p;        
        return 0;
    };

    int bt = test_off_p( beg_off_p );
    if ( bt == 0 )
        return;
    if ( bt > 0 ) {
        gc.beg_density_ind = beg_density_ind;        
        gc.beg_p += beg_off_p;        
        return;
    }

    int et = test_off_p( end_off_p );
    if ( et == 0 )
        return;
    if ( et < 0 ) {
        gc.beg_density_ind = beg_density_ind;        
        gc.beg_p += end_off_p;        
        return;
    }


    //
    for( TF b = beg_off_p, e = end_off_p, c = 0; ; ++c ) {
        ASSERT( c < 20 * density_xs.size() );

        const TF mid_off_p = ( beg_off_p + end_off_p ) / 2;
        int mt = test_off_p( mid_off_p );
        if ( mt == 0 ) {
            // PI n = 0;
            // PI ind = gc.beg_ind, density_ind = beg_density_ind;
            // const TF base_beg_p = gc.beg_p + ( n + 0 ) * cell_size;
            // TF d0 = dirac_xs[ dirac_nums[ ind ] ] - p_to_x( density_ind, base_beg_p );

            // n = gc.nb_cells() - 1;
            // ind = gc.end_ind - 1;
            // const TF base_end_p = gc.beg_p + ( n + 1 ) * cell_size;
            // while ( density_ind + 2 < density_ps.size() && density_ps[ density_ind + 1 ] <= base_end_p )
            //     ++density_ind;
            // TF d1 = p_to_x( density_ind, base_end_p ) - dirac_xs[ dirac_nums[ ind ] ];

            // if ( ( d1 - d0 ) / ( d1 ) > 1e-6 ) {
            //     if ( gc.beg_p > 1e-3 && end_p( gc ) < 1 - 1e-3 ) {
            //         P( d0, d1, ( d1 - d0 ) / ( d1 ), gc.nb_cells(), gc.beg_p );
            //         TF minc = cost( gc, 0, 0 );
            //         for( TF o = -1e-2; o < 1e2; o += 1e-4 )
            //             if ( gc.beg_p + o >= 0 and end_p( gc ) + o <= 1 and minc - 1e-6 > cost( gc, 0, o ) ) {
            //                 P( minc, cost( gc, 0, o ), gc, o );
            //                 ERROR( "" );
            //             }
            //     }
            // }
            return;
        }

        TF &r = mt < 0 ? beg_off_p : end_off_p;
        r = mid_off_p;
    }
}

InitialSolution make_initial_solution( Span<InitialSolution::TF> dirac_coords, Span<InitialSolution::TF> density_coords, Span<InitialSolution::TF> density_values, InitialSolution::TF mass_ratio ) {
    using TF = InitialSolution::TF;

    // check inputs
    ASSERT( density_coords.size() == density_values.size() + 1 );

    // base init
    InitialSolution res;
    res.density_xs = density_coords;
    res.density_vs = density_values;
    res.dirac_xs   = dirac_coords;
    res.mass_ratio = mass_ratio;

    // basic computed stuff
    res.cell_size = mass_ratio / res.nb_diracs();

    // dirac_nums 
    res.dirac_nums = { FromSizeAndFunctionOnIndex(), res.nb_diracs(), []( auto i ) { return i; } };
    std::sort( res.dirac_nums.begin(), res.dirac_nums.end(), [&]( PI a, PI b ) {
        return res.dirac_xs[ a ] < res.dirac_xs[ b ];
    } );

    // density_ps (normalized primitive of density, `p` standing as primitive)
    TF acc = 0;
    res.density_ps.resize( res.density_xs.size() );
    for( PI n = 0; n < res.density_xs.size(); ++n ) {
        res.density_ps[ n ] = acc;
        if ( n + 1 < res.density_xs.size() )
            acc += ( res.density_xs[ n + 1 ] - res.density_xs[ n + 0 ] ) * res.density_vs[ n + 0 ];
    }
    for( PI n = 0; n < res.density_ps.size(); ++n )
        res.density_ps[ n ] /= res.density_ps.back();

    // start with the assumption that the cells do not intersect
    GluedCells *first_gc = nullptr;
    GluedCells *last_gc = nullptr;
    PI density_ind = 0;
    TF curr_p = 0;
    for( PI ind = 0, d = 0; ind < res.dirac_nums.size(); ++ind ) {
        GluedCells *gc = res.pool.create<GluedCells>();
        gc->num_moving_phase = 0;

        if ( ind ) {
            last_gc->next_gc_to_test = gc;
            last_gc->next_gc = gc;
        } else
            first_gc = gc;
  
        gc->prev_gc = last_gc;
        last_gc = gc;

        // push the new cell in an inexact position
        gc->beg_density_ind = density_ind;
        gc->beg_ind = ind + 0;
        gc->end_ind = ind + 1;
        gc->beg_p = curr_p;

        // find the optimal solution for this cell without taking care of the intersections
        res.optimize_position( *gc, 0, 1 - res.cell_size );
    }

    res.first_glued_cells = first_gc;
    if ( last_gc ) {
        last_gc->next_gc_to_test = nullptr;
        last_gc->next_gc = nullptr;
    }

    // while we can find new aglomerates
    GluedCells *first_gc_to_test = first_gc; // cells that have moved during the previous iteration + the prev ones if not already in the set
    for( PI num_moving_phase = 0; ; ++num_moving_phase ) {
        // first pass: aglomerate the touching sets
        GluedCells *last_gc_to_test = nullptr;
        for( GluedCells *curr_gc = std::exchange( first_gc_to_test, nullptr ); curr_gc; curr_gc = curr_gc->next_gc_to_test ) {
            // off_update = how much this new aggregate can be moved
            curr_gc->off_update = 0;

            GluedCells *next_gc = curr_gc->next_gc;
            TF end_p = res.end_p( *curr_gc );
            if ( next_gc && end_p >= next_gc->beg_p ) {
                // register the agregate in the list of sets to be tested during the next iteration of the loop
                if ( first_gc_to_test )
                    last_gc_to_test->next_gc_to_test = curr_gc;
                else
                    first_gc_to_test = curr_gc;
                last_gc_to_test = curr_gc;

                // agregate curr_gc with next_gc
                do {
                    curr_gc->off_update += end_p - next_gc->beg_p;
                    curr_gc->end_ind = next_gc->end_ind;

                    if ( curr_gc->next_gc_to_test == next_gc )
                        curr_gc->next_gc_to_test = curr_gc->next_gc_to_test->next_gc_to_test;
                
                    end_p = res.end_p( *next_gc );

                    next_gc = next_gc->next_gc;
                    curr_gc->next_gc = next_gc;
                } while ( next_gc && end_p >= next_gc->beg_p );

                if ( next_gc )
                    next_gc->prev_gc = curr_gc;
            }
        }

        // if no new intersection, break the loop
        if ( ! first_gc_to_test )
            break;

        // else, close the list of moved sets of cells
        last_gc_to_test->next_gc_to_test = nullptr;

        // update the positions
        for( GluedCells *gc = first_gc_to_test; gc; gc = gc->next_gc_to_test )
            res.optimize_position( *gc, - gc->off_update, 0 );

        // register to gc_to_test list for the next phase
        last_gc_to_test = nullptr;
        for( GluedCells *gc = std::exchange( first_gc_to_test, nullptr ); gc; gc = gc->next_gc_to_test ) {
            //
            auto add_item = [&]( GluedCells *item ) {
                if ( item->num_moving_phase == num_moving_phase + 1 )
                    return;
                item->num_moving_phase = num_moving_phase + 1;

                if ( first_gc_to_test )
                    last_gc_to_test->next_gc_to_test = item;
                else
                    first_gc_to_test = item;
                last_gc_to_test = item;
            };

            //
            if ( gc->prev_gc )
                add_item( gc->prev_gc );
            add_item( gc );
        }
    }
 
    res.sanity_check();

    return res;
}
