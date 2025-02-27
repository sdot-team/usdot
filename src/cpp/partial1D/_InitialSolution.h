#pragma once

#include <tl/support/memory/BumpPointerPool.h>
#include "CellInInitialSolution.h"
#include "GluedCells.h"

/**
 * Piecewise constant density
 */
class InitialSolution {
public:
    using                  Cell                  = CellInInitialSolution;
    using                  TF                    = Cell::TF;
    using                  PI                    = Cell::PI;
     
    void                   for_each_cell         ( const GluedCells &gc, const auto &func ) const; ///< func( cell )
    void                   for_each_cell         ( const auto &func ) const; ///< func( cell )
    void                   sanity_check          () const;
    Vec<TF>                barycenters           () const;
    TF                     barycenter            ( const Cell &cell ) const;
    PI                     nb_diracs             () const { return dirac_xs.size(); }
    void                   display               ( Displayer &ds ) const;

    friend InitialSolution make_initial_solution ( Span<InitialSolution::TF> dirac_coords, Span<InitialSolution::TF> density_coords, Span<InitialSolution::TF> density_values, InitialSolution::TF mass_ratio  );

private:
    void                   push_touching_cells   ();
    void                   optimize_position     ( GluedCells &gc, TF beg_p, TF end_p );
    void                   add_off_x             ( Cell &cell, TF off_x ) const;
    void                   add_off_p             ( Cell &cell, TF off_p ) const;
    Vec<TF,3>              der_cost              ( const GluedCells &gc, TF off_x );
    TF                     end_p                 ( const GluedCells &gc ) const { return gc.beg_p + cell_size * gc.nb_cells(); }
    TF                     cost                  ( const GluedCells &gc, TF off_x, TF off_p ) const;

    TF                     p_to_x                ( PI density_ind, TF p ) const;
    TF                     x_to_p                ( PI density_ind, TF x ) const;

    // input   
    TF                     mass_ratio;           ///<
    Span<TF>               density_xs;           ///<
    Span<TF>               density_vs;           ///<
    Span<TF>               dirac_xs;             ///<
       
    // output   
    GluedCells*            first_glued_cells;    ///<
    Vec<PI>                dirac_nums;           ///< num dirac vs num cell
    Vec<TF>                density_ps;           ///< primitive value
    TF                     cell_size;            ///< in the p axis
    BumpPointerPool        pool;                 ///<
};

void InitialSolution::for_each_cell( const GluedCells &gc, const auto &func ) const {
    PI density_ind = 0;
    for( PI n = 0, t = gc.nb_cells(); n < t; ++n ) {
        Cell cell;

        cell.dirac_num = dirac_nums[ gc.beg_ind + n ];
        cell.ind = gc.beg_ind + n;

        cell.beg_p = gc.beg_p + ( n + 0 ) * cell_size;
        cell.end_p = gc.beg_p + ( n + 1 ) * cell_size;

        while ( density_ind + 2 < density_ps.size() && density_ps[ density_ind + 1 ] <= cell.beg_p )
            ++density_ind;
        cell.first_density_ind = density_ind;
        cell.beg_x = p_to_x( density_ind, cell.beg_p );

        while ( density_ind + 2 < density_ps.size() && density_ps[ density_ind + 1 ] <= cell.end_p )
            ++density_ind;
        cell.last_density_ind = density_ind;
        cell.end_x = p_to_x( density_ind, cell.end_p );

        func( cell );
    }
}

void InitialSolution::for_each_cell( const auto &func ) const {
    PI density_ind = 0;
    for( const GluedCells *gc = first_glued_cells; gc; gc = gc->next_gc ) {
        for( PI n = 0, t = gc->nb_cells(); n < t; ++n ) {
            Cell cell;

            cell.dirac_num = dirac_nums[ gc->beg_ind + n ];
            cell.ind = gc->beg_ind + n;

            cell.beg_p = gc->beg_p + ( n + 0 ) * cell_size;
            cell.end_p = gc->beg_p + ( n + 1 ) * cell_size;

            while ( density_ind + 2 < density_ps.size() && density_ps[ density_ind + 1 ] <= cell.beg_p )
                ++density_ind;
            cell.first_density_ind = density_ind;
            cell.beg_x = p_to_x( density_ind, cell.beg_p );

            while ( density_ind + 2 < density_ps.size() && density_ps[ density_ind + 1 ] <= cell.end_p )
                ++density_ind;
            cell.last_density_ind = density_ind;
            cell.end_x = p_to_x( density_ind, cell.end_p );

            func( cell );
        }
    }
}

InitialSolution make_initial_solution( Span<InitialSolution::TF> dirac_coords, Span<InitialSolution::TF> density_coords, Span<InitialSolution::TF> density_values, InitialSolution::TF mass_ratio = 1 );
