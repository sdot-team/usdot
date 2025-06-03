#pragma once

#include "common_macros.h"
#include "common_types.h"
#include <vector>

namespace usdot {

/**
*/
template<class TF>
class IndexedPositions {
public:
    struct Index           { PI n; TF f, x0, x1; T_T void display( T &ds ) const { ds.start_array(); ds << n << f << x0 << x1; ds.end_array(); } };
    using  PInds           = std::array<Index,2>;
    using  VF              = std::vector<TF>;
    using  VI              = std::vector<PI>;
 
    /**/   IndexedPositions( const VF &positions = {}, TF mul_nb_indices = 1 );
  
    #ifdef TL_DISPLAYER_IS_DEFINED
    void   display         ( Displayer &ds ) const;
    #endif

    Index  bounded_index   ( TF x ) const;
    TF     index_f         ( TF x ) const;
    Index  index           ( TF x ) const;
    TF     min_x           () const { return positions.front(); }
    TF     max_x           () const { return positions.back(); }

    VF     positions;      ///<
    VI     beg_inds;       ///<
    TF     beg_x;          ///<
    TF     end_x;          ///<
    TF     mul_x;          ///<
};


#define DTP template<class TF>
#define UTP IndexedPositions<TF>

DTP UTP::IndexedPositions( const VF &positions, TF mul_nb_indices ) : positions( positions ) {
    using namespace std;
    if ( positions.empty() )
        return;

    const PI nb_inds = PI( ceil( positions.size() * mul_nb_indices ) );
    beg_x = positions.front();
    end_x = positions.back();

    const TF dx = end_x - beg_x;
    mul_x = nb_inds / dx;

    beg_inds.resize( nb_inds, positions.size() );
    for( PI i = 1; i < positions.size(); ++i ) {
        const PI x0 = PI( floor( ( positions[ i - 1 ] - beg_x ) * mul_x ) );
        const PI x1 = PI( ceil( ( positions[ i - 0 ] - beg_x ) * mul_x ) );
        for( PI n = x0; n < std::min( x1, nb_inds ); ++n )
            beg_inds[ n ] = std::min( beg_inds[ n ], i );
    }
}

DTP typename UTP::Index UTP::index( TF x ) const {
    if ( x < beg_x )
        return { 0, -1 };
    if ( x > end_x )
        return { 0, +2 };

    const PI ox = PI( ( x - beg_x ) * mul_x );
    if ( ox >= beg_inds.size() )
        return { 0, +2 };

    for( PI i = beg_inds[ ox ]; ; ++i ) {
        const TF x1 = positions[ i - 0 ];
        if ( x1 >= x ) {
            const TF x0 = positions[ i - 1 ];
            if ( const TF dx = x1 - x0 )
                return { i, ( x - x0 ) / dx, x0, x1 };
            return { i, 0, x0, x1 };
        }
    }
}

DTP TF UTP::index_f( TF x ) const {
    if ( x < beg_x )
        return 0;
    if ( x > end_x )
        return positions.size() - 1;

    const PI ox = PI( ( x - beg_x ) * mul_x );
    if ( ox >= beg_inds.size() )
        return positions.size() - 1;

    for( PI i = beg_inds[ ox ]; ; ++i ) {
        const TF x1 = positions[ i - 0 ];
        if ( x1 >= x ) {
            const TF x0 = positions[ i - 1 ];
            if ( const TF dx = x1 - x0 )
                return i - 1 + ( x - x0 ) / dx;
            return i - 0.5;
        }
    }
}

DTP typename UTP::Index UTP::bounded_index( TF x ) const {
    if ( x < beg_x )
        return { 0, 0 };
    if ( x > end_x )
        return { positions.size() - 1, 1 };

    const PI ox = PI( ( x - beg_x ) * mul_x );
    if ( ox >= beg_inds.size() )
        return { positions.size() - 1, 1 };

    for( PI i = beg_inds[ ox ]; ; ++i ) {
        const TF x1 = positions[ i - 0 ];
        if ( x1 >= x ) {
            const TF x0 = positions[ i - 1 ];
            if ( const TF dx = x1 - x0 )
                return { i, ( x - x0 ) / dx, x0, x1 };
            return { i, 0, x0, x1 };
        }
    }
}

#undef DTP
#undef UTP

} // namespace usdot
