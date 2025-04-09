#include <tl/support/operators/argmin.h>
#include <tl/support/containers/Vec.h>
#include <tl/support/P.h>
#include <matplotlibcpp.h>
#include <cmath>

using TF = FP80;

struct Integ {
    // Integ(  )    

    PI nb_cells() const {
        return dirac_ys.size();
    }

    TF cost( TF offset ) {
        Vec<TF> dirac_xs;
        for( TF y : dirac_ys ) 
            dirac_xs << y_to_x( y );

        TF res = 0;
        for( PI num_cell = 0; num_cell < nb_cells(); ++num_cell ) {
            TF b = offset + mass_ratio * ( num_cell + 0 ) / nb_cells();
            TF e = offset + mass_ratio * ( num_cell + 1 ) / nb_cells();
            PI step = 1e6;
            for( TF i = 0; i < step; ++i ) {
                TF d = ( e - b ) / step;
                TF r = ( i + 0.5 ) / step;
                TF y = ( 1 - r ) * b + r * e;
                res += std::pow( y_to_x( y ) - dirac_xs[ num_cell ], 2 ) * density( y_to_x( y ) ) * d;
            }
        }
        return res / 100;
    }

    TF density( TF x ) const {
        return ( x < 0.3 ) || ( x > 0.7 );
    }

    TF x_to_y( TF x ) const {
        if ( x < 0.3 ) 
            return x;
        if ( x < 0.7 ) 
            return 0.3;
        return 0.3 + ( x - 0.7 );
    }

    TF y_to_x( TF y ) const {
        if ( y < 0.3 )
            return y;
        return y + 0.4;
    }

    // TF density( TF x ) const {
    //     return x / 2;
    //     // return ( x < 0.3 ) || ( x > 0.7 );
    // }

    // TF x_to_y( TF x ) const {
    //     // if ( x < 0.3 ) 
    //     //     return x;
    //     // if ( x < 0.7 ) 
    //     //     return 0.3;
    //     // return 0.3 + ( x - 0.7 );
    //     return x * x;
    // }

    // TF y_to_x( TF y ) const {
    //     // return std::sqrt( y );
    //     // x + x * x = 2 * y => 
    //     // return 0.5 * ( std::sqrt( 8 * y + 1 ) - 1 );
    //     // if ( y < 0.3 )
    //     //     return y;
    //     // return y + 0.4;
    //     return std::sqrt( y );
    // }

    // TF density( TF x ) const {
    //     return 1;
    //     // return ( x < 0.3 ) || ( x > 0.7 );
    // }

    // TF x_to_y( TF x ) const {
    //     return x;
    // }

    // TF y_to_x( TF y ) const {
    //     return y;
    // }

    Vec<TF> dirac_ys{  .3, .4, .5 }; 
    TF mass_ratio = 0.5;
};

void add_f( std::vector<TF> ds, Str l ) {
    Integ integ;

    integ.dirac_ys = ds;
    std::vector<TF> x0, y0;
    for( TF offset = 0.0; offset < 0.47; offset += 0.01 ) {
        x0.push_back( offset );
        y0.push_back( integ.cost( offset + 1e-6 ) - integ.cost( offset ) );
    }
    // P( argmin( y0 ) );

    matplotlibcpp::plot( x0, y0, l );
}

int main() {
    add_f( {  .3, .4, .5 }, "-" );
    add_f( {  .3, .35, 0.45, .5 }, "--" );
    add_f( {  .3, .35, 0.4, 0.45, .5 }, "." );
    matplotlibcpp::show();
}

