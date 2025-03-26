#pragma once

#include <tl/support/string/to_string.h>
#include <tl/support/containers/Vec.h>
#include <fstream>

template<class TF>
void glot( Vec<TF> xs, auto &&...funcs ) {
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

template<class TF>
void glot_vec( Vec<TF> xs, auto... ys ) {
    std::ofstream fs( "glot.py" );
    fs << "from matplotlib import pyplot\n";
    auto pf = [&]( auto &&y ) {
        fs << "pyplot.plot( " << to_string( xs ) << ", " << to_string( y ) << " )\n";
    };
    ( pf( ys ), ... );
    fs << "pyplot.show()\n";
}

