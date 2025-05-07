#pragma once

// #include <tl/support/string/to_string.h>
// #include <tl/support/containers/Vec.h>
#include <fstream>
#include <vector>

namespace usdot {
    
template<class Func>
void glot_stream( const Func &func ) {
    std::ofstream fs( "glot.py" );
    fs << "from matplotlib import pyplot\n";
    func( fs );
    fs << "pyplot.show()\n";
}

template<class TF,class ...Funcs>
void glot( std::vector<TF> xs, Funcs &&...funcs ) {
    std::ofstream fs( "glot.py" );
    fs << "from matplotlib import pyplot\n";
    auto pf = [&]( auto &&func ) {
        fs << "pyplot.plot( [ ";
        for( auto x : xs )
            fs << x << ", ";
        fs << " ], [ ";
        for( auto x : xs )
            fs << func( x ) << ", ";
        fs << " ] )\n";
    };
    ( pf( funcs ), ... );
    fs << "pyplot.show()\n";
}

// template<class TF>
// void glot_vec( std::vector<TF> xs, auto... ys ) {
//     std::ofstream fs( "glot.py" );
//     fs << "from matplotlib import pyplot\n";
//     std::size_t cpt = 0;
//     auto pf = [&]( auto &&y ) {
//         fs << "pyplot.plot( " << to_string( xs ) << ", " << to_string( y ) << ", label='" << cpt++ << "' )\n";
//     };
//     ( pf( ys ), ... );
//     fs << "pyplot.legend()\n";
//     fs << "pyplot.show()\n";
// }

template<class... Funcs>
void glot_vec_ys( Funcs &&...ys ) {
    std::ofstream fs( "glot.py" );
    fs << "from matplotlib import pyplot\n";
    int cpt = 0;
    auto pf = [&]( auto &&y ) {
        fs << "pyplot.plot( [ ";
        for( std::size_t i = 0; i < y.size(); ++i )
            fs << i << ", ";
        fs << " ], [ ";
        for( auto v : y )
            fs << v << ", ";
        fs << " ], label='" << cpt++ << "' )\n";
    };
    ( pf( ys ), ... );
    if ( cpt > 1 )
        fs << "pyplot.legend()\n";
    fs << "pyplot.show()\n";
}

} // namespace usdot
