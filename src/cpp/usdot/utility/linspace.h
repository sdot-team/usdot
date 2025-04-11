#pragma once

#include <vector>

namespace usdot {

template<class TF> std::vector<TF> linspace( TF beg, TF end, std::size_t size, bool with_end = true ) {
    std::vector<TF> res( size );
    for( std::size_t i = 0; i < size; ++i )
        res[ i ] = beg + ( end - beg ) * i / ( size - with_end );
    return res;
}

template<class TF> std::vector<TF> cellspace( TF beg, TF end, std::size_t size, TF pos_in_cell = 0.5 ) {
    std::vector<TF> res( size );
    for( std::size_t i = 0; i < size; ++i ) {
        TF a = beg + ( end - beg ) * ( i + 0 ) / size;
        TF b = beg + ( end - beg ) * ( i + 1 ) / size;
        res[ i ] = a + ( b - a ) * pos_in_cell;
    }
    return res;
}

template<class TF> std::vector<TF> fill( std::size_t size, TF value ) {
    return std::vector<TF>( size, value );
}

template<class TF> std::vector<TF> iota( std::size_t size ) {
    std::vector<TF> res( size );
    for( std::size_t i = 0; i < size; ++i )
        res[ i ] = i;
    return res;
}

} // namespace usdot
