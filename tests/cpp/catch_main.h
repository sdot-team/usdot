#include <catch2/benchmark/catch_benchmark_all.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>
#include <catch2/catch_test_macros.hpp>

// #include <vfs/support/make_Vec.h> // IWYU pragma: export
// #include <tl/support/string/to_string.h> // IWYU pragma: export
// #include <tl/support/P.h> // IWYU pragma: export

#define CHECK_REPR( A, B ) \
    CHECK( to_string( A, { .always_display_delimiters = true } ) == to_string( B, { .always_display_delimiters = true } ) )

#define CHECK_PROX( A, B, TOL ) \
    REQUIRE_THAT( double( A ), Catch::Matchers::WithinAbs( double( B ), TOL ) )

#define CHECK_PROX_IF_NOT_EMPTY( A, B, TOL ) \
    do { if ( ( A ).size() && ( B ).size() ) { \
        REQUIRE( ( A ).size() == ( B ).size() ); \
        for( PI i = 0, s = ( A ).size(); i < s; ++i ) \
            CHECK_THAT( double( A[ i ] ), Catch::Matchers::WithinAbs( double( B[ i ] ), TOL ) ); \
    } } while( false )

//using namespace Vfs;

//#define CHECK_REPR_F( A, B, F ) \
//    CHECK( to_string( A, { .filter = F, .rec = 1 } ) == to_string( B, { .filter = F, .rec = 1 } ) )

