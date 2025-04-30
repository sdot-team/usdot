#include "../../src/cpp/usdot/utility/TridiagonalPlus1LineSymmetricMatrix.h"
// #include "../../src/cpp/usdot/utility/TridiagonalSymmetricMatrix.h"
#include <tl/support/P.h> // IWYU pragma: export
#include "catch_main.h"
using namespace usdot;

TEST_CASE( "TridiagonalPlus1LineSymmetricMatrix", "" ) {
    TridiagonalPlus1LineSymmetricMatrix<double> M( 5 );

    for( PI i = 0; i < M.tridiag_values.size(); ++i )
        M.tridiag_values[ i ] = double( rand() ) / RAND_MAX + ( i % 2 == 0 );
    for( PI i = 0; i < M.full_matrix_size(); ++i )
        M.line_value( i ) = i < M.tridiag_size();

    std::vector<double> X( M.tridiag_size() );
    for( PI i = 0; i < X.size(); ++i )
        X[ i ] = double( rand() ) / RAND_MAX;
    P( M );
    P( X );

    M.inplace_ldlt_decomposition();
    P( M );

    P( "[ 0.64050615 -0.00557929  0.03081179  0.82131004 -0.48704869 ]" );
    P( "0.29491433" );
    P( M.solve_using_ldlt( X, 1 ) );

    // Vec<double> v( FromSizeAndItemValue(), m.nb_rows(), 1.0 );
    // Vec<double> s = m.solve( v );
    // P( s );
}
