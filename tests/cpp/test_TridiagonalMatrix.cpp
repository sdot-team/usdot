#include <partial1D/TridiagonalMatrix.h>
#include "catch_main.h"

TEST_CASE( "TridiagonalMatrix", "" ) {
    TridiagonalMatrix<double> m( 5 );
    for( PI i = 0; i < m.data.size(); ++i )
        m.data[ i ] = i;
    P( m );
    //     // P( m );

    //     Vec<double> v( FromSizeAndItemValue(), m.nb_rows(), 1.0 );
    //     Vec<double> s = m.solve( v );
    //     // P( s );

    //     P( m.mul( s ) );
    // }
}
