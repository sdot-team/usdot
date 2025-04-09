#include <partial1D/BandMatrix.h>
#include "catch_main.h"

TEST_CASE( "BandMatrix", "" ) {
    for( PI i = 0; i < 5; ++i ) {
        BandMatrix<double> m( 5 );
        for( PI i = 0; i < m.data.size(); ++i )
            m.data[ i ] = double( rand() ) / RAND_MAX + ( i % 2 == 0 );
        // P( m );

        Vec<double> v( FromSizeAndItemValue(), m.nb_rows(), 1.0 );
        Vec<double> s = m.solve( v );
        // P( s );

        P( m.mul( s ) );
    }
}
