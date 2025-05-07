#include <tl/support/P.h>

#include <usdot/ConvolutedDensity.h>
#include <usdot/utility/glot.h>
#include "catch_main.h"
// #include "usdot/utility/linspace.h"

using namespace usdot;
using namespace std;
using TF = double;

TEST_CASE( "ConvolutedDensity", "" ) {
    std::vector<double> values( 10, 0 );
    for( PI i = 0; i < values.size(); ++i )
        values[ i ] = ( i > values.size() / 2 );

    ConvolutedDensity<TF,2> gd( values );
    

    glot( linspace<TF>( -1.0, 11.0, 100 ), 
        [&]( TF x ) { return gd.value( x, 3 ); },
        [&]( TF x ) { return ( gd.value( x, 2 ) - ge.value( x, 2 ) ) / eps; }
    );
}
