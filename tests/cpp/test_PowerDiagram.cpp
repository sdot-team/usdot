#include <partial1D/PowerDiagram.h>
#include <partial1D/Lebesgue.h>
#include <tl/support/ASSERT.h>
#include "catch_main.h"
using TF = FP64;

TEST_CASE( "PowerDiagram areas", "" ) {
    using Pd = PowerDiagram<FP64,Lebesgue<FP64>>;

    Vec<TF> positions, weights;
    for( PI i = 0; i < 10; ++i ) {
        const TF x = ( i + 0.5 ) / 10;
        positions << ( x - 0.5 ) * 10;
        weights << 1;
    }

    // weights[ 5 ] = 1.05;
    Pd pd( positions, weights, { 0, 1, /*h*/1 } );
    pd.update_weights();

    // P( pd.masses() );

    // P( pd.barycenters_ap() );
    // P( pd.barycenters() );

    // Vec<TF> xs, vs;
    // for( TF x = -2; x <= 3; x += 0.25/8 ) {
    //     vs << pd.density_value( x );
    //     xs << x;
    // }
    // P( xs, vs );
}

