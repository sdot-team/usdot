#include <partial1D/Density/PiecewiseConstantDensity.h>
#include <partial1D/Density/ConvolutedDiracs.h>

#include "partial1D/Convolution/InvX2Convolution.h"

#include <partial1D/Solver.h>
#include <fstream>

#include "catch_main.h"

using namespace usdot;
using namespace std;
using TF = double;

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

TEST_CASE( "Density", "" ) {
    auto cd = RcPtr<ConvolutedDiracs<TF>>::New( Vec<TF>{ 0, 1 }, Vec<TF>{ 1, 1 }, new GaussianConvolution<TF>( 0.1 ) );
    auto cb = RcPtr<PiecewiseConstantDensity<TF>>::New( Vec<TF>::linspace( 0.0, 1.0, 4 ), Vec<TF>{ 1, 0, 1 } );
    auto cv = cb->convoluted( new GaussianConvolution<TF>( 0.1 ) );

    // for( TF x : Vec<TF>::linspace( 1, 2, 200 ) ) {
    //     CHECK_PROX( cd->value( x ), cd->integral( x, x + 1e-6 ) / 1e-6, 1e-4 );
    //     CHECK_PROX( cv->value( x ), cv->integral( x, x + 1e-6 ) / 1e-6, 1e-4 );
    // }

    glot( Vec<TF>::linspace( -1, 2, 200 ),
        // [&]( TF x ) { return 1 + cd->value( x ); }, [&]( TF x ) { return 1 + cd->integral( x, x + 1e-2 ) / 1e-2; }
        // [&]( TF x ) { return 0 + cb->value( x ); }, [&]( TF x ) { return 0 + cb->integral( x, x + 1e-4 ) / 1e-4; }
        [&]( TF x ) { return 2 + cv->value( x ); }, [&]( TF x ) { return 2 + cv->integral( x, x + 1e-2 ) / 1e-2; }
    );
}
