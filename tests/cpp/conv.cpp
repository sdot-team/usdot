// #include "boost/multiprecision/fwd.hpp"
// #include <partial1D/ConvolutedLebesgue.h>
// #include <partial1D/ConvolutedImage.h>
#include <partial1D/ConvolutedDiracs.h>
#include <partial1D/SumOfDensities.h>
#include <partial1D/BoundedDensity.h>
#include <partial1D/PowerDiagram.h>
#include <partial1D/Lebesgue.h>

#include <tl/support/string/to_string.h>
#include <tl/support/operators/argmax.h>
#include <tl/support/operators/min.h>
#include <tl/support/operators/abs.h>
#include <eigen3/Eigen/Dense>
#include <fstream>

using namespace std;

// #include <partial1D/gmp.h>
// using namespace boost::multiprecision;
// using TF = number<backends::cpp_bin_float<40>>;

using TF = FP64;


Vec<FP64> inside_boundaries( const Vec<FP64> &bp ) {
    Vec<FP64> lp;
    for( PI n = 1; n + 1 < bp.size(); n += 2 )
        lp << bp[ n ];
    return lp;
}

void glot( TF mi, TF ma, PI nd, auto &&...funcs ) {
    std::ofstream fs( "glot.py" );
    fs << "from matplotlib import pyplot\n";
    auto pf = [&]( auto &&func ) {
        Vec<FP64> xs, ys;
        for( PI n = 0; n <= nd; ++n ) {
            TF x = mi + ( ma - mi ) * n / nd;
            ys << func( x );
            xs << x;
            fs << "pyplot.plot( " << to_string( xs ) << ", " << to_string( ys ) << " )\n";
        }
    };
    ( pf( funcs ), ... );
    fs << "pyplot.show()\n";
}

void glot_barycenters( const Vec<Vec<FP64>> &bp, const Vec<Vec<Vec<FP64>>> &bpi ) {
    const PI beg = bp.size() * 0 / 2;
    const PI end = bp.size() * 2 / 2;

    std::ofstream fs( "glot.py" );
    fs << "from matplotlib import pyplot\n";
    Vec<FP64> ys;
    for( PI nw = beg; nw < end; ++nw )
        ys << ( nw - beg ) / TF( end - beg - 1 );
    for( PI nd = 0; nd < bp[ 0 ].size(); ++nd ) {
        Vec<FP64> xs;
        for( PI nw = beg; nw < end; ++nw )
            xs << bp[ nw ][ nd ];
        fs << "pyplot.plot( " << to_string( xs ) << ", " << to_string( ys ) << " )\n";

        for( PI i = 0; i < bpi.size(); ++i ) {
            const auto &bp = bpi[ i ];
            if ( bp.empty() )
                continue;
            Vec<FP64> xs;
            for( PI nw = beg; nw < end; ++nw )
                xs << bp[ nw ][ nd ];
            Vec<Str> colors{ "red", "blue", "grey", "green", "orange" };
            if ( nd == 0 )
                fs << "pyplot.plot( " << to_string( xs ) << ", " << to_string( ys ) << ", '--', label = '" << i << "', color = '" << colors[ i - 2 ] << "' )\n";
            else
                fs << "pyplot.plot( " << to_string( xs ) << ", " << to_string( ys ) << ", '--', color = '" << colors[ i - 2 ] << "' )\n";
        }
    }
    fs << "pyplot.legend()\n";
    fs << "pyplot.show()\n";
}

void glot_erreurs( const Vec<Vec<FP64>> &bp, const Vec<Vec<Vec<FP64>>> &bpe, const Vec<Vec<FP64>> &ws, const Vec<Vec<Vec<FP64>>> &wse ) {
    std::ofstream fs( "glot.py" );
    fs << "from matplotlib import pyplot\n";
    Vec<FP64> ys;
    for( PI ne = 0; ne < bp.size(); ++ne )
        ys << ne / TF( bp.size() - 1 );

    for( PI o = 0; o < bpe.size(); ++o ) {
        if ( bpe[ o ].empty() )
            continue;

        Vec<FP64> xbs;
        Vec<FP64> xws;
        for( PI ne = 0; ne < bp.size(); ++ne ) {
            FP64 xb = 0;
            FP64 xw = 0;
            for( PI nd = 0; nd < bp[ 0 ].size(); ++nd ) {
                xb = max( xb, abs( bp[ ne ][ nd ] - bpe[ o ][ ne ][ nd ] ) );
                xw = max( xw, abs( ws[ ne ][ nd ] - wse[ o ][ ne ][ nd ] ) );
            }
            xbs << xb;
            xws << xw;
        }
        // fs << "pyplot.plot( " << to_string( ys ) << ", " << to_string( xbs ) << ", label = 'b" << o << "' )\n";
        fs << "pyplot.plot( " << to_string( ys ) << ", " << to_string( xws ) << ", label = 'w" << o << "' )\n";
        P( o, argmax( xws > 1e-4 ) );
        P( o, argmax( xbs > 1e-3 ) );
    }
    fs << "pyplot.legend()\n";
    fs << "pyplot.show()\n";
}

// using TF = cpp_bin_float_100;
int main() {
    using Dp = RcPtr<Density<TF>>;
    using Pd = PowerDiagram<TF>;

    // auto cc = RcPtr<ConvolutedLebesgue<TF,IntGauss<TF>>>::New( 0, 1, /*h*/1, IntGauss<TF>{ 0.2 } );
    // auto c0 = RcPtr<ConvolutedLebesgue<TF,IntInvX2<TF>>>::New( -2, -1, /*h*/1, IntInvX2<TF>{ 0.5 } );
    // auto c1 = RcPtr<ConvolutedLebesgue<TF,IntInvX2<TF>>>::New( +1, +2, /*h*/1, IntInvX2<TF>{ 0.5 } );
    // auto c0 = RcPtr<Lebesgue<TF>>::New( -4, -1 );
    // auto c1 = RcPtr<Lebesgue<TF>>::New( +1, +4 );
    // auto c0 = RcPtr<ConvolutedLebesgue<TF,IntInvX2<TF>>>::New( -4, -2, /*h*/1, IntInvX2<TF>{ 1e-3 } );
    // auto c1 = RcPtr<ConvolutedLebesgue<TF,IntInvX2<TF>>>::New( +2, +4, /*h*/1, IntInvX2<TF>{ 1e-3 } );
    // auto c0 = RcPtr<ConvolutedLebesgue<TF,IntGauss<TF>>>::New( -4, -2, /*h*/1, IntGauss<TF>{ 1e-3 } );
    // auto c1 = RcPtr<ConvolutedLebesgue<TF,IntGauss<TF>>>::New( +2, +4, /*h*/1, IntGauss<TF>{ 1e-3 } );
    // const TF base_width = 2.0;
    // auto cb = RcPtr<SumOfDensities<TF>>::New( Vec<Dp>{ c0, c1 } );
    // // auto cb = RcPtr<BoundedDensity<TF>>::New( -10, +10, cs );

    // auto cb = RcPtr<ConvolutedLebesgue<TF,IntGauss<TF>>>::New( -4, -2, /*h*/1, IntGauss<TF>{ 1e-3 } );
    // auto c0 = RcPtr<ConvolutedImage<TF,IntInvX2<TF>>>::New( Vec<TF>{ 1, 0, 1 } );
    auto c0 = RcPtr<ConvolutedDiracs<TF,IntGauss<TF>>>::New( Vec<TF>{ 0, 1 }, Vec<TF>{ 1, 1 }, IntGauss<TF>{ 1e-1 } );
    auto cb = RcPtr<BoundedDensity<TF>>::New( -3, +2, c0 );
    const TF base_width = 0.5;

    // glot( -3, +6, 500, 
    //     [&]( TF x ) { return cb->value( x ); },
    //     [&]( TF x ) { return cb->integral( x, x + 1e-3 ) / 1e-3; }
    // );

    PI nb_cells = 4;
    Vec<TF> positions, weights;
    for( PI i = 0; i < nb_cells; ++i ) {
        const TF x = TF( i + 0.5 ) / nb_cells;
        positions << x * 1;
        weights << 1;
    }

    Pd pd( positions, weights, cb );
    pd.target_inf_mass_error = numeric_limits<TF>::epsilon() * 1e4;
    pd.relax = 1e-2;

    const TF intermediate_width = 0.9 * base_width; 

    //
    TF wo = base_width;
    auto interp_width_0 = [&]( TF w1, TF nd = 100 ) {
        for( PI d = 0; d <= nd; ++d ) {
            const TF width = wo + ( w1 - wo ) * pow( d / nd, 1 );
            c0->int_conv.width = width;
            pd.update_weights();
        }
        wo = w1;
    };
    interp_width_0( intermediate_width, 40 );
    
    // on fait un version localisée pour trouver les dérivées
    Vec<Vec<TF>> wsi;
    const TF epsilon = 1e-3;
    for( PI d = 0; d <= 6; ++d ) {
        c0->int_conv.width = intermediate_width - epsilon * d;
        pd.update_weights();
        wsi << pd.sorted_seed_weights;
    }

    // on va jusqu'au bout
    Vec<Vec<Vec<TF>>> bpe( FromSize(), 10 ); // barycentres avec poids extrapolés
    Vec<Vec<Vec<TF>>> wse( FromSize(), 10 ); // poids extrapolés
    Vec<Vec<FP64>> bp;
    Vec<Vec<TF>> ws;
    auto interp_width_1 = [&]( TF w1, TF nd = 100 ) {
        for( PI d = 0; d <= nd; ++d ) {
            const TF width = wo + ( w1 - wo ) * pow( d / nd, 1 );
            c0->int_conv.width = width;

            P( width );
            pd.update_weights();
    
            bp << pd.barycenters_ap( 1000 );
            ws << pd.sorted_seed_weights;

            // on sort une version extrapolée
            // for( PI o = 2; o <= 5; ++o ) {
            //     pd.sorted_seed_weights = extrapolation( wsi, ( intermediate_width - width ) / epsilon, o + 1 );
            //     bpe[ o ] << pd.barycenters_ap( 10000 );
            //     wse[ o ] << pd.sorted_seed_weights;
            //     pd.sorted_seed_weights = ws.back();
            // }
        }
        wo = w1;
    };
    // for( TF x = 0.5; x > 1e-4; x *= 0.9 )        
    interp_width_1( 1e-3, 20000 );

    glot_barycenters( bp, bpe );
    // glot_erreurs( bp, bpe, ws, wse );

    // for( PI d = 0, v = 20; d <= v; ++d ) {
    //     pd.sorted_seed_weights = extrapolation( ws, d * nd / v / epsilon );
    //     bp << pd.cell_boundaries();
    // }

    // P( inside_boundaries( bp.back() ) );

    // pd.density.std = 1e-2;
    // pd.update_weights();
    // P( inside_boundaries( pd.cell_boundaries() ) );
}

