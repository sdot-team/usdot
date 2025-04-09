// test de procédure avec convolution décroissante
#include <partial1D/ConvolutedImage.h>
#include <partial1D/BoundedDensity.h>
#include <partial1D/PowerDiagram.h>

#include <tl/support/string/to_string.h>
#include <tl/support/operators/argmax.h>
#include <tl/support/operators/min.h>
#include <tl/support/operators/abs.h>
#include <eigen3/Eigen/Dense>
#include <fstream>

// using namespace boost::multiprecision;
using namespace std;

// using TF = number<backends::cpp_bin_float<40>>;
using TF = FP64;

/// assuming ws = value at 0, 1, ..., ws.size() - 1
Vec<TF> extrapolation( Vec<Vec<TF>> &ws, TF x, int n = -1 ) {
    if ( n < 0 )
        n = ws.size();

    using TM = Eigen::Matrix<TF,Eigen::Dynamic,Eigen::Dynamic>;
    using TV = Eigen::Matrix<TF,Eigen::Dynamic,1>;
    using std::pow;

    TM M( n, n );
    TV V( n );
    for( PI r = 0; r < n; ++r )
        for( PI c = 0; c < n; ++c )
            M.coeffRef( r, c ) = pow( TF( r ), c );

    Vec<TF> res( FromReservationSize(), ws.size() );
    Eigen::FullPivLU<TM> lu( M );
    for( PI i = 0; i < ws[ 0 ].size(); ++i ) {
        for( PI c = 0; c < n; ++c )
            V[ c ] = ws[ c ][ i ];
        auto R = lu.solve( V );

        TF v = 0;
        for( PI c = 0; c < n; ++c )
            v += pow( x, c ) * R[ c ];
        res << v;

    }

    return res;
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
    std::ofstream fs( "glot.py" );
    fs << "from matplotlib import pyplot\n";
    Vec<FP64> ys;
    for( PI nw = 0; nw < bp.size(); ++nw )
        ys << nw / TF( bp.size() - 1 );
    for( PI nd = 0; nd < 3/*bp[ 0 ].size()*/; ++nd ) {
        Vec<FP64> xs;
        for( PI nw = 0; nw < bp.size(); ++nw )
            xs << bp[ nw ][ nd ];
        fs << "pyplot.plot( " << to_string( xs ) << ", " << to_string( ys ) << " )\n";

        for( PI i = 0; i < bpi.size(); ++i ) {
            const auto &bp = bpi[ i ];
            if ( bp.empty() )
                continue;
            Vec<FP64> xs;
            for( PI nw = 0; nw < bp.size(); ++nw )
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
        fs << "pyplot.plot( " << to_string( ys ) << ", " << to_string( xbs ) << ", label = 'b" << o << "' )\n";
        // fs << "pyplot.plot( " << to_string( ys ) << ", " << to_string( xws ) << ", label = 'w" << o << "' )\n";
        P( o, argmax( xws > 1e-4 ) );
        P( o, argmax( xbs > 1e-3 ) );
    }
    fs << "pyplot.legend()\n";
    fs << "pyplot.show()\n";
}

struct System {
    void iteration() {
        // on fait une version localisée
        Vec<Vec<TF>> wsi;
        const TF epsilon = 1e-3;
        for( PI d = 0; d <= 6; ++d ) {
            c0->int_conv.width = base_width - epsilon * d;
            pd->update_weights();

            wsi << pd->sorted_seed_weights;
        }

        // 
        
    }
        
    RcPtr<ConvolutedImage<TF,IntInvX2<TF>>> c0;
    RcPtr<BoundedDensity<TF>> cb;
    RcPtr<PowerDiagram<TF>> pd;
    TF base_width;
};


// using TF = cpp_bin_float_100;
int main() {
    System system;

    PI nb_cells = 20;
    Vec<TF> positions, weights;
    for( PI i = 0; i < nb_cells; ++i ) {
        const TF x = TF( i + 0.5 ) / nb_cells;
        positions << x * 1;
        weights << 1;
    }

    // auto cb = RcPtr<ConvolutedLebesgue<TF,IntGauss<TF>>>::New( -4, -2, /*h*/1, IntGauss<TF>{ 1e-3 } );
    system.c0 = RcPtr<ConvolutedImage<TF,IntInvX2<TF>>>::New( Vec<TF>{ 1, 0, 1 } );
    system.cb = RcPtr<BoundedDensity<TF>>::New( -3, +6, system.c0 );
    system.pd = RcPtr<PowerDiagram<TF>>::New( positions, weights, system.cb );
    system.pd->target_inf_mass_error = numeric_limits<TF>::epsilon() * 1e4;
    system.base_width = 0.5;



    //
    Vec<Vec<Vec<TF>>> bpe( FromSize(), 10 ); // barycentres avec poids extrapolés
    Vec<Vec<Vec<TF>>> wse( FromSize(), 10 ); // poids extrapolés
    Vec<Vec<FP64>> bp;
    Vec<Vec<TF>> ws;
    TF wo = base_width;
    auto interp_width = [&]( TF w1, TF nd = 100 ) {
        for( PI d = 0; d <= nd; ++d ) {
            const TF width = wo + ( w1 - wo ) * pow( d / nd, 1 );
            c0->int_conv.width = width;

            P( width );
            bool res = pd.update_weights();
            ASSERT( res );
    
            bp << pd.barycenters_ap( 1000 );
            ws << pd.sorted_seed_weights;

            // on sort une version extrapolée
            for( PI o = 2; o <= 5; ++o ) {
                pd.sorted_seed_weights = extrapolation( wsi, ( base_width - width ) / epsilon, o + 1 );
                bpe[ o ] << pd.barycenters_ap( 10000 );
                wse[ o ] << pd.sorted_seed_weights;
                pd.sorted_seed_weights = ws.back();
            }
        }
        wo = w1;
    };
    interp_width( 1.0e-3, 80 );
    // interp_width( 1e-3, 100 );

    // glot_barycenters( bp, bpe );
    glot_erreurs( bp, bpe, ws, wse );

    // for( PI d = 0, v = 20; d <= v; ++d ) {
    //     pd.sorted_seed_weights = extrapolation( ws, d * nd / v / epsilon );
    //     bp << pd.cell_boundaries();
    // }

    // // P( inside_boundaries( bp.back() ) );

    // // pd.density.std = 1e-2;
    // // pd.update_weights();
    // // P( inside_boundaries( pd.cell_boundaries() ) );

}

