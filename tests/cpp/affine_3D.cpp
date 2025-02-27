#include <tl/support/string/va_string.h>
#include "FindAffineTransformation.h"
#include "partial1D/VtkOutput.h"
#include <fstream>

// #undef P
// #include "CImg.h"
// using namespace cimg_library;

static constexpr int dim = 3;

using FA = FindAffineTransformation<dim>;
using Tr = AffineTransformation<dim>;
using TF = Tr::TF;
using Pt = Tr::Pt;

void load_points( Vec<Pt> &res, std::string filename ) {
    std::ifstream is( filename );
    while ( true ) {
        TF x, y, z;
        char t;
        is >> t >> x >> y >> z;
        if ( ! is )
            break;
        res << Pt{ x, y, z };
    }
}

// void disp( const CImg<unsigned char> &img ) {
//     CImgDisplay main_disp( img, "Click a point" );
//     while ( ! main_disp.is_closed() ) {
//         main_disp.wait();
//         // if ( ma)
//         if ( main_disp.button() | main_disp.key() ) {
//             // const int y = main_disp.mouse_y();
//             // visu.fill(0).draw_graph(image.get_crop(0,y,0,0,image.width()-1,y,0,0),red,1,1,0,255,0);
//             // visu.draw_graph(image.get_crop(0,y,0,1,image.width()-1,y,0,1),green,1,1,0,255,0);
//             // visu.draw_graph(image.get_crop(0,y,0,2,image.width()-1,y,0,2),blue,1,1,0,255,0).display(draw_disp);
//             break;
//         }
//     }
// }

void disp( FA f, PI num ) {
    VtkOutput vo_target;
    f.display_target_vtk( vo_target );
    vo_target.save( va_string( "results/mumble_target_$0.vtk", 1000 + num ) );

    VtkOutput vo_diracs;
    f.display_diracs_vtk( vo_diracs );
    vo_diracs.save( va_string( "results/mumble_diracs_$0.vtk", 1000 + num ) );

    VtkOutput vo_new_diracs;
    f.display_new_diracs_vtk( vo_new_diracs );
    vo_new_diracs.save( va_string( "results/mumble_new_diracs_$0.vtk", 1000 + num ) );

}

int main() {
    FA f;
    f.nb_bins = 400;
    f.nb_dirs = 10;

    // load points
    load_points( f.target, "../spot/Datasets/Pointsets/3D/mumble_sitting_100000.pts" );

    P( min( f.target ) );
    P( max( f.target ) );

    Vec<Pt> diracs;
    load_points( diracs, "../spot/Datasets/Pointsets/3D/mumble_sitting_80000.pts" );
    // load_points( diracs, "../spot/Datasets/Pointsets/3D/mumble_sitting_3000.pts" );
    for( PI n = 0; n < diracs.size(); ++n ) {
        Pt dir{ 0.42546746810417513, 0.35717031177892833, -0.8315087503861676 };
        Pt ori{ -26.23, 1.35, 50.83 };
        if ( sp( dir, diracs[ n ] ) < sp( dir, ori ) )
            f.diracs << diracs[ n ];
    }
    // for( PI n = 0; n < diracs.size(); ++n ) {
    //     Pt dir{ 0.932775485884615, -0.03429470970197802, -0.3588227498629058 };
    //     Pt ori{ -23.23316619310342, 0.76218257836255, 34.762438946432965 };
    //     if ( sp( dir, diracs[ n ] ) < sp( dir, ori ) )
    //         f.diracs << diracs[ n ];
    // }

    f.mass_ratio = f.diracs.size() / 80000.0;
    P( f.mass_ratio );
    // load_points( f.diracs, "../spot/Datasets/Pointsets/3D/mumble_sittingcut_80000.pts" );
    // f.mass_ratio = 0.8;

    Tr rot = Tr::rotation( 4 * M_PI / 180 ); // * Tr::translation( { -2, 1, 1 } );
    for( auto &d : f.diracs )
        d = rot.apply( d );

    // iterations
    for( PI num_iter = 0; num_iter < 1000; ++num_iter ) {
        // f.compute_with_rand_dir();
        f.compute_new_diracs( 5 );

        disp( f, num_iter );

        if ( f.iteration_SVD() < 1e-4 )
            break;
    }

    PE( f.dirac_trans );
}
