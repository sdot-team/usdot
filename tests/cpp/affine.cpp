#include <tl/support/operators/norm_2.h>
#include "FindAffineTransformation.h"

#undef P
#include "CImg.h"
using namespace cimg_library;

static constexpr int dim = 2;

using FA = FindAffineTransformation<dim>;
using Tr = AffineTransformation<dim>;
using TF = Tr::TF;
using Pt = Tr::Pt;

void add_disc( Vec<Pt> &pts, Pt C, TF R, PI nb_points ) {
    for( PI n = 0; n < nb_points; ) {
        Pt p{
            ( 2.0 * rand() / RAND_MAX - 1.0 ),
            ( 2.0 * rand() / RAND_MAX - 1.0 ),
        };
        if ( norm_2_p2( p ) <= 1 ) {
            pts << C + R * p;
            ++n;
        }
    }
}

int main() {
    FA f;
    f.nb_dirs = 20;

    Vec<std::pair<Pt,TF>> discs{
        { Pt{ 0.0, 0.0 }, 1.0 },
        { Pt{ 1.5, 0.0 }, 0.5 },
    }; 

    // target
    for( auto &disc : discs )
        add_disc( f.target, disc.first, disc.second, 10000 * pow( disc.second, 2 ) );

    // diracs
    Tr trans = Tr::rotation( 30 * M_PI / 180 ); // * Tr::translation( { 1, 0 } );
    // Tr trans = Tr::translation( { 1, 0 } );
    for( auto &disc : discs )
        add_disc( f.diracs, trans.apply( disc.first ), disc.second, 1000 * pow( disc.second, 2 ) );

    // iterations
    for( PI num_iter = 0; num_iter < 20; ++num_iter ) {
        f.compute_new_diracs( 10 );

        //
        CImg<unsigned char> img( 800, 800, 1, 3, 0 );
        unsigned char white[] = { 255, 255, 255 };
        img.draw_rectangle( 0, 0, img.width(), img.height(), white );
        f.display_points( img, Tr::scale( 200 ) * Tr::translation( { 300, 400 } ) );

        std::string fn = std::format( "results/2discs_{}.png", 100 + num_iter );
        img.save_png( fn.c_str() );

        // CImgDisplay main_disp( img, "Click a point" );
        // while ( ! main_disp.is_closed() ) {
        //     main_disp.wait();
        //     // if ( ma)
        //     if ( main_disp.button() | main_disp.key() ) {
        //         // const int y = main_disp.mouse_y();
        //         // visu.fill(0).draw_graph(image.get_crop(0,y,0,0,image.width()-1,y,0,0),red,1,1,0,255,0);
        //         // visu.draw_graph(image.get_crop(0,y,0,1,image.width()-1,y,0,1),green,1,1,0,255,0);
        //         // visu.draw_graph(image.get_crop(0,y,0,2,image.width()-1,y,0,2),blue,1,1,0,255,0).display(draw_disp);
        //         break;
        //     }
        // }

        f.iteration_SVD();
    }
}
