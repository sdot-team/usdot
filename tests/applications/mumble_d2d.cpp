#include <partial1D/ShapeRegistration.h>

using namespace usdot;
using namespace std;

using TF = double;
using Sr = ShapeRegistration<TF,3>;

int main() {
    // shape    
    Sr sr;
    Sr::load_points( sr.shape_points, "../spot/Datasets/Pointsets/3D/mumble_sitting_100000.pts" );

    // diracs
    Sr::load_points( sr.diracs, "../spot/Datasets/Pointsets/3D/mumble_sitting_3000.pts" );
    TF min_x = min( sr.diracs )[ 0 ];
    TF max_x = max( sr.diracs )[ 0 ];
    // sr.diracs = sr.diracs.filtered( [&]( const auto &p ) {
    //     return p[ 0 ] < min_x + 0.6 * ( max_x - min_x );
    // } );

    // start from a false position
    sr.transformation = Sr::Tr::rotation( 15 * M_PI / 180 ) * Sr::Tr::translation( { 4, -4, 0 } );
    sr.mass_ratio = TF( sr.diracs.size() ) / 3000;

    //  
    sr.update_displacements();

    // ouptut
    VtkOutput vo;
    sr.display_shape_vtk( vo );
    sr.display_diracs_vtk( vo );
    sr.display_displacements_vtk( vo );
    vo.save( "results/mumble_sitting_diracs.vtk" );
}
