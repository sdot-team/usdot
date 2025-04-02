#include <partial1D/ConvGridSolver.h>
#include <partial1D/glot.h>
#include "catch_main.h"

// #include <partial1D/gmp.h>
// using namespace boost::multiprecision;
// using TF = number<backends::cpp_bin_float<40>>;

using namespace usdot;
using namespace std;
using TF = double;

// struct System {
//     using TM = Eigen::Matrix<TF,Eigen::Dynamic,Eigen::Dynamic>;
//     using TV = Eigen::Matrix<TF,Eigen::Dynamic,1>;

//     System() {
//         dirac_positions = Vec<TF>::cellspace( 0.0, 0.5, 5 );
//         dirac_masses = Vec<TF>::fill( dirac_positions.size(), TF( 1 ) / dirac_positions.size() );
//         dirac_weights = Vec<TF>::fill( dirac_positions.size(), pow( TF( 1 ) / dirac_positions.size(), 2 ) / 2 );
//     }

//     PI nb_cells() const {
//         return dirac_positions.size();
//     }

//     TF _cut( PI n1 ) const {
//         const TF d0 = dirac_positions[ n1 - 1 ];
//         const TF d1 = dirac_positions[ n1 - 0 ];
//         const TF w0 = dirac_weights[ n1 - 1 ];
//         const TF w1 = dirac_weights[ n1 - 0 ];
//         return ( d0 + d1 + ( w0 - w1 ) / ( d1 - d0 ) ) / 2;
//     }

//     TF c0( PI n ) const {
//         return n ? _cut( n ) : numeric_limits<TF>::lowest();
//     }

//     TF c1( PI n ) const {
//         return n + 1 < dirac_positions.size() ? _cut( n + 1 ) : numeric_limits<TF>::max();
//     }

//     TF error() {
//         TF err = 0;
//         for( TF v : d0_error() )
//             err += pow( v, 2 );
//         return err;
//     }

//     Vec<Vec<TF,2>> cell_boundaries() const {
//         Vec<Vec<TF,2>> res( FromSize(), nb_cells() );
//         for( PI num_cell = 0; num_cell < dirac_positions.size(); ++num_cell ) {
//             TF rd = sqrt( dirac_weights[ num_cell ] );
//             TF pd = dirac_positions[ num_cell ];
 
//             TF p0 = max( c0( num_cell ), pd - rd );
//             TF p1 = min( c1( num_cell ), pd + rd );

//             res[ num_cell ][ 0 ] = p0;
//             res[ num_cell ][ 1 ] = p1;
//         }
//         return res;
//     }
    
//     Vec<TF> d0_error( TF mass_ratio = 1 ) {
//         Vec<TF> res( FromSize(), nb_cells() );
//         for( PI num_cell = 0; num_cell < dirac_positions.size(); ++num_cell ) {
//             TF rd = sqrt( dirac_weights[ num_cell ] );
//             TF pd = dirac_positions[ num_cell ];
 
//             TF p0 = max( c0( num_cell ), pd - rd );
//             TF p1 = min( c1( num_cell ), pd + rd );

//             res[ num_cell ] = mass_ratio * dirac_masses[ num_cell ] - ( p1 - p0 );
//         }
//         return res;
//     }

//     TM d1_error( TF eps = 1e-6 ) {
//         const Vec<TF> ref_d0_error = d0_error();
//         TM res( nb_cells(), nb_cells() );
//         for( PI n = 0; n < nb_cells(); ++n ) {
//             TF &ref_weight = dirac_weights[ n ];
//             TF old_weight = ref_weight;
//             ref_weight += eps;
    
//             const Vec<TF> new_d0_error = d0_error();
//             for( PI m = 0; m < nb_cells(); ++m )
//                 res.coeffRef( m, n ) = ( new_d0_error[ m ] - ref_d0_error[ m ] ) / eps;
    
//             ref_weight = old_weight;
//         }

//         const PI e = nb_cells() - 1;
//         res.coeffRef( 0, 0 ) += 1e10;
//         res.coeffRef( e, e ) += 1e10;

//         return res;
//     }

//     TV to_TV( auto &&v ) {
//         TV res( v.size() );
//         for( PI m = 0; m < v.size(); ++m )
//             res[ m ] = v[ m ];
//         return res;
//     }

//     void solve() {
//         // solve 
//         TF c0 = dirac_positions.front() - sqrt( dirac_weights.front() );
//         TF c1 = dirac_positions.back() + sqrt( dirac_weights.back() );
//         TF cm = c1 - c0;
//         P( cm );

//         TM M = d1_error();
//         Eigen::FullPivLU<TM> lu( M );
//         Vec<TF> X( lu.solve( to_TV( d0_error( cm ) ) ) );
        
//         P( X );
//         dirac_weights -= X;
//     }

//     Vec<TF> dirac_positions;
//     Vec<TF> dirac_weights;
//     Vec<TF> dirac_masses;
// };

TEST_CASE( "Conv Grid solver", "" ) {
    // System system;
    // P( system.error() );
    // system.solve();
    // P( system.error() );

    // for( auto v : system.cell_boundaries() )
    //     P( v[ 1 ] - v[ 0 ] );

    ConvGridSolverInput<TF> si;
    si.dirac_positions = Vec<TF>::cellspace( 0, 1, 8 );
    si.starting_filter_value = 0.0; //99;
    si.target_filter_value = 0.0;
    si.global_mass_ratio = 2.0 / 3.0;

    si.density_values = { 1, 1 };
    si.beg_x_density = -1;
    si.end_x_density = 2;
 
    ConvGridSolver<TF> solver( std::move( si ) );
    // P( solver.normalized_cell_boundaries() );
    P( solver.sorted_dirac_weights );
    solver.solve();
    solver.plot();

    // glot( Vec<TF>::linspace( -7, 12, 1000 ), 
    //     [&]( TF x ) { return solver.density_value( x ); }
    // );
}
// TF cp( TF mi ) {
//     ConvGridSolverInput<TF> si;
//     si.dirac_positions = { 0, 0.5, 1 };
//     si.starting_filter_value = 1e-2;
//     si.target_filter_value = 0.0;

//     si.density_values = { mi, mi, 1, 1 };
//     si.beg_x_density = 0;
//     si.end_x_density = 1;
 
//     ConvGridSolver<TF> solver( std::move( si ) );
//     auto d = solver.newton_dir().first;
//     for( PI i = 0; i < solver.nb_diracs(); ++i )
//         solver.sorted_dirac_weights[ i ] += d[ i ];

//     auto bnds = solver.normalized_cell_boundaries();
//     P( d, bnds );

//     return bnds[ 0 ][ 1 ];
// }

// TEST_CASE( "Conv Grid solver", "" ) {
//     for( TF mi = 1; mi > 1e-3; mi /= 2 )
//         P( mi, cp( mi ) );
// }

