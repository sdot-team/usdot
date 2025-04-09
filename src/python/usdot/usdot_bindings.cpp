// cppimport
#include "../../cpp/usdot/System.h"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace usdot;
using TF = double;

using Array = pybind11::array_t<TF, pybind11::array::c_style>;

struct OtResult {
    Array norm_2_residual_history;
    Str   error_message;
    Array barycenters;
    Array boundaries;
    Array weights;
    Array masses;
};

struct OtParms {
    TF max_mass_ratio_error_target = 1e-6;
    int verbosity = 1;
    TF epsilon = 0;
};

static Array array_from_vec( const std::vector<TF> &v ) {
    Array res( std::array<PI,1>{ v.size() } );
    for( PI i = 0; i < v.size(); ++i )
        res.mutable_at( i ) = v[ i ];
    return res;
}

template<PI N>
static Array array_from_vec( const std::vector<std::array<TF,N>> &v ) {
    Array res( std::array<PI,2>{ v.size(), PI( N ) } );
    for( PI i = 0; i < v.size(); ++i )
        for( PI d = 0; d < N; ++d )
            res.mutable_at( i, d ) = v[ i ][ d ];
    return res;
}

static std::vector<TF> vec_from_array( const Array &array ) {
    std::vector<TF> res( array.size() );
    for( PI r = 0; r < res.size(); ++r )
        res[ r ] = array.at( r );
    return res;
}

template<std::size_t N> 
static std::vector<std::array<TF,N>> vec_of_vec_from_array( const Array &array ) {
    std::vector<std::array<TF,N>> res( array.shape( 0 ) );
    for( PI r = 0; r < res.size(); ++r )
        for( PI c = 0; c < N; ++c )
            res[ r ][ c ] = array.at( r, c );
    return res;
}

// OtResult ot_solve( Array &dirac_positions, TF global_mass_ratio, Array &relative_mass_ratios, RcPtr<Density<TF>> de, OtParms &parms ) {
//     // init power diagram
//     auto pd = RcPtr<PowerDiagram<TF>>::New( vec_from_array( dirac_positions ), Vec<TF>::fill( dirac_positions.size(), 1 ) );
//     pd->epsilon = parms.epsilon;
 
//     // init solver
//     Solver<TF> solver( pd, de );
//     solver.max_mass_ratio_error_target = parms.max_mass_ratio_error_target;
//     solver.relative_mass_ratios = vec_from_array( relative_mass_ratios );
//     solver.global_mass_ratio = global_mass_ratio;
//     solver.verbosity = parms.verbosity;

//     // run
//     solver.solve();

//     // result summary
//     Vec<Vec<TF,2>> bnds = pd->cell_boundaries( false );
//     TF mi = de->min_x();
//     TF ma = de->max_x();
//     for( Vec<TF,2> &b : bnds )
//         for( TF &v : b )
//             v = min( ma, max( mi, v ) );

//     OtResult res;
//     res.norm_2_residual_history = array_from_vec( solver.norm_2_residual_history );
//     res.error_message = ""; // solver.nb_errors ? "error" : 
//     res.barycenters = array_from_vec( pd->barycenters( *de, false ) );
//     res.boundaries = array_from_vec( bnds );
//     res.weights = array_from_vec( pd->get_weights( false ) );
//     res.masses = array_from_vec( pd->masses( *de, false ) );
//     return res;
// }

// OtResult ot_diracs_to_piecewise_polynomial( Array &dirac_positions, TF global_mass_ratio, Array &target_mass_ratios, Array &density_positions, Array &density_values, OtParms &parms ) {
//     PI nc = density_values.ndim() >= 2 ? density_values.shape( 1 ) : 1;
//     RcPtr<Density<TF>> de;
//     if ( nc == 1 ) {
//         de = RcPtr<PiecewiseConstantDensity<TF>>::New( 
//             vec_from_array( density_positions ),
//             vec_from_array( density_values )
//         );
//     } else if ( nc == 2 ) {
//         de = RcPtr<PiecewiseAffineDensity<TF>>::New( 
//             vec_from_array( density_positions ),
//             vec_of_vec_from_array<2>( density_values )
//         );
//     } else {
//         throw std::runtime_error( "TODO: high order polynomials" );
//     }

//     return ot_solve( dirac_positions, global_mass_ratio, target_mass_ratios, de, parms );
// }

OtResult from_p1_grid( Array &dirac_positions, TF global_mass_ratio, Array &relative_mass_ratios, Array density_values, TF density_beg, TF density_end, OtParms &parms ) {
    // normalize
    std::vector<TF> dp = vec_from_array( dirac_positions );
    std::vector<TF> dv = vec_from_array( density_values );
    TF pos_ratio = 1;
    if ( density_beg != 0 || density_end != dv.size() - 1 ) {
        pos_ratio = ( dv.size() - 1 ) / ( density_end - density_beg );
        for( TF &v : dp )
            v = ( v - density_beg ) * pos_ratio;
    }

    // prepare the system
    GridDensity<TF> gd( std::move( dv ) );

    System<TF,GridDensity<TF>> si;
    si.verbosity = parms.verbosity;
    si.stream = &std::cout;

    si.set_relative_dirac_masses( vec_from_array( relative_mass_ratios ) );
    si.set_global_mass_ratio( global_mass_ratio );
    si.set_dirac_positions( dp );
    si.set_density( &gd );

    // solve it
    si.solve();

    // make a summary
    auto barycenters = si.cell_barycenters();
    auto boundaries = si.cell_boundaries();
    auto weights = si.dirac_weights();
    auto masses = si.cell_masses();

    if ( density_beg != 0 || density_end != dv.size() - 1 ) {
        for( TF &v : barycenters )
            v = density_beg + v / pos_ratio;

        for( std::array<TF,2> &b : boundaries )
            for( TF &v : b )
                v = density_beg + v / pos_ratio;

        for( TF &v : weights )
            v /= pow( pos_ratio, 2 );
    
        for( TF &v : masses )
            v /= pos_ratio;
    }

    // make a summary
    OtResult res;
    res.norm_2_residual_history = array_from_vec( {} ); // solver.norm_2_residual_history
    res.error_message = ""; // solver.nb_errors ? "error" : 
    res.barycenters = array_from_vec( barycenters );
    res.boundaries = array_from_vec( boundaries );
    res.weights = array_from_vec( weights );
    res.masses = array_from_vec( masses );
    return res;
}

PYBIND11_MODULE( usdot_bindings, m ) {
    pybind11::class_<OtResult>( m, "OtResult" )
        .def_readwrite( "norm_2_residual_history", &OtResult::norm_2_residual_history, "" )
        .def_readwrite( "error_message", &OtResult::error_message, "" )
        .def_readwrite( "barycenters", &OtResult::barycenters, "" )
        .def_readwrite( "boundaries", &OtResult::boundaries, "" )
        .def_readwrite( "weights", &OtResult::weights, "" )
        .def_readwrite( "masses", &OtResult::masses, "" )
        .def( "__repr__",
            []( const OtResult &a ) {
                return "{ error_message: '" + a.error_message + "', barycenters: [...], boundaries: [...], weights: [...] }";
            }
        )
        ;

    pybind11::class_<OtParms>( m, "OtParms" )
        .def( py::init<>() )
        .def_readwrite( "max_mass_ratio_error_target" , &OtParms::max_mass_ratio_error_target, "" )
        .def_readwrite( "verbosity" , &OtParms::verbosity, "" )
        .def_readwrite( "epsilon" , &OtParms::epsilon, "" )
        .def( "__repr__",
            []( const OtParms &a ) {
                return "{ max_mass_ratio_error_target: " + std::to_string( a.max_mass_ratio_error_target ) + " }";
            }
        )
        ;;

    m.def( "from_p1_grid", &from_p1_grid );
} 
/* 
<%
setup_pybind11(cfg)
cfg['include_dirs'] += [ '../../../ext' ]
cfg['extra_compile_args'] = ['-std=c++20']

# cfg['sources'] = ['extra_source1.cpp', 'extra_source2.cpp']

# import logging
# logging.basicConfig()

cfg['dependencies'] = [
    '../../cpp/usdot/WeightInitializer.cxx',
    '../../cpp/usdot/WeightInitializer.h',
    '../../cpp/usdot/WeightUpdater.cxx',
    '../../cpp/usdot/WeightUpdater.h',
    '../../cpp/usdot/System.cxx',
    '../../cpp/usdot/System.h',
]
%>
*/