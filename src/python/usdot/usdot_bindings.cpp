// cppimport
#include "../../cpp/partial1D/Density/PiecewiseConstantDensity.h"
#include "../../cpp/partial1D/Solver.h"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <tl/support/P.h>

#include <tl/support/display/DisplayItem_Pointer.cpp>
#include <tl/support/display/DisplayItem_String.cpp>
#include <tl/support/display/DisplayItem_Number.cpp>
#include <tl/support/display/DisplayItem_List.cpp>
#include <tl/support/display/DisplayItem.cpp>
#include <tl/support/Displayer.cpp>

#include <tl/support/string/read_arg_name.cpp>

namespace py = pybind11;
using namespace usdot;
using TF = double;

using Array = pybind11::array_t<TF, pybind11::array::c_style>;

struct OtResult {
    Array barycenters;
    Array boundaries;
    Array weights;
    std::string error_message;
};

struct OtParms {
    TF max_mass_ratio_error_target = 1e-6;
};

static Array array_from_vec( const Vec<TF> &v ) {
    Array res( std::array<PI,1>{ v.size() } );
    for( PI i = 0; i < v.size(); ++i )
        res.mutable_at( i ) = v[ i ];
    return res;
}

static Vec<TF> vec_from_array( const Array &array ) {
    Vec<TF> res( FromSize(), array.size() );
    for( PI r = 0; r < res.size(); ++r )
           res[ r ] = array.at( r );
    return res;
}

OtResult ot_solve( Array &dirac_positions, Array &target_mass_ratios, RcPtr<Density<TF>> de, OtParms &parms ) {
    auto pd = RcPtr<PowerDiagram<TF>>::New(
        vec_from_array( dirac_positions ),
        Vec<TF>::fill( dirac_positions.size(), 1 )
    );
 
    Solver<TF> solver( pd, de );
    solver.max_mass_ratio_error_target = parms.max_mass_ratio_error_target;
    solver.target_mass_ratios = vec_from_array( target_mass_ratios );
    
    solver.update_weights();

    // 
    Vec<TF> bnds;
    TF mi = de->min_x();
    TF ma = de->max_x();
    for( TF p : pd->cell_boundaries() )
        bnds << min( ma, max( mi, p ) );

    OtResult res;
    res.barycenters = array_from_vec( pd->barycenters( *de ) );
    res.boundaries = array_from_vec( bnds );
    res.weights = array_from_vec( pd->get_weights() );
    return res;
}

OtResult ot_diracs_to_piecewise_constant( Array &dirac_positions, Array &target_mass_ratios, Array &density_positions, Array &density_values, OtParms &parms ) {
    auto de = RcPtr<PiecewiseConstantDensity<TF>>::New( 
        vec_from_array( density_positions ),
        vec_from_array( density_values )
    );

    return ot_solve( dirac_positions, target_mass_ratios, de, parms );
}

PYBIND11_MODULE(usdot_bindings, m) {
    pybind11::class_<OtResult>( m, "OtResult" )
        .def_readwrite( "error_message" , &OtResult::error_message, "" )
        .def_readwrite( "barycenters" , &OtResult::barycenters, "" )
        .def_readwrite( "boundaries" , &OtResult::boundaries, "" )
        .def_readwrite( "weights" , &OtResult::weights, "" )
        .def( "__repr__",
            []( const OtResult &a ) {
                return "{ error_message: '" + a.error_message + "', barycenters: [...], boundaries: [...], weights: [...] }";
            }
        )
        ;

    pybind11::class_<OtParms>( m, "OtParms" )
        .def( py::init<>() )
        .def_readwrite( "max_mass_ratio_error_target" , &OtParms::max_mass_ratio_error_target, "" )
        .def( "__repr__",
            []( const OtParms &a ) {
                return "{ max_mass_ratio_error_target: " + std::to_string( a.max_mass_ratio_error_target ) + " }";
            }
        )
        ;

    m.def( "ot_diracs_to_piecewise_constant", &ot_diracs_to_piecewise_constant );
}
/*
<%
setup_pybind11(cfg)
cfg['include_dirs'] += [ '../../../ext' ]
cfg['extra_compile_args'] += ['-std=c++20']
%>
*/