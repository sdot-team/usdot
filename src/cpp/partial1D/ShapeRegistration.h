#pragma once

#include <tl/support/containers/Opt.h>
#include "AffineTransformation.h"
#include "PowerDiagram.h"
#include "VtkOutput.h"

namespace usdot {

    
/**
*/
template<class TF,int dim>
class ShapeRegistration {
public:
    using       Tr                       = AffineTransformation<TF,dim>;
    using       Pt                       = Vec<TF,dim>;
    struct      Triangle                 { Vec<Pt,3> points; };
         
    void        display_displacements_vtk( VtkOutput &vo );
    void        display_diracs_vtk       ( VtkOutput &vo );
    void        display_shape_vtk        ( VtkOutput &vo );
     
    static void load_triangles           ( Vec<Triangle> &triangles, Str filename );
    static void load_points              ( Vec<Pt> &points, Str filename );
    void        load_shape               ( Str filename );
     
    void        update_displacements     (); ///<
    auto        projected_density        ( Pt proj_dir ) -> std::pair<Vec<TF>,Vec<TF>>;
    Vec<Pt>     get_proj_dirs            (); ///<
     
    Vec<TF>     delta_for_dir            ( Pt dir );
    static void glot                     ( Vec<TF> xs, auto &&...funcs );    

    Vec<Pt>     displacements;           ///< displacements of the diracs
     
    TF          shape_thickness          = 1e-6;
    Vec<Pt>     shape_triangles;     
    Vec<Pt>     shape_points;     
     
    Tr          transformation;     
    TF          mass_ratio               = 1;
    Vec<Pt>     diracs;     
     
    PI          nb_proj_dirs             = 100;
    PI          nb_bins                  = 300;
};


} // namespace usdot

#include "ShapeRegistration.cxx"
