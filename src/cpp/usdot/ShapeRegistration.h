#pragma once

#include <tl/support/containers/Opt.h>
#include "AffineTransformation.h"
#include "VtkOutput.h"

namespace usdot {

    
/**
*/
template<class TF,int dim>
class ShapeRegistration {
public:
    using         Tr                       = AffineTransformation<TF,dim>;
    using         Pt                       = Vec<TF,dim>;
    struct        Triangle                 { Vec<Pt,3> points; };
           
    void          display_barycenters_vtk  ( VtkOutput &vo );
    void          display_diracs_vtk       ( VtkOutput &vo );
    void          display_shape_vtk        ( VtkOutput &vo );
    void          display_vtk              ( Str filename, PI num = 0 );
  
    static void   load_triangles           ( Vec<Triangle> &triangles, Str filename );
    static void   load_points              ( Vec<Pt> &points, Str filename );
    void          load_shape               ( Str filename );
       
    auto          projected_density        ( Pt proj_dir ) -> std::tuple<TF,TF,Vec<TF>>;
    Vec<Pt>       projection_dirs          (); ///<
      
    void          compute_new_diracs       ( PI nb_iter = 50 ); ///<
    Vec<TF>       delta_for_dir            ( Pt dir );
  
    static void   glot                     ( Vec<TF> xs, auto &&...funcs );    

    Vec<PI>       solve_steps;
    Vec<Pt>       new_diracs;              ///< new dirac positions
       
    TF            shape_thickness          = 1e-6;
    Vec<Triangle> shape_triangles;     
    Vec<Pt>       shape_points;     
       
    Tr            transformation;     
    TF            mass_ratio               = 1;
    Vec<Pt>       diracs;     
       
    PI            nb_projection_dirs       = 25;
    PI            nb_bins                  = 100;
};  


} // namespace usdot

#include "ShapeRegistration.cxx"
