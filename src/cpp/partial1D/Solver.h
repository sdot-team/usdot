#pragma once

#include <tl/support/containers/Opt.h>
#include "PowerDiagram.h"

namespace usdot {

    
/**
*/
template<class TF>
class Solver {
public:
    using                   ConvolutionFactory              = std::function<RcPtr<Convolution<TF>>( TF width )>;             
    using                   IterationCallback               = std::function<bool()>;

    struct                  UpdateConvexHullDensityRatioPrm { Opt<TF> target_value; PI polynomial_order = 3; TF epsilon = 1e-4; };
    struct                  UpdateConvolutionWidthPrm       { Opt<TF> target_value; PI polynomial_order = 3; };
    struct                  InitializeWeightsPrm            { bool exact_position = false; };
    struct                  UpdateWeightsPrm                { TF relaxation = 1; PI max_newton_iterations = 300, min_newton_iterations = 0; };
    struct                  State                           { TF convex_hull_density_ratio, convolution_width; Vec<TF> weights; };
    struct                  Error                           {  };
                 
    /**/                    Solver                          ( RcPtr<PowerDiagram<TF>> power_diagram, RcPtr<Density<TF>> density, TF target_mass_ratio = 1 );

    // parameter update
    void                    update_convex_hull_density_ratio( UpdateConvexHullDensityRatioPrm parms = {} ); ///< must be done after an update_weights
    void                    update_convolution_width        ( UpdateConvolutionWidthPrm parms = {} ); ///< must be done after an update_weights
    
    // initialization
    void                    find_approx_weights_using_cdf   ( InitializeWeightsPrm parms = {} ); ///< make a rough guess for the weight using the cdf or an approximation of the cdf
    void                    find_exact_weights_using_cdf    ( InitializeWeightsPrm parms = {} ); ///< get the exact weights using cdf (easy to compute if global_mass_ratio == 1)

    // iterative procedures
    int                     update_weights_using_newton     ( UpdateWeightsPrm parms = {} ); ///< newton

    // all-in-one procedures
    void                    solve                           ();

    //
    RcPtr<Density<TF>>      current_density                 ();
    State                   get_state                       () const;
    void                    set_state                       ( const State &state );

    bool                    converged                       () const;

    // parameters
    TF                      max_mass_ratio_error_target     = 1e-8;
    ConvolutionFactory      convolution_factory             = {};
    IterationCallback       iteration_callback              = {};
    bool                    throw_if_error                  = 1;
    int                     verbosity                       = 1;
       
    // outputs
    Vec<TF>                 max_mass_ratio_error_history    = {};
    Vec<TF>                 norm_2_residual_history         = {};
    Vec<TF>                 norm_2_rhs_history              = {};
    PI                      solve_step                      = 0;
             
    // intermediate data     
    TF                      _convoluted_density_width       = 0; ///< 
    RcPtr<Density<TF>>      _convoluted_density;            ///< 
    SymmetricBandMatrix<TF> cholesky;                       ///<
   
    // state (input data that can be modified)
    TF                      convex_hull_density_ratio       = 0; ///< density = this->density + ( convex_hull_density - this->density ) * convex_hull_density_ratio
    TF                      convolution_width               = 0; ///<
             
    // input             
    Vec<TF>                 relative_mass_ratios            = {}; ///< relative mass ratio for each dirac (each dirac mass = global_mass_ratios * relative_mass_ratios[ num_dirac ])
    TF                      global_mass_ratio               = 1; ///< mass ratio for each dirac
    RcPtr<PowerDiagram<TF>> power_diagram;                  ///<
    RcPtr<Density<TF>>      density;                        ///<
};


} // namespace usdot

#include "Solver.cxx"
