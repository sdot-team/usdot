#pragma once

#include <vector>

namespace usdot {

template<class TF>
struct GradientSolverInput {
    using  TV                             = std::vector<TF>;
    
    TV     dirac_positions                = {}; ///< mandatory
    TV     dirac_weights                  = {}; ///< mandatory
    TV     dirac_masses                   = {}; ///< optional

    TV     density_values                 = {}; ///< mandatory
    TF     beg_x_density                  = 0; ///< 
    TF     end_x_density                  = 1; ///< 

    TF     global_mass_ratio              = 1;

    TV     filter_values                  = { 1e-3, 1e-2, 1e-1, 0.5 };

    TF     min_dirac_separation           = 1e-6;
    TF     target_l2_error                = 1e-5;

    bool   throw_if_error                 = 1;
    bool   multithread                    = 0;
    int    verbosity                      = 1;
};

} // namespace usdot

