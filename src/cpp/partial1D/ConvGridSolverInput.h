#pragma once

#include <vector>

namespace usdot {

template<class TF>
struct ConvGridSolverInput {
    using   Vec                             = std::vector<TF>;
    
    Vec     dirac_positions                 = {}; ///< mandatory
    Vec     dirac_weights                   = {}; ///< mandatory
    Vec     dirac_masses                    = {}; ///< optional

    Vec     density_values                  = {}; ///< mandatory
    TF      beg_x_density                   = 0; ///< 
    TF      end_x_density                   = 1; ///< 

    TF      global_mass_ratio               = 1;

    TF      starting_filter_value           = 1e-6;
    TF      target_filter_value             = 1e-6;

    TF      min_dirac_separation            = 1e-6;
    TF      target_l2_error                 = 1e-5;

    bool    throw_if_error                  = 1;
    bool    multithread                     = 0;
    int     verbosity                       = 1;
};

} // namespace usdot

