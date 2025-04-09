#pragma once

#include <tl/support/containers/Vec.h>
#include <tl/support/containers/Opt.h>

namespace usdot {

template<class TF>
struct FastGridSolverInput {
    Vec<TF> dirac_positions                 = {}; ///< mandatory
    Vec<TF> dirac_weights                   = {}; ///< mandatory
    Vec<TF> dirac_masses                    = {};

    Vec<TF> density_values                  = {}; ///< mandatory
    Opt<TF> beg_x_density                   = {}; ///< 0 by default
    Opt<TF> end_x_density                   = {}; ///< 1 by default

    TF      global_mass_ratio               = 1;

    TF      starting_contrast_ratio         = 1e-6;
    TF      target_contrast_ratio           = 1e-6;
    TF      min_dirac_separation            = 1e-6;
    TF      target_l2_error                 = 1e-5;
    bool    throw_if_error                  = 1;
    bool    multithread                     = 0;
    int     verbosity                       = 1;
};

} // namespace usdot

