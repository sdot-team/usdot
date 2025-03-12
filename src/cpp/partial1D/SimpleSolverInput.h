#pragma once

#include <tl/support/containers/Vec.h>
#include <tl/support/containers/Opt.h>

namespace usdot {

template<class TF>
struct SimpleSolverInput {
    Vec<TF> dirac_positions                 = {}; ///< mandatory
    Vec<TF> dirac_masses                    = {};

    Vec<TF> density_values                  = {};
    Opt<TF> beg_x_density                   = {};
    Opt<TF> end_x_density                   = {};

    TF      global_mass_ratio               = 1;

    TF      max_mass_ratio_error_target     = 1e-5;
    TF      relative_dirac_separation       = 1e-6;
    TF      starting_contrast_ratio         = 1e-2;
    bool    throw_if_error                  = 1;
    int     verbosity                       = 1;
};

} // namespace usdot

