import cppimport.import_hook
from . import usdot_bindings
import numpy as np

def d2p( dirac_positions, density_positions, density_values, mass_ratio = None, dirac_mass_ratios = None, ot_parms = None ):
    """ optimal transport plan with src = diracs, dst = piecewise polynomial density

       polynomial order is determined by `density_values.shape[ 1 ] - 1`. Each row gives the coefficients of a polynomial with x = 0 at the beginning
         of the interval and x = 1 at the end of the interval (intervals are determined by `density_positions`).
    """
    if ot_parms is None:
        ot_parms = usdot_bindings.OtParms()

    if dirac_mass_ratios is None:
        if mass_ratio is None:
           mass_ratio = 1.0
        dirac_mass_ratios = np.ones( len( dirac_positions ) ) / len( dirac_positions )

    return usdot_bindings.ot_diracs_to_piecewise_constant( 
        np.ascontiguousarray( dirac_positions ),
        np.ascontiguousarray( dirac_mass_ratios ),
        np.ascontiguousarray( density_positions ),
        np.ascontiguousarray( density_values ),
        ot_parms
    )

