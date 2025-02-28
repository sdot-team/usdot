import cppimport.import_hook
from . import usdot_bindings
import numpy as np

OtResult = usdot_bindings.OtResult
OtParms = usdot_bindings.OtParms

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
        dirac_mass_ratios = np.full( len( dirac_positions ), mass_ratio / len( dirac_positions ) )

    return usdot_bindings.ot_diracs_to_piecewise_constant( 
        np.ascontiguousarray( dirac_positions ),
        np.ascontiguousarray( dirac_mass_ratios ),
        np.ascontiguousarray( density_positions ),
        np.ascontiguousarray( density_values ),
        ot_parms
    )

def plot( ot_result ):
    """ plot the optimal transport result """
    import matplotlib.pyplot as plt

    dx = max( ot_result.boundaries ) - min( ot_result.boundaries )
    y0 = - dx / 30
    y1 = + dx / 30
    for n in range( 0, ot_result.boundaries.size, 2 ):
        x0, x1 = ot_result.boundaries[ n + 0 ], ot_result.boundaries[ n + 1 ]
        xs = [ x0, x1, x1, x0, x0 ]
        ys = [ y1, y1, y0, y0, y1 ]
        plt.plot( xs, ys )

    plt.plot( ot_result.barycenters, np.zeros_like( ot_result.barycenters ), '+' )

    plt.axis( 'equal' )
    plt.show()
