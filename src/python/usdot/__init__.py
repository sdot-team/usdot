import cppimport.import_hook
from . import usdot_bindings
import numpy as np

OtResult = usdot_bindings.OtResult
OtParms = usdot_bindings.OtParms

def from_p1_grid( dirac_positions, density_positions, density_values, global_mass_ratio = 1, relative_mass_ratios = None, ot_parms = None ):
    """ 
        optimal transport plan with src = diracs, dst = piecewise polynomial density, with values on a grid between density_beg and density_end
    """
    if ot_parms is None:
        ot_parms = usdot_bindings.OtParms()

    if relative_mass_ratios is None:
        relative_mass_ratios = np.ones( len( dirac_positions ) )

    return usdot_bindings.from_p1_grid( 
        np.ascontiguousarray( dirac_positions ),
        global_mass_ratio,
        np.ascontiguousarray( relative_mass_ratios ),
        np.ascontiguousarray( density_positions ),
        np.ascontiguousarray( density_values ),
        ot_parms
    )

def plot( ot_result ):
    """ plot the optimal transport result """
    import matplotlib.pyplot as plt

    dx = np.ptp( ot_result.boundaries )
    y0 = - dx / 30
    y1 = + dx / 30
    for x0, x1 in ot_result.boundaries:
        xs = [ x0, x1, x1, x0, x0 ]
        ys = [ y1, y1, y0, y0, y1 ]
        plt.plot( xs, ys )

    plt.plot( ot_result.barycenters, np.zeros_like( ot_result.barycenters ), '+' )

    plt.axis( 'equal' )
    plt.show()
