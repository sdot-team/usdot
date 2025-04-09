import cppimport.import_hook
from . import usdot_bindings
import numpy as np

OtResult = usdot_bindings.OtResult
OtParms = usdot_bindings.OtParms

# polynomial order is determined by `density_values.shape[ 1 ] - 1`. Each row gives the coefficients of a polynomial where x = 0 at the beginning
# of the interval and x = 1 at the end of the interval (intervals are determined by `density_positions`).
# def d2cap( dirac_positions, density_positions, density_values, global_mass_ratio = 1, relative_mass_ratios = None, ot_parms = None ):
#     """ optimal transport plan with src = diracs, dst = continuous piecewise affine density

#          `density_positions` are the boundaries of the intervals where the density is defined. 
#          `density_values` are the values of the density for each value of `density_positions`.
#     """
#     assert len( density_positions ) == len( density_values )
#     assert len( density_positions )

#     if not isinstance( density_values, np.ndarray ):
#         density_values = np.array( density_values )

#     dv = np.empty( [ len( density_positions ) - 1, 2 ] )
#     dv[ :, 0 ] = density_values[ :-1 ]
#     dv[ :, 1 ] = density_values[ 1: ] - density_values[ :-1 ]

#     return d2p( dirac_positions, density_positions, dv, global_mass_ratio = global_mass_ratio, relative_mass_ratios = relative_mass_ratios, ot_parms = ot_parms )


def from_p1_grid( dirac_positions, density_values, density_beg = None, density_end = None, global_mass_ratio = 1, relative_mass_ratios = None, ot_parms = None ):
    """ 
        optimal transport plan with src = diracs, dst = piecewise polynomial density, with values on a grid between density_beg and density_end
    """
    if ot_parms is None:
        ot_parms = usdot_bindings.OtParms()

    if relative_mass_ratios is None:
        relative_mass_ratios = np.ones( len( dirac_positions ) )

    if density_beg is None:
        density_beg = 0
    if density_end is None:
        density_end = density_beg + len( density_values ) - 1

    return usdot_bindings.from_p1_grid( 
        np.ascontiguousarray( dirac_positions ),
        global_mass_ratio,
        np.ascontiguousarray( relative_mass_ratios ),
        np.ascontiguousarray( density_values ),
        density_beg, 
        density_end,
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
