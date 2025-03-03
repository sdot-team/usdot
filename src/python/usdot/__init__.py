import cppimport.import_hook
from . import usdot_bindings
import numpy as np

OtResult = usdot_bindings.OtResult
OtParms = usdot_bindings.OtParms

def d2cap( dirac_positions, density_positions, density_values, mass_ratio = None, dirac_mass_ratios = None, ot_parms = None ):
    """ optimal transport plan with src = diracs, dst = continuous piecewise affine density

         `density_positions` are the boundaries of the intervals where the density is defined. 
         `density_values` are the values of the density for each value of `density_positions`.
    """
    assert len( density_positions ) == len( density_values )
    assert len( density_positions )

    if not isinstance( density_values, np.ndarray ):
        density_values = np.array( density_values )

    dv = np.empty( [ len( density_positions ) - 1, 2 ] )
    dv[ :, 0 ] = density_values[ :-1 ]
    dv[ :, 1 ] = density_values[ 1: ] - density_values[ :-1 ]

    print( density_positions, dv )

    return d2p( dirac_positions, density_positions, dv, mass_ratio = mass_ratio, dirac_mass_ratios = dirac_mass_ratios, ot_parms = ot_parms )


def d2p( dirac_positions, density_positions, density_values, mass_ratio = None, dirac_mass_ratios = None, ot_parms = None ):
    """ optimal transport plan with src = diracs, dst = piecewise polynomial density

       polynomial order is determined by `density_values.shape[ 1 ] - 1`. Each row gives the coefficients of a polynomial where x = 0 at the beginning
         of the interval and x = 1 at the end of the interval (intervals are determined by `density_positions`).
    """
    if ot_parms is None:
        ot_parms = usdot_bindings.OtParms()

    if dirac_mass_ratios is None:
        if mass_ratio is None:
           mass_ratio = 1.0
        dirac_mass_ratios = np.full( len( dirac_positions ), mass_ratio / len( dirac_positions ) )

    return usdot_bindings.ot_diracs_to_piecewise_polynomial( 
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
