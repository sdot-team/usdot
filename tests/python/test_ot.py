import numpy as np
import usdot

res = usdot.d2p( dirac_positions = np.linspace( 0, 1, 10 ), density_positions = [ 0, 0.4, 0.6, 1 ], density_values = [ 1, 0, 1 ], global_mass_ratio = 1 )
# print( "hi", res.norm_2_residual_history )
print( "ba", res.barycenters )
print( "bo", res.boundaries )
# print( "we", res.weights )
# print( "ma", res.masses )
usdot.plot( res )
