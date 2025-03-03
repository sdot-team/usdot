import numpy as np
import usdot

ot_parms = usdot.OtParms()
# ot_parms.epsilon = 1e-3
ot_parms.verbosity = 2

diracs = np.linspace( 0, 1.3, 20 )

dens_p = [ 0, 1, 2, 3 ]
dens_v = [ 1, 1e-2, 1 ]

res = usdot.d2p( diracs, dens_p, dens_v, mass_ratio = 0.5, ot_parms = ot_parms )
# print( "hi", res.norm_2_residual_history )
# print( "ba", res.barycenters )
print( "bo", res.boundaries )
# print( "we", res.weights )
# print( "ma", res.masses )
usdot.plot( res )

# from matplotlib import pyplot
# pyplot.plot( 1 / ( res.boundaries[ :, 1 ] - res.boundaries[ :, 0 ] ) )
# pyplot.show()
