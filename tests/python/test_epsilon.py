import numpy as np
import usdot

ot_parms = usdot.OtParms()
# ot_parms.epsilon = 1e-3

diracs = np.linspace( 0, 1, 30 )

dens_p = [ 0, 0.5, 1 ]
dens_v = [ 0, 1, 0 ]

res = usdot.d2cap( diracs, dens_p, dens_v, mass_ratio = 1, ot_parms = ot_parms )
print( "hi", res.norm_2_residual_history )
print( "ba", res.barycenters )
print( "bo", res.boundaries )
print( "we", res.weights )
print( "ma", res.masses )
usdot.plot( res )
