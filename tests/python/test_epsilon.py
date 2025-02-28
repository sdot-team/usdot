import numpy as np
import usdot

ot_parms = usdot.OtParms()
# ot_parms.epsilon = 1e-3

dip = np.linspace( 0, 1, 30 )
dep = [ 0, 0.4, 0.6, 1 ]
dev = [ 1, 0.001, 1 ]

res = usdot.d2p( dip, dep, dev, mass_ratio = 0.6, ot_parms = ot_parms )
print( "hi", res.norm_2_residual_history )
print( "ba", res.barycenters )
print( "bo", res.boundaries )
print( "we", res.weights )
print( "ma", res.masses )
usdot.plot( res )
