import numpy as np
import usdot

res = usdot.d2p( np.linspace( 0, 1, 10 ), np.linspace( 0, 1, 4 ), [ 1, 1, 1 ] )
print( res.barycenters )
