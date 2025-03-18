import sys, os
sys.path.append( os.path.join( os.path.dirname( __file__ ), "src", "python" ) )

import matplotlib.pyplot as plt
import numpy as np
import usdot

# print( np.random.get_state() )

xs = np.linspace( 0, 1, 1000 )
ys = np.full( 1000, 1e-3 )

s0 = 0.3
s1 = 0.9
for n, x in enumerate( xs ):
    if x >= s0 and x <= s1:
        ys[ n ] = 0.5 - 0.5 * np.cos( ( x - s0 ) / ( s1 - s0 ) * 2 * np.pi )
    if x >= .09 and x <= .11:
        ys[ n ] = 1
    if x >= .19 and x <= .21:
        ys[ n ] = 1

nb_diracs = 10
diracs = np.random.random( nb_diracs )

# ot with a continous piecewise affine function
res = usdot.d2cap( diracs, xs, ys )

# for x0, x1 in res.boundaries:
#     xs = [ x0, x1, x1, x0, x0 ]
#     ys = [ y1, y1, y0, y0, y1 ]
#     plt.plot( xs, ys )
for n in range( nb_diracs ):
    # plt.plot( res.barycenters, np.zeros_like( res.barycenters ), '+' )
    b0 = res.boundaries[ n ][ 0 ]
    b1 = res.boundaries[ n ][ 1 ]
    ba = res.barycenters[ n ]

    lxs = []
    lys = []
    for x, y in zip( xs, ys ):
        if x >= b0 and x <= b1:
            lxs.append( x )
            lys.append( y )
    lxs = [ lxs[ 0 ] ] + lxs + [ lxs[ -1 ] ]
    lys = [ 0 ] + lys + [ 0 ]

    P0 = np.array( [ diracs[ n ], 1.5 ] )
    P1 = np.array( [ ba, lys[ np.argmin( abs( ba - lxs ) ) ] ] )
    plt.annotate( "", xy = P1, xytext = P0, arrowprops = dict( 
        facecolor = 'lightgrey', 
        width = 0.001, 
        shrink = 0.04, 
        headwidth = 4, 
        # linestyle = (0, (1, 10)) 
    ) )
    # plt.arrow( x0[ 0 ], x1[ 0 ], x0[ 1 ] - x0[ 0 ], x1[ 1 ] - x1[ 0 ], head_width = 0.015 )
    # plt.plot( [ diracs[ n ], ba ], [ 1.45, lys[ np.argmin( abs( ba - lxs ) ) ] ] )

    plt.plot( lxs, lys )

plt.plot( diracs, np.zeros_like( res.barycenters ) + 1.5, '+', markersize=10, markeredgewidth=2 )
# plt.plot( res.barycenters, np.zeros_like( res.barycenters ), '+' )
# plt.axis( 'equal' )
plt.show()
