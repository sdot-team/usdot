from pysdot.domain_types import ConvexPolyhedraAssembly
from pysdot.radial_funcs import RadialFuncInBall
from pysdot import OptimalTransport
import numpy as np

with open( "../spot/Datasets/Pointsets/3D/mumble_sitting_100000.pts" ) as f:
    a = []
    for l in f.readlines():
        t, x, y, z = l.strip().split( ' ' )
        p = np.array( [ float( x ), float( y ), float( z ) ] )
        d = np.array( [ 0.42546746810417513, 0.35717031177892833, -0.8315087503861676 ] )
        o = np.array( [ -26.23, 1.35, 50.83 ] )
        if np.dot( d, p ) < np.dot( d, o ):
            a.append( [ float( x ), float( y ), float( z ) ] )
    base_pts = np.array( a )

base_pts = ( base_pts + 50 ) / 150
print( np.min( base_pts ) )
print( np.max( base_pts ) )

target_radius = 1e-4

# diracs
domain = ConvexPolyhedraAssembly()
domain.add_box( [ 0, 0, 0 ], [ 1, 1, 1 ] )

ot = OptimalTransport( domain )
# ot.set_weights( np.ones( base_pts.shape[ 0 ] ) * target_radius ** 3 )
# ot.set_masses( np.ones( base_pts.shape[ 0 ] ) * np.pi * target_radius ** 3 )
ot.set_positions( base_pts )
ot.verbosity = 2

# solve
ot.adjust_weights( relax = 0.5 )

# # display
# # ot.display_vtk( "results/pd.vtk" )

# # print( ot.pd.display_html() )

