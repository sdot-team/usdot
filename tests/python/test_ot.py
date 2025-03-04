import numpy as np
import usdot

def test_affine():
    res = usdot.d2cap( np.linspace( 0, 1, 100 ), [ 0, 1, 2 ], [ 0, 1, 0 ] )
    # print( "ba", res.barycenters )
    # print( "bo", res.boundaries )
    # print( "we", res.weights )
    # print( "ma", res.masses )
    # usdot.plot( res )

    p = np.linspace( 0, 2, int( 1e6 ) )
    v = np.minimum( p, 2 - p )[ :-1 ]
    ref = usdot.d2p( np.linspace( 0, 1, 100 ), p, v )

    assert np.max( np.abs( res.barycenters - ref.barycenters ) ) < 3e6

def test_partial():
    res = usdot.d2p( np.linspace( 0, 5, 20 ), [ 0, 1, 2, 3, 4, 5 ], [ 1, 0, 10, 0, 1 ], 0.85 )
    # print( "re", res.norm_2_residual_history )
    # print( "ba", res.barycenters )
    print( "bo", res.boundaries )
    # print( "we", res.weights )
    print( "ma", res.masses )
    usdot.plot( res )


# test_affine()
test_partial()
