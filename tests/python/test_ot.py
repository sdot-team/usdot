import numpy as np
import unittest
import usdot

class TestOt( unittest.TestCase ):
    def test_basic( self ):
        """ simple [0,1] ot """
        res = usdot.from_p1_grid( dirac_positions = np.linspace( 0, 1, 10 ), density_values = [ 1, 1 ] )

        self.assertTrue( np.allclose( res.boundaries[ :, 0 ], np.linspace( 0.0, 0.9, 10 ), atol = 1e-6 ) )
        self.assertTrue( np.allclose( res.boundaries[ :, 1 ], np.linspace( 0.1, 1.0, 10 ), atol = 1e-6 ) )
        self.assertTrue( np.allclose( res.barycenters, np.linspace( 0.05, 0.95, 10 ), atol = 1e-6 ) )
        self.assertTrue( np.allclose( res.masses, 0.1, atol = 1e-6 ) )

    def test_scale( self ):
        """ scaled density positions """
        res = usdot.from_p1_grid( dirac_positions = np.linspace( 0, 1, 10 ), density_values = [ 1, 1 ], density_beg = -1, density_end = +1 )

        self.assertTrue( np.allclose( res.boundaries[ :, 0 ], np.linspace( -1.0, +0.8, 10 ), atol = 1e-6 ) )
        self.assertTrue( np.allclose( res.boundaries[ :, 1 ], np.linspace( -0.8, +1.0, 10 ), atol = 1e-6 ) )
        self.assertTrue( np.allclose( res.barycenters, np.linspace( -0.9, +0.9, 10 ), atol = 1e-6 ) )
        self.assertTrue( np.allclose( res.masses, 0.2, atol = 1e-6 ) )

    def test_same_pos( self ):
        """ scaled density positions """
        res = usdot.from_p1_grid( dirac_positions = [ 0, 0.5, 0.5, 1 ], density_values = [ 1, 1 ], density_beg = 0, density_end = 1 )
 
        self.assertTrue( np.allclose( res.boundaries[ :, 0 ], [ 0, 0.25, 0.25, 0.75 ], atol = 1e-6 ) )
        self.assertTrue( np.allclose( res.boundaries[ :, 1 ], [ 0.25, 0.75, 0.75, 1 ], atol = 1e-6 ) )
        self.assertTrue( np.allclose( res.barycenters, [ 0.125, 0.5, 0.5, 0.875 ], atol = 1e-6 ) )
        self.assertTrue( np.allclose( res.masses, 0.25, atol = 1e-6 ) )

    # def test_affine( self ):
    #     # ot with a continous piecewise affine function
    #     res = usdot.d2cap( np.linspace( 0, 1, 100 ), [ 0, 1, 2 ], [ 0, 1, 0 ] )

    #     # same thing with a piecewise constant approximation
    #     p = np.linspace( 0, 2, int( 1e6 ) )
    #     v = np.minimum( p, 2 - p )[ :-1 ]
    #     ref = usdot.d2p( np.linspace( 0, 1, 100 ), p, v )

    #     # comparison
    #     self.assertTrue( np.allclose( res.barycenters, ref.barycenters, atol = 3e-6 ) )

    # def test_partial( self ):
    #     # partial ot (just to test the convergence with "obstacles" and null densities)
    #     res = usdot.d2p( np.linspace( 0, 5, 20 ), [ 0, 1, 2, 3, 4, 5 ], [ 1, 0, 10, 0, 1 ], 0.85 )
    #     self.assertTrue( np.allclose( res.masses, 0.51, atol = 1e-6 ) )

    # def test_smurf( self ):
    #     # convex_hull_density_ratio: 0.0724313
    #     # global_mass_ratio: 0.984
    #     # current_density(): {
    #     #   sub_densities: [
    #     #     [ 0.0724313, { x0: -43.7729, x1: 34.6651, h: 2000 } ]

    #     # for i in range(1,len(d)):
    #     #     if d[ i ] < d[ i - 1 ] + 1e-2:
    #     #          d[ i ] = d[ i - 1 ] + 1e-2

    #     d = np.array( d )

    #     import matplotlib.pyplot as plt
    #     plt.hist( d, bins = 128 )
    #     plt.plot( xs, np.array(ys + [ ys[-1] ])/30 )
    #     plt.show()

    #     # partial ot (just to test the convergence with "obstacles" and null densities)
    #     # p = usdot.OtParms()
    #     # p.verbosity = 2

    #     # res = usdot.d2p( d, xs, ys, 0.984, ot_parms = p )
    #     # print( res.norm_2_residual_history )
    #     # self.assertTrue( np.allclose( res.masses, 0.51, atol = 1e-6 ) )

    # def test_epsilon( self ):
    #     ot_parms = usdot.OtParms()
    #     ot_parms.epsilon = 1e-3
    #     ot_parms.verbosity = 2

    #     diracs = np.linspace( 0, 1.3, 20 )

    #     dens_p = [ 0, 1, 2, 3 ]
    #     dens_v = [ 1, 1e-2, 1 ]

    #     res = usdot.d2p( diracs, dens_p, dens_v, mass_ratio = 0.5, ot_parms = ot_parms )
    #     # print( "hi", res.norm_2_residual_history )
    #     # print( "ba", res.barycenters )
    #     print( "bo", res.boundaries )
    #     # print( "we", res.weights )
    #     # print( "ma", res.masses )
    #     usdot.plot( res )

if __name__ == '__main__':
    unittest.main()

