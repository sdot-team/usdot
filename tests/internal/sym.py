from matplotlib import pyplot as plt
from scipy.interpolate import pade
import numpy as np
import sympy


class System:
    def __init__( self, coords, density, error_function = None ):
        if error_function is None:
            error_function = lambda a, b: sympy.log( a / b )
        self.position_symbol = sympy.Symbol( "x" )
        self.error_function = error_function
        self.nb_cells = len( coords )
        self.coords = coords

        self.density = density( self.position_symbol )

        self.weights_symbols = sympy.Matrix( [ sympy.Symbol( f'w_{ n }' ) for n in range( self.nb_cells ) ] )
        self.weights_values = sympy.ones( self.nb_cells, 1 ) * 0

        self.exp_masses = sympy.ones( self.nb_cells, 1 ) * sympy.integrate( self.density, ( self.position_symbol, -sympy.oo, +sympy.oo ) ) / self.nb_cells

    def cells( self ):
        res = []
        res.append( -sympy.oo )

        for i in range( 1, self.nb_cells ):
            mid = self.coords[ i - 1 ] + self.coords[ i - 0 ]
            dlt = ( self.weights_symbols[ i - 1 ] - self.weights_symbols[ i - 0 ] ) / ( self.coords[ i - 0 ] - self.coords[ i - 1 ] )
            res.append( ( mid + dlt ) / 2 )
            res.append( ( mid + dlt ) / 2 )

        res.append( +sympy.oo )
        return sympy.Matrix( res ).reshape( self.nb_cells, 2 )

    def subs( self, value ):
        for s, d in zip( self.weights_symbols, self.weights_values ):
            value = value.subs( s, d )
        return value

    def errors( self ) -> tuple[ sympy.Matrix, int ]:
        res = []
        nb_arcs = 0
        cells = self.cells()
        for num_cell in range( self.nb_cells ):
            cur_area = sympy.integrate( self.density, ( self.position_symbol, cells[ num_cell, 0 ], cells[ num_cell, 1 ] ) )
            exp_area = self.exp_masses[ num_cell ]
            # res.append( sympy.log( sympy.Max( 1e-10, cur_area / exp_area ) ) )
            # x = cur_area / exp_area
            # res.append( x - 1 )
            res.append( self.error_function( cur_area, exp_area ) )
        return ( sympy.Matrix( res ), nb_arcs )

    def error( self ) -> sympy.Expr:
        res = 0
        for e in self.errors()[ 0 ]:
            res += e**2
        return res

    def newton_dir( self ):
        E, nb_arcs = self.errors()
        # print( E )

        M = np.zeros( [ self.nb_cells, self.nb_cells ] )
        V = np.zeros( [ self.nb_cells ] )
        for r in range( self.nb_cells ):
            for c in range( self.nb_cells ):
                M[ r, c ] = self.subs( E[ r ].diff( self.weights_symbols[ c ] ) )
            V[ r ] = - self.subs( E[ r ] )

        # print( M )
        # print( V )
        if nb_arcs == 0:
            n = np.argmin( np.abs( M.diagonal() ) )
            M[ n, n ] += 1

        print( np.linalg.eigvals( M ) )

        return sympy.Matrix( np.linalg.solve( M, V ) )

def poly_coeffs( func: sympy.Expr, x: sympy.Expr, deg: int, c: sympy.Expr = 0 ):
    d = func
    coeffs = []
    for order in range( deg + 1 ):
        coeffs.append( d.subs( x, c ) / sympy.factorial( order ) )
        d = d.diff( x )
    print( coeffs )
    return coeffs

def poly_approx( func: sympy.Expr, x: sympy.Expr, deg: int, c: sympy.Expr = 0 ):
    res = 0
    for order, coeff in enumerate( poly_coeffs( func, x, deg, c ) ):
        res += coeff * ( x - c ) ** order
    return res

def pade_approx( func: sympy.Expr, x: sympy.Expr, deg_num: int, deg_den: int, c: sympy.Expr = 0 ):
    coeffs = [ float( c ) for c in poly_coeffs( func, x, deg_num + deg_den, c ) ] 
    p, q = pade( coeffs, deg_num )

    num = 0
    for order, coeff in enumerate( p.coefficients[ ::-1 ] ):
        num += coeff * ( x - c ) ** order

    den = 0
    for order, coeff in enumerate( q.coefficients[ ::-1 ] ):
        den += coeff * ( x - c ) ** order

    return num / den

def test_relax( density_function, error_function, max_x = 1, nb_cells = 7 ):
    """
        On prend 1 direction de Newton, et on regarde comment se comporte l'erreur réelle
        On compare ensuite avec le taylor de err( w + a * dir )
    """
    
    sys = System( np.linspace( -2, 2, nb_cells ), density_function, error_function )
    dir = sys.newton_dir()
    print( [ float( d ) for d in dir ] )

    alp = sympy.Symbol( "a" )
    err = sys.error().subs( zip( sys.weights_symbols, sys.weights_values + alp * dir ) )

    xs = np.concat( [
        np.linspace( 0.0, max_x, 30 ),
        # np.linspace( 1.1, 2.33, 10 )
    ] )
    plt.plot( xs, [ err.subs( alp, a ) for a in xs ], "--", label = "ref" )

    for num, den in [ ( 5, 2 ), ]:
        print( num, den )
        ape = pade_approx( err, alp, num, den )
        plt.plot( xs, [ ape.subs( alp, a ) for a in xs ], "+", label = f"p_{ num }_{ den }" )

    for order in [ 2, 4, 5, 6 ]:
        print( order )
        ape = poly_approx( err, alp, order )

        plt.plot( xs, [ ape.subs( alp, a ) for a in xs ], label = f"o{ order }" )
        print( order )

    plt.legend()
    plt.show()

# Prop: on fait une convolution avec 1/((x/std)^2+1) pour avoir un comportement un peu plus linéaire
#

# error_function = lambda a, b: sympy.log( sympy.Max( 1e-5, a / b ) )
test_relax(
    lambda x: sympy.exp( - ( x - 2 ) ** 2 ) + sympy.exp( - ( x + 2 ) ** 2 ),
    # lambda a, b: sympy.log( sympy.Max( 1e-5, a / b ) )
    lambda a, b: a - b,
    0.4
)

