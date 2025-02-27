from tabulate import tabulate
from sympy import * 
# init_printing(use_unicode=False)
# x = Symbol('x')
# y = Symbol('y')

# def integ( p, density, x_to_y, y_to_x ):
#     m = 0.5
#     a = Symbol( 'a' )

#     res = 0
#     for i in range( len( p ) ):
#         b = a + ( i + 0 ) * m / len( p )
#         e = a + ( i + 1 ) * m / len( p )
#         res += integrate( ( y_to_x( y ) - y_to_x( p[ i ] ) )**2 * density( y_to_x( y ) ), ( y, b, e ) )

#     print( diff( res, a ) )

# def test( density, x_to_y, y_to_x ):
#     integ( [ 0.5            ], density, x_to_y, y_to_x )
#     integ( [ 0.4, 0.5       ], density, x_to_y, y_to_x )
#     integ( [ 0.4, 0.45, 0.5 ], density, x_to_y, y_to_x )

# test(
#     lambda x: x / 2,
#     lambda x: x * x,
#     lambda y: y ** 0.5,
# )

def subs( e, a, b ):
    return Subs( e, a, b ).doit()

def err( diracs, mass_per_dirac ):
    dofs = [ Symbol( f'd_{ i }' ) for i in range( len( diracs ) ) ]
    vals = [ 2 * ( i == 0 ) for i in range( len( diracs ) ) ]

    weights = [ sum( dofs[ : i + 1 ] ) for i in range( len( diracs ) ) ]

    areas = []
    for i in range( len( diracs ) ):
        d1 = diracs[ i ]
        w1 = weights[ i ]
        b1 = d1 - w1 ** 0.5
        e1 = d1 + w1 ** 0.5
        if i:
            d0 = diracs[ i - 1 ]
            w0 = weights[ i - 1 ]
            e0 = d0 + w0 ** 0.5
            if Subs( e0 - b1, dofs, vals ).doit() >= 0:
                b1 = ( d0 + d1 ) / 2 + ( w0 - w1 ) / ( d1 - d0 ) / 2
        if i + 1 < len( diracs ):
            d2 = diracs[ i + 1 ]
            w2 = weights[ i + 1 ]
            b2 = d2 - w2 ** 0.5
            if Subs( e1 - b2, dofs, vals ).doit() >= 0:
                e1 = ( d1 + d2 ) / 2 + ( w1 - w2 ) / ( d2 - d1 ) / 2

        areas.append( e1 - b1 )

    for area in areas:
        print( subs( area, dofs, vals ) )
    m = []
    for area in areas:
        r = []
        for d in dofs:
            r.append( float( subs( area.diff( d ), dofs, vals ) ) )
        m.append( r ) 
    print( tabulate( m ) )

err( [ 0, 1, 2, 3, 4, 10, 11, 12 ], 2 )


