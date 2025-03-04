from sympy import *

I0 = Symbol( "I0" )
l1 = Symbol( "l1" )
l2 = Symbol( "l2" )
c0 = Symbol( "c0" ) # position du 1er dirac
d1 = Symbol( "d1" ) # écart entre les diracs
c1 = c0 + d1
I1 = I0 + l1
I2 = I1 + l2

E = + I2 * ( I2 - 2 * c1 ) \
    + I0 * ( 2 * c0 - I0 ) \
    + 2 * I1 * ( c1 - c0 )

E = expand( E )
# E = collect( E, I0 )
# E = collect( E, l1 )
# E = collect( E, l2 )
print( E )

E = 2 * I0 * ( l1 + l2 ) - 2 * c0 * ( l1 + l2 ) - 2 * d1 * l2 + ( l1 + l2 ) * ( l1 + l2 )
E = ( l1 + l2 ) * ( 2 * ( I0 - c0 ) + l1 + l2 ) - 2 * d1 * l2

E( c0 ) = ( l1 + l2 ) * ( l1 + l2 ) - 2 * d1 * l2

# E = 2 * I0 * ( l1 + l2 ) - 2*c0*l1 - 2*c1*l2 + ( l1 + l2 ) ^ 2

# E = 2 * ( c1 - c0 ) * o1 + o2 * ( o2 + 2 * ( I0 - c1 ) )

# 2 * o1 *( c1 - c0 ) + o2 * ( o2 * 2*I0 - 2*c1)
