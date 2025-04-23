import numpy as np

def mat( t, s = 10 ):
    res = np.zeros( [ s, s ] )
    for i in range( s - 1 ):
        res[ i, i + 0 ] = 1
        res[ i, i + 1 ] = t - 1
        res[ s - 1, i ] = 1 - t
    res[ s - 1, s - 1 ] = 1
    return res

def sol( t, D ):
    s = len( D )
    # M = mat( t, s )

    a = 1 / ( 1 - t )
    # l = M[ s - 1, : ]
    x = 1 - t + t * D[ s - 1 ]
    for i in range( s - 2, 0, -1 ):
        # l += a * M[ i, : ]
        x += a * D[ i ]
        a = 1 + a / ( 1 - t )
    # l += a * M[ 0, : ]
    x += a * D[ 0 ]

    R = np.empty( [ s ] )
    R[ 0 ] = x / ( a + 1 - t )
    for i in range( 1, s ):
        R[ i ] = ( D[ i - 1 ] - R[ i - 1 ] ) / ( t - 1 )
    return R


s = 10
t = 0.1
M = mat( t, s )
D = np.random.random( [ s ] )

E = D.copy()
E[ -1 ] = 1 - t + t * D[ -1 ]

X = np.linalg.solve( M, E )
print( "X:", X )
print( "S:", sol( t, D ) )
