import matplotlib.pyplot as plt
import numpy as np
import usdot

parms = usdot.OtParms()
parms.verbosity = 2 

res = usdot.from_p1_grid( dirac_positions = np.random.random( [ 20 ] ), density_positions = np.linspace( 10, 11, 3, endpoint = True ), density_values = [1., .1, 1.], global_mass_ratio = 0.95, ot_parms = parms )
usdot.plot( res )
print( res )

# def system_with_extension( D, t, coeff_extension = 2 ):
#     s = len( D )
#     e = int( coeff_extension * s )
#     n = s + 2 * e
    
#     m = np.mean( D )

#     mat = np.zeros( [ n, n ] )
#     vec = np.zeros( [ n ] )
#     for i in range( n ):
#         c = t * s**1

#         # lapl
#         if i:
#             mat[ i, i - 1 ] -= c
#             mat[ i, i - 0 ] += c
#         if i + 1 < n:
#             mat[ i, i + 0 ] += c
#             mat[ i, i + 1 ] -= c
        
#         # id
#         if i >= e and i < e + s:
#             mat[ i, i ] += 1 - t
#             vec[ i ] = D[ i - e ] * ( 1 - t )

#     print( mat )
    
#     Y = np.linalg.solve( mat, vec )
#     return Y[ e : -e ]

# def system_without_bnd_eqs( D, t ):
#     s = len( D )
    
#     # m = np.mean( D )

#     mat = np.zeros( [ s, s ] )
#     vec = np.zeros( [ s ] )
#     for i in range( s ):
#         c = t * s**1

#         mat[ i, i ] = 1 - t + 2 * c
#         if i:
#             mat[ i, i - 1 ] = - c
#         if i + 1 < s:
#             mat[ i, i + 1 ] = - c
        
#         vec[ i ] = D[ i - e ] * ( 1 - t )

#     Y0 = np.linalg.solve( mat, vec )
#     G0 = Y0[ 1: ] - Y0[ :-1 ]
 
#     dvec = vec.copy()
#     dvec[ 0 ] += 1
#     Y1 = np.linalg.solve( mat, dvec )
#     G1 = Y1[ 1: ] - Y1[ :-1 ]

#     dvec = vec.copy()
#     dvec[ -1 ] += 1
#     Y2 = np.linalg.solve( mat, dvec )
#     G2 = Y2[ 1: ] - Y2[ :-1 ]

#     nmat = np.empty( [ 2, 2 ] )
#     nvec = np.empty( [ 2 ] )
#     nmat[ 0, 0 ] = np.dot( G1, G1 )
#     nmat[ 0, 1 ] = np.dot( G1, G2 )
#     nmat[ 1, 0 ] = np.dot( G2, G1 )
#     nmat[ 1, 1 ] = np.dot( G2, G2 )
#     nvec[ 0 ] = np.dot( G0, G1 )
#     nvec[ 1 ] = np.dot( G0, G2 )

#     nsol = np.linalg.solve( nmat, nvec )
#     dvec = vec.copy()
#     dvec[ +0 ] -= nsol[ 0 ]
#     dvec[ -1 ] -= nsol[ 1 ]

#     return np.linalg.solve( mat, dvec )


# s = 3
# D = np.zeros( [ s ] )
# b = int( 0 * s / 10 )
# e = int( 8 * s / 10 )
# D[ b : e ] = 1
# D /= np.sum( D )
# plt.plot( system_with_extension( D, 1 ) )

# # for t in np.linspace( 0.1, 0.99, 10, endpoint = True ):
# #     plt.plot( system_with_extension( D, t ) )
# #     print( np.sum( system_with_extension( D, t ) ) )
# # plt.show()


# # def sol_from_D( t, D, b, e ):
# #     s = len( D )
# #     E = np.zeros( [ s + 1 ] )
# #     E[ 0:s ] = ( 1 - t ) * D
# #     E[ s ] = np.sum( D[ b : e] ) - 0.5 * ( D[ b ] + D[ e - 1 ] )
# #     return np.linalg.solve( mat( t, s, b, e ), E )[ :-1 ]
# #     # return sol( t, E )

# # def test_filter( D, l = 1 ):
# #     r = []
# #     o = D[ 0 ]
# #     tau = pow( 0.5, 1 / ( len( D ) * l ) )
# #     print( tau )
# #     for i in range( len( D ) ):
# #         r.append( o )
# #         o = ( 1 - tau ) * D[ i ] + tau * o

# #     s = []
# #     for i in range( len( D ) ):
# #         s.append( o )
# #         o = ( 1 - tau ) * D[ i ] + tau * o

# #     plt.plot( s )


# # s = 50
# # D = np.zeros( [ s ] )
# # b = int( 2 * s / 10 )
# # e = int( 8 * s / 10 )
# # D[ b : e ] = 1
# # test_filter( D )
# # plt.show()
# # # D[ 0 ] = s

# #     Y = sol_from_D( t, D, b, e )
# #     plt.plot( Y )
# # plt.show()
# # M = np.array([
# #     [  0,  0,  0 ],
# #     [ -1,  2, -1 ],
# #     [  0,  0,  0 ],
# # ])
# # M = np.array([
# #     [  0,  0,  0,  0 ],
# #     [ -1,  2, -1,  0 ],
# #     [  0, -1,  2, -1 ],
# #     [  0,  0,  0,  0 ],
# # ])
# # print( M )
# # print( M.T @ M )

# # def sol( t, E ):
# #     s = len( E )

# #     a = 1 / ( 1 - t )
# #     # x = 1 - t + t * D[ s - 1 ]
# #     x = E[ s - 1 ]
# #     for i in range( s - 2, 0, -1 ):
# #         x += a * E[ i ]
# #         a = 1 + a / ( 1 - t )
# #     x += a * E[ 0 ]

# #     R = np.empty( [ s ] )
# #     R[ 0 ] = x / ( a + 1 - t )
# #     for i in range( 1, s ):
# #         R[ i ] = ( E[ i - 1 ] - R[ i - 1 ] ) / ( t - 1 )
# #     return R

# # def ap_der( t, D, eps = 1e-4, N = 1 ):
# #     if N == 1:
# #         return ( sol_from_D( t + eps, D ) - sol_from_D( t, D ) ) / eps
# #     return ( ap_der( t + eps, D, eps, N - 1 ) - ap_der( t, D, eps, N - 1 ) ) / eps

# # def der( X, t, D ):
# #     s = len( X )
# #     F = D.copy()
# #     for i in range( s - 1 ):
# #         F[ i ] -= X[ i + 1 ]
# #         F[ -1 ] += X[ i ]
# #     F[ -1 ] -= 1
# #     return sol( t, F )

# # def der_next( X, t, N ):
# #     s = len( X )
# #     F = np.zeros( [ s ] )
# #     for i in range( s - 1 ):
# #         F[ i ] = - N * X[ i + 1 ]
# #         F[ -1 ] += N * X[ i ]
# #     return sol( t, F )


# # s = 5
# # t = 0.1
# # D = np.random.random( [ s ] )
# # D = np.array( [ s ] )
# # D /= np.sum( D )

# # M = mat( t, s )
# # E = t * D
# # E[ -1 ] += 1 - t
# # X = np.linalg.solve( M, E )

# # print( "X:", X )
# # print( "S:", sol( t, E ) )

# # d1 = der( X, t, D )
# # print( "D1:", ap_der( t, D ) )
# # print( "d1:", d1 )

# # d2 = der_next( d1, t, 2 )
# # print( "D2:", ap_der( t, D, N = 2 ) )
# # print( "d2:", d2 )

# # d3 = der_next( d2, t, 3 )
# # print( "D3:", ap_der( t, D, N = 3 ) )
# # print( "d3:", d3 )

# # for t in np.linspace( 0.0, 0.999, 10 ):
# #     plt.plot( sol_from_D( t, D ) )
# #     print( np.sum( sol_from_D( t, D ) ) )
# # plt.show()

# # n = 10
# # A = np.zeros( [ n + 1, n + 1 ] )
# # for i in range( n ):
# #     if i and i + 1 < n:
# #         A[ i, i ] = 2
# #     else:
# #         A[ i, i ] = 1
# #     if i:
# #         A[ i, i - 1 ] = -1
# #     if i + 1 < n:
# #         A[ i, i + 1 ] = -1
# #     A[ -1, i ] = 1
# #     A[ i, -1 ] = 1
# # A[ -1,  0 ] = 0.5
# # A[ -1, -2 ] = 0.5
# # A[  0, -1 ] = 0.5
# # A[ -2, -1 ] = 0.5

# # # E =np.array( [    0.111105,    0.111106,    0.111108,    0.111111,    0.111115,    -0.88889,    0.111116,    0.111113,    0.111111,     0.11111, 0 ] )
# # E =np.array( [ -0.111109,  -0.0779788,   0.0792141,  0.00441568,   0.0230846,  -0.0559587,  -0.0992613,   0.0598843,    0.059993,    0.124323, 0 ] )
# # print( np.linalg.solve( A, E ) )

# # def an( n ):
# #     return n * n / 2

# # def bn( n ):
# #     return n * n / 2 + 1

# # def solve( M, Y ):
# #     a0 = 0
# #     a1 = a0 + 0.5
# #     b0 = 1
# #     b1 = b0 + 0.5
# #     la = a0 * M[ 0, : ] + a1 * M[ 1, : ]
# #     lb = b0 * M[ 0, : ] + b1 * M[ 1, : ]
# #     ay = a0 * Y[ 0 ] + a1 * Y[ 1 ]
# #     by = b0 * Y[ 0 ] + b1 * Y[ 1 ]
# #     for n in range( 2, 10 ):
# #         a2 = 2 * a1 + 1 - a0
# #         b2 = 2 * b1 + 1 - b0
# #         la += a2 * M[ n, : ]
# #         lb += b2 * M[ n, : ]
# #         ay += a2 * Y[ n ]
# #         by += b2 * Y[ n ]

# #         a0 = a1
# #         a1 = a2
# #         b0 = b1
# #         b1 = b2

# #     la += M[ -1, : ]
# #     lb += M[ -1, : ]
# #     ay += Y[ -1 ]
# #     by += Y[ -1 ]
# #     print( "la", la )
# #     print( "lb", lb )
# #     # n, 
# #     # print( ( n - 2 ) * ( n - 1 ) * ( 2 * n - 3 ) / 12 + ( n - 1 ) * ( n - 1 ) / 4 )
# #     # print( ( n - 2 ) * ( n - 1 ) * ( 2 * n - 3 ) / 12 + ( n - 1 ) * ( n - 1 ) / 4 + ( n - 1 ) )
# #     # 
# #     # print( ay, by )
# #     n = M.shape[ 0 ] - 1
# #     print( n )
# #     Z = ( n - 2 ) * ( n - 1 ) * ( 2 * n - 3 ) / 12 + ( n - 1 ) * ( n - 1 ) / 4
    
# #     print( "sol:", np.linalg.solve( M, Y ) )
    
# #     lam = ( by - ay ) / ( n - 1 )
# #     print( "lam", lam )
# #     print( "Z", Z )

# #     X = np.zeros( [ n ] )
# #     X[ n - 1 ] = ( ay - Z * lam ) / ( n - 1 )
# #     X[ n - 2 ] = X[ n - 1 ] + lam / 2 - Y[ n - 1 ]
# #     for i in range( n - 3, -1, -1 ):
# #         X[ i ] = 2 * X[ i + 1 ] - X[ i + 2 ] + lam - Y[ i + 1 ]
# #     print( X )


# # # solve( A, E )
