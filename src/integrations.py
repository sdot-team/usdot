from sympy import *

def scollect( expr, s ):
    loc = Symbol( 'loc', real = True )
    expr = expr.subs( s, loc )
    expr = collect( expr, loc )
    expr = expr.subs( loc, s )
    return expr

def lcollect( expr, l ):
    loc = Symbol( 'loc', real = True )
    for s in l:
        expr = expr.subs( s, loc )
    expr = collect( expr, loc )
    expr = expr.subs( loc, l[ 0 ] )
    return expr

def prim_area( c, w, x0, x1, v0, v1, a ):
    """ primitive pour calcul l'aire """
    cx = Symbol( "c_x")
    va = v0 + ( v1 - v0 ) / ( x1 - x0 ) * ( cx + sqrt( w ) * cos( a ) )
    area = integrate( - 2 * w * sin( a ) ** 2 * va, a )

    area = trigsimp( area )
    # area = collect( area, cx )
    area = collect( area, 6 * a * cx )
    area = collect( area, 6 * a * v0 )
    area = collect( area, 3 * cx * sin( 2 * a ) )
    area = collect( area, 3 * sqrt( w ) * sin( a ) )
    area = collect( area, sqrt( w ) * sin( 3 * a ) )
    area = collect( area, 3 * v0 * sin( 2 * a ) )
    area = collect( area, v1 - v0 )
    area = collect( area, x1 - x0 )
    area = collect( area, x1 - x0 )
    area = collect( area, 3 * cx )
    area = collect( area, 3 * v0 )
    area = collect( area, sqrt( w ) )

    n, d = fraction( area )
    n /= 6
    d /= 6

    # dx = Symbol( "d_x")
    # area = area.subs( x1 - x0, dx )
    # n = n.subs( x1 - x0, dx )
    # d = d.subs( x1 - x0, dx )
    # # area = area.subs( - x0 + x1, dx )
    # # area = area.subs( x0 - x1, - dx )
    # area = n / d

    return area.subs( cx, c - x0 )

def disp_area():
    center = Symbol( 'center', real = True )
    weight = Symbol( 'weight', real = True )
    x0 = Symbol( 'x0', real = True )
    x1 = Symbol( 'x1', real = True )
    v0 = Symbol( 'v0', real = True )
    v1 = Symbol( 'v1', real = True )
    a = Symbol( 'a', real = True )

    # primitve
    pa = prim_area( center, weight, x0, x1, v0, v1, a )
    print( ccode( pa ) )

def disp_darea():
    epsilon = Symbol( 'epsilon', real = True )
    center = Symbol( 'center', real = True )
    weight = Symbol( 'weight', real = True )
    coeff = Symbol( 'coeff', real = True )
    x0 = Symbol( 'x0', real = True )
    x1 = Symbol( 'x1', real = True )
    v0 = Symbol( 'v0', real = True )
    v1 = Symbol( 'v1', real = True )
    a = Symbol( 'a', real = True )
    x = Symbol( 'x', real = True )

    # primitve
    pa = prim_area( center, weight, x0, x1, v0, v1, a )

    #
    def disp_cond( cond, expr ):
        print( "if ( a_event_types[ num_event ] == EventType::" + cond + " ) {" )
        print( "    V[ na ] += " + ccode( coeff * expr ) + ";" )
        print( "    M( na, na ) += " + ccode( coeff * expr.diff( weight ) ) + ";" )
        print( "    return;" )
        print( "}" )

    #
    # DensityCut, InterCell
    disp_cond( "NewDensity", pa.subs( a, acos( ( x - center ) / sqrt( weight ) ) ) )
    disp_cond( "StartFlat" , pa.subs( a, pi - asin( epsilon / sqrt( weight ) ) ) )
    disp_cond( "EndFlat"   , pa.subs( a, asin( epsilon / sqrt( weight ) ) ) )
    disp_cond( "StartArc"  , pa.subs( a, pi ) )
    disp_cond( "EndArc"    , pa.subs( a, 0 ) )
    print( "P( a_event_types[ num_event ] );" )
    print( "TODO;" )


# def deps_darea_flat():
#     p = c - sqrt( weight - epsilon**2 )    

disp_darea()
