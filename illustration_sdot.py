import sys, os
sys.path.append( os.path.join( os.path.dirname( __file__ ), "src", "python" ) )

import matplotlib.pyplot as plt
import numpy as np
# import usdot

def avec_fleches():
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

def gaussian( x, c, s ):
    return np.exp( - ( x - c )**2 / ( 2 * s**2 ) )
                  
def discretisation( x, a, n ):
    d = np.linspace( x[ 0 ], x[ -1 ], n )
    for i in range( 3 ):
        res = usdot.from_p1_grid( d, x, a )
        d = res.barycenters
    return d

def affichage_dirac_comme_fonction( d, l , w, c ):
    y = 1 / ( d[ 1: ] - d[ :-1 ] )
    y /= len( y )
    plt.plot( d[ 1: ], y, l, linewidth = w, color = c )

def animation_F2F():
    n = 3000
    a = np.zeros( [ n ] )
    b = np.zeros( [ n ] )
    x = np.linspace( 0, 1, n, endpoint = True )
    # a += gaussian( x, 0.2, 0.05 )
    # b += gaussian( x, 0.8, 0.05 )
    a += .001 + 2 * gaussian( x, 0.3, 0.1 )
    a += .001 + 1 * gaussian( x, 0.7, 0.2 )
    b += .001 + 2 * gaussian( x, 0.4, 0.1 )
    b += .001 + 1 * gaussian( x, 0.5, 0.2 )
    a /= np.sum( a )
    b /= np.sum( b )

    da = discretisation( x, a, n )
    db = discretisation( x, b, n )

    # ax.set_ylim([ymin, ymax]
    x1 = list( np.linspace( 0, 1, 20 ) )
    x2 = list( np.linspace( 1, 0, 20 ) )
    for i, t in enumerate( x1 + x2 ):
        ax = plt.gca()
        ax.set_xlim( [ 0, 1 ] )
        # ax.set_ylim( [ 0.5, 2 ] )
        affichage_dirac_comme_fonction( da, '--', 1, "red" )
        affichage_dirac_comme_fonction( db, '-.', 1, "blue" )
        affichage_dirac_comme_fonction( ( 1 - t ) * da + t * db, "-", 2, "black" )
        # plt.show()
        plt.savefig( f'f2f_{ 100 + i }.png', dpi = 300 )
        plt.cla()

    # parms = usdot.OtParms()
    # parms.verbosity = 2 

    # res = usdot.from_p1_grid( x, x, b, relative_mass_ratios = a, ot_parms = parms )

    # t = 1
    # p = ( 1 - t ) * x + t * res.barycenters
    # y = ( 1 - t ) * a + t * a / ( res.boundaries[ :, 1 ] - res.boundaries[ :, 0 ] )
    # plt.plot( p, y )
    # plt.plot( x, a )
    # plt.plot( x, b )

def animation_F2D( only_the_curve = 1 ):
    def func( x ):
        return .001 + 2 * gaussian( x, 0.3, 0.1 ) + 1 * gaussian( x, 0.7, 0.2 )
    n = 3000
    x = np.linspace( 0, 1, n, endpoint = True )
    a = np.array( list( map( func, x ) ) )
    a /= np.sum( a )

    d = np.linspace( x[ 0 ], x[ -1 ], 8 + 2, endpoint = True )[ 1: -1 ]
    res = usdot.from_p1_grid( d, x, a )

    if only_the_curve:
        t = 0
        for n in range( len( d ) ):
            xs = []
            ys = []
            bx = res.boundaries[ n ][ 0 ]
            ex = res.boundaries[ n ][ 1 ]
            for cx in np.linspace( bx, ex, endpoint = 1 ):
                ax = plt.gca()
                ax.set_xlim( [ 0, 1 ] )
                ax.set_ylim( [ 0, 2.5 ] )
                # affichage_dirac_comme_fonction( da, '--', 1, "red" )
                # affichage_dirac_comme_fonction( db, '-.', 1, "blue" )
                # affichage_dirac_comme_fonction( ( 1 - t ) * da + t * db, "-", 2, "black" )
                ys.append( func( cx ) * ( 1 - t ) + t * 1.5 )
                xs.append( ( 1 - t ) * cx + t * d[ n ] )
            ys.insert( 0, 0 )
            xs.insert( 0, xs[ 0 ] )
            ys.append( 0 )
            xs.append( xs[ -1 ] )
            # , color = "black"
            plt.plot( xs, ys, linewidth = 2 )
    
        plt.savefig( f'f2d.png', dpi = 300 )
        plt.show()


    eps = 1e-3
    x1 = list( np.linspace( 0, 1 - eps, 20 ) )
    x2 = list( np.linspace( 1 - eps, 0, 20 ) )
    for i, t in enumerate( x1 + x2 ):
        for n in range( len( d ) ):
            plt.plot( [ d[ n ], d[ n ] ], [ 0, 3 ], '--', linewidth = 1, color = "red" )
        plt.plot( x, np.array( list( map( func, x ) ) ), '--', linewidth = 1, color = "blue" )

        # # ax.set_ylim([ymin, ymax]
        for n in range( len( d ) ):
            xs = []
            ys = []
            bx = res.boundaries[ n ][ 0 ]
            ex = res.boundaries[ n ][ 1 ]
            for cx in np.linspace( bx, ex, endpoint = 1 ):
                ax = plt.gca()
                ax.set_xlim( [ 0, 1 ] )
                ax.set_ylim( [ 0, 2.5 ] )
                # affichage_dirac_comme_fonction( da, '--', 1, "red" )
                # affichage_dirac_comme_fonction( db, '-.', 1, "blue" )
                # affichage_dirac_comme_fonction( ( 1 - t ) * da + t * db, "-", 2, "black" )
                ys.append( func( cx ) * ( 1 - t ) + t * 1.5 )
                xs.append( ( 1 - t ) * cx + t * d[ n ] )
            if t:
                ys.insert( 0, 0 )
                xs.insert( 0, xs[ 0 ] )
                ys.append( 0 )
                xs.append( xs[ -1 ] )
                
            plt.plot( xs, ys, linewidth = 2, color = "black" )
    
        # plt.show()
        plt.savefig( f'f2f_{ 100 + i }.png', dpi = 300 )
        plt.cla()

def pd():
    from sdot import PowerDiagram
    import numpy as np

    # to get a nice picture
    np.random.seed( 357 )

    # we add some 2D seeds. By default, weights are equal to 1, leading to a voronoi diagram
    dp = np.random.random( [ 40, 2 ] )
    pd = PowerDiagram( dp )

    # when display context is not specified, PowerDiagram.plot uses pyplot
    pd.plot( plt )
    for p in dp:
        print( p )
        plt.plot( [ p[ 0 ] ], [ p[ 1 ] ], "+", color = "red" )
    # plt.savefig( f'f2d.png', dpi = 300 )
    plt.show()

def delaunay():
    from scipy.spatial import Delaunay
    import matplotlib.pyplot as plt
    import numpy as np

    points = np.array([[0, 0], [0, 1.], [1, 0], [1, 1]])
    for i in range( 12 ):
        tri = Delaunay( points )
        plt.triplot( points[:,0], points[:,1], tri.simplices )
        plt.triplot( points[:,0], points[:,1], tri.simplices )
        plt.plot( points[:,0], points[:,1], 'o' )

        plt.axis('off')

        plt.savefig( f'f2f_{ 100 + i }.png', dpi = 300, bbox_inches='tight' )
        plt.show()
        plt.cla()


        points = list( points )
        points.append( np.random.random( [ 2 ] ) )
        points = np.array( points )
        print( points )
        # plt.savefig( f'f2f_{ 100 + i }.png', dpi = 300 )

def cell_ctor():
    from sdot import Cell

    cell = Cell( ndim = 2 )

    cell.cut( [ -1,  0 ], 0 )
    cell.cut( [  0, -1 ], 0 )
    cell.cut( [ +1,  0 ], 1 )
    cell.cut( [  0, +1 ], 1 )

    c = np.array( [ 0.5, 0.5 ] )
    np.random.seed( 355 )

    ds = []
    for i in range( 10 ):
        ds.append( np.random.random( [ 2 ] ) )
    for i in range( 10 ):
        plt.axis( 'off' )
        ax = plt.gca()
        ax.set_xlim( [ -0.1, 1.1 ] )
        ax.set_ylim( [ -0.1, 1.1 ] )

        d = ds[ i ]
        e = d - c
        cell.cut( e, np.dot( e, ( d + c ) / 2 ) )
        cell.plot( ax )

        ax.plot( c[ 0 ], c[ 1 ], 'x', color = "blue" )
        for o in ds:
            ax.plot( o[ 0 ], o[ 1 ], 'x', color = "black" )
        ax.plot( d[ 0 ], d[ 1 ], 'o', linewidth = 2, color = "red" )

        plt.savefig( f'f2f_{ 100 + i }.png', dpi = 300, bbox_inches='tight' )
        plt.show()
        plt.cla()

def points_on_circle( C, R, n, a_off, na ):
    from pysdot.domain_types import ConvexPolyhedraAssembly
    from pysdot import OptimalTransport
    import numpy as np

    domain = ConvexPolyhedraAssembly()
    domain.add_convex_polyhedron([ [
        C[ 0 ] + R * np.cos( a ), # point X
        C[ 1 ] + R * np.sin( a ), # point Y
        np.cos( a ), # normal X
        np.sin( a )  # normal Y
    ] for a in a_off + np.linspace( 0, 2 * np.pi, na, endpoint=False ) ])

    # diracs
    P = C + R / 4 * np.random.random( [ n, 2 ] )
    for i in range( 5 ):
        ot = OptimalTransport(domain=domain)
        ot.set_positions( P )

        ot.adjust_weights()

        P = ot.get_centroids()

    # display
    # ot.display_vtk( "results/pd.vtk" )
    return P


def pot():
    import numpy as np
    import matplotlib.pylab as pl
    import ot, time
    import ot.plot

    n = 100  # nb samples

    mu_s = np.array([0, 0])
    cov_s = np.array([[1, 0], [0, 1]])

    mu_t = np.array([4, 4])
    cov_t = np.array([[1, -0.8], [-0.8, 1]])

    xs = ot.datasets.make_2D_samples_gauss(n, mu_s, cov_s)
    xt = ot.datasets.make_2D_samples_gauss(n, mu_t, cov_t)

    a, b = np.ones((n,)) / n, np.ones((n,)) / n  # uniform distribution on samples

    # loss matrix
    M = ot.dist(xs, xt)

    t0 = time.time()
    # G0 = ot.emd( a, b, M, numItermax = 1000000 )
    lambd = 1e-1

    G0 = ot.sinkhorn( a, b, M, lambd )
    t1 = time.time()

    print( t1 - t0 )
    ot.plot.plot2D_samples_mat(xs, xt, G0, c=[0.5, 0.5, 1])
    # pl.show()
    
def plot_D2D( n = 50 ):
    print( n )
    import ot, time
    import ot.plot

    np.random.seed( 357 )
    P0 = points_on_circle( [ 0.0, 0.0 ], 1.2, n, 0.0, 4 )
    P1 = points_on_circle( [ 3.5, 1.5 ], 1.0, n, 0.1, 5 )

    M = ot.dist( P0, P1 )
    a, b = np.ones((n,)) / n, np.ones((n,)) / n  # uniform distribution on samples
    t0 = time.time()
    # G0 = ot.emd( a, b, M, numItermax = 1000000 )
    G0 = ot.sinkhorn( a, b, M, 1e-1 )
    # lambd = 1e-1
    # G0 = ot.bregman.empirical_sinkhorn( P0, P1, lambd )
    t1 = time.time()

    # plt.gca().set_aspect('equal')
    # ot.plot.plot2D_samples_mat( P0, P1, G0, c=[0.5, 0.5, 1])
    # plt.plot( P0[ :, 0 ], P0[ :, 1 ], 'o', color = "red" )
    # plt.plot( P1[ :, 0 ], P1[ :, 1 ], 'o', color = "green" )
    # plt.axis( 'off' )    
    # plt.show()
    print( t1 - t0 )
    return t1 - t0

def plot_SD( n = 50 ):
    print( n )
    from pysdot.domain_types import ConvexPolyhedraAssembly
    from pysdot import OptimalTransport
    import ot, time
    import ot.plot

    np.random.seed( 357 )
    C0 = [ 0, 0 ]
    R0 = 1.1
    C1 = [ 3.5, 1.5 ]
    R1 = 1.0

    # P0 = points_on_circle( [ 0.0, 0.0 ], 1.2, n, 0.0, 4 )
    # P1 = points_on_circle( [ 3.5, 1.5 ], 1.0, n, 0.1, 5 )

    # P0 = points_on_circle( C, R, n, 0.0, 4 )
    # P1 = points_on_circle( [ 2.5, 1.2 ], 1.0, n, 0.1, 5 )
    P1 = points_on_circle( [ 0, 0 ], 0.8, n, 0.1, 5 )

    domain = ConvexPolyhedraAssembly()
    domain.add_convex_polyhedron([ [
        C0[ 0 ] + R0 * np.cos( a ), # point X
        C0[ 1 ] + R0 * np.sin( a ), # point Y
        np.cos( a ), # normal X
        np.sin( a )  # normal Y
    ] for a in np.linspace( 0, 2 * np.pi, 4, endpoint=False ) ])

    # diracs
    ot = OptimalTransport(domain=domain)
    ot.set_positions( P1 )
    ot.verbosity = 2

    t0 = time.time()
    ot.adjust_weights()
    t1 = time.time()

    # P0 = ot.get_centroids()

    # # plt.gca().set_aspect('equal')
    # # ot.plot.plot2D_samples_mat( P0, P1, G0, c=[0.5, 0.5, 1])
    # # plt.plot( P0[ :, 0 ], P0[ :, 1 ], 'o', color = "red" )
    # plt.plot( C1[ 0 ] + R1 * P1[ :, 0 ], C1[ 1 ] + R1 * P1[ :, 1 ], 'o', color = "green" )
    # for i in range( n ):
    #     plt.plot( [ C1[ 0 ] + R1 * P1[ i, 0 ], P0[ i, 0 ] ], [ C1[ 1 ] + R1 * P1[ i, 1 ], P0[ i, 1 ] ], color = "blue" )

    # # xs = [+0.338, -0.514, -0.313, 0.629, +0.338]
    # # ys = [-0.594, -0.324, +0.540, 0.300, -0.594]
    # xs = []
    # ys = []
    # for a in np.linspace( 0, 2 * np.pi, 4 + 1, endpoint=1 ) + np.pi / 4:
    #     xs.append( 1.5 * np.cos( a ) )
    #     ys.append( 1.5 * np.sin( a ) )
    # plt.plot( xs, ys, color = "black" )

    # # plt.axis( 'off' )
    # plt.show()
    return t1 - t0

def plot_scaling():
    xs = [ 8, 56 ]
    n = 15
    ys = [ 150 + 46 * n, 24 + 8 * n ]
    zs = [ 53 * ( n + 1 ), 48 * ( n + 1 ) ]

    plt.loglog( xs, ys, label = "SDOT" )
    plt.loglog( xs, zs, label = "Geogram" )
    plt.ylabel( "execution time (3D case)" )
    plt.xlabel( "# cores" )
    plt.legend()
    plt.show()

# def plot_laby( man = False, nb_diracs = 1000 ):
#     from pysdot.domain_types import ScaledImage
#     from pysdot import OptimalTransport
#     import scipy, cv2

#     np.random.seed( 357 )

#     img = cv2.imread( "labyrinthe.png" )[ :, :, 0 ] * 1.0
#     img /= np.mean(img)
#     old = img.copy()

#     for len in np.linspace( 800, 1, 10 ):
#         t = np.linspace( -1, 1, 2 * len )
#         x, y = np.meshgrid( t, t )
#         gau = np.exp( - ( x**2 + y**2 ) / len**2 )

#         img = scipy.signal.convolve2d( old, gau, fillvalue=1 )
#         img /= np.mean(img)

#         plt.imshow( img )

#         # # domain
#         # ot = OptimalTransport()
#         # # ot.set_positions( points_on_circle( [ 100, 500 ], 100, nb_diracs, 0.0, 20 ) )
#         # ot.set_positions( np.random.random( [ nb_diracs, 2 ] ) * np.array( [ 200, 600 ] ) )
#         # ot.set_domain( ScaledImage( [0, 0], [800, 600], img ) )
#         # ot.verbosity = 2
#         # ot.max_iter = 200000

#         # # display
#         # cpt = 0
#         # def cb( ot ):
#         #     nonlocal cpt
#         #     if cpt % 5 == 0:
#         #         ot.pd.display_vtk( f"results/pd_{ 100 + int( cpt / 5 ) }.vtk" )
#         #     print( cpt )
#         #     cpt += 1
#         # ot.adjust_weights( cb = cb, relax = 1e-1 )
def plot_laby( man = False, nb_diracs = 1000 ):
    from pysdot.domain_types import ScaledImage
    from pysdot import OptimalTransport
    import scipy, cv2

    np.random.seed( 357 )

    img = cv2.imread( "labyrinthe.png" )[ :, :, 0 ] * 1.0
    img /= np.mean(img)
    old = img.copy()

    # domain
    ot = OptimalTransport()
    ot.set_positions( np.random.random( [ nb_diracs, 2 ] ) * np.array( [ 200, 600 ] ) )
    ot.max_iter = 200000
    ot.verbosity = 2

    lens = list( np.linspace( 200, 1, 40 ) ) # + list( np.linspace( 30, 1, 10 ) )
    for cpt, len in enumerate( lens ):
        print( len )
        t = np.linspace( -1, 1, 2 * int( len ) )
        x, y = np.meshgrid( t, t )
        gau = np.exp( - ( x**2 + y**2 ) / len**2 )

        img = scipy.signal.convolve2d( old, gau, mode = 'same', fillvalue=1 )
        img /= np.mean(img)

        # plt.imshow( img )
        # plt.show()

        ot.set_domain( ScaledImage( [0, 0], [800, 600], img ) )

        # display
        ot.pd.display_vtk( f"results/pd_{ 100 + cpt }.vtk" )
        ot.adjust_weights()

# animation_F2F()
# animation_F2D()
# pd()
# delaunay()
# cell_ctor()
# pot()
# plot_D2D()
# plot_laby()
# plot_scaling()
def timings():
    import matplotlib
    # dt = [ plot_SD ( int( n ) ) for n in [ 1e2, 1e3, 1e4, 1e5, 1e6 ] ]
    bx = list( np.logspace( np.log10( 500 ), np.log10( 30000 ), 10, endpoint=True ) )
    sx = bx + [ 1e5 ]
    dx = bx
    
    # dt = [ plot_D2D( int( n ) ) for n in dx ] # [ 0.0020368099212646484, 0.006168842315673828, 0.02243781089782715, 0.6325011253356934, 2.4354560375213623, 9.787442922592163, 21.975932836532593 ]
    # st = [ plot_SD( int( n ) ) for n in sx ] # [ 0.0027718544006347656, 0.011106729507446289, 0.021039247512817383, 0.3116772174835205, 7.801603078842163 ] # , 280.72287487983704
    dt = [0.006140232086181641, 0.012012720108032227, 0.03782773017883301, 0.10511994361877441, 0.25022101402282715, 0.5880172252655029, 1.4581818580627441, 3.5607998371124268, 8.894384860992432, 26.059863090515137]
    st = [0.009979963302612305, 0.012831926345825195, 0.022614240646362305, 0.03346419334411621, 0.051442861557006836, 0.09084796905517578, 0.17048406600952148, 0.2967109680175781, 0.49010419845581055, 1.0261833667755127, 5.328494310379028]
    print( dt )    
    print( st )    

    font = {
        # 'family' : 'normal',
        # 'weight' : 'bold',
        'size'   : 15}

    matplotlib.rc('font', **font)

    plt.loglog( sx, st, label = "semi", linewidth = 2 )
    plt.loglog( dx, dt, label = "D2D", linewidth = 2 )
    plt.ylabel( "execution time (s)" )
    plt.xlabel( "# diracs" )
    plt.legend( loc = 'upper left' )
    plt.show()

def lloyd_patat( nb_diracs = 10000 ):
    from pysdot.domain_types import ScaledImage
    from pysdot import OptimalTransport
    import scipy, cv2

    np.random.seed( 357 )

    img = cv2.imread( "patatoide.png" )[ :, :, 0 ] * 1.0
    img /= np.mean( img )
    old = np.flip( img, axis = 0 )

    # plt.imshow( img )
    # plt.show()
    print( img.shape )
    shape = [ 474, 355 ]

    # domain
    ot = OptimalTransport()
    ot.set_positions( np.random.random( [ nb_diracs, 2 ] ) * np.array( shape ) )
    ot.max_iter = 200000
    ot.verbosity = 2

    lens = list( np.linspace( 100, 14, 10 ) ) + list( np.linspace( 14, 4, 10, endpoint=True ) )
    # lens = [ 4 ]
    for cpt, len in enumerate( lens ):
        print( len )
        t = np.linspace( -1, 1, 2 * int( len ) )
        x, y = np.meshgrid( t, t )
        gau = np.exp( - ( x**2 + y**2 ) / len**2 )

        img = scipy.signal.convolve2d( old, gau, mode = 'same', fillvalue=old[ 0, 0 ] )
        img /= np.mean(img)

        # plt.imshow( img )
        # plt.show()

        ot.set_domain( ScaledImage( [ 0, 0 ], shape, img ) )
        ot.adjust_weights()

        ot.pd.display_vtk( f"results/pd_{ 100 }.vtk" )

    # np.save( "weight", ot.get_weights() )
    # w = np.load( "weight.npy" )
    # ot.set_weights( w )
    # ot.adjust_weights()

    plt.style.use( 'dark_background' )
    for i in range( 8 ):
        print( i )
        c = ot.get_centroids()
        # ot.pd.display_vtk( f"results/pd_{ 100 + i }.vtk", centroids=True )
        # plt.scatter( c[ :, 0 ], c[ :, 1 ], s = 0.5, c="white" )
        # plt.axis( 'equal' )
        # plt.axis( 'off' )
        # plt.savefig( f'pat_{ 100 + i }.png', dpi = 600, bbox_inches='tight' )
        # plt.show()
        # plt.cla()

        tau = 1.0
        ot.set_positions( ( 1 - tau ) * ot.get_positions() + tau * ot.get_centroids() )
        ot.set_weights( np.ones( nb_diracs ) )
        ot.adjust_weights()

        ot.pd.display_vtk( f"results/pd_{ 100 + i }.vtk" )

    return ot.get_centroids()

def registration( nb_diracs = 500 ):
    from pysdot.domain_types import ScaledImage
    from pysdot import OptimalTransport
    import scipy
    import cv2

    np.random.seed( 357 )

    img = cv2.imread( "patatoide.png" )[ :, :, 0 ] * 1.0
    img = np.flip( img, axis = 0 )
    img = np.max( img ) - img
    img /= np.mean( img )

    shape = np.array( [ img.shape[ 1 ], img.shape[ 0 ] ] )

    # domain
    # ot = OptimalTransport()
    # ot.set_positions( np.random.random( [ nb_diracs, 2 ] ) * np.array( shape ) )
    # ot.max_iter = 200000
    # ot.verbosity = 2

    # lens = list( np.linspace( 100, 14, 10 ) ) + list( np.linspace( 14, 4, 10, endpoint=True ) )
    # old = img.copy()
    # for cpt, len in enumerate( lens ):
    #     print( len )
    #     t = np.linspace( -1, 1, 2 * int( len ) )
    #     x, y = np.meshgrid( t, t )
    #     gau = np.exp( - ( x**2 + y**2 ) / len**2 )

    #     img = scipy.signal.convolve2d( old, gau, mode = 'same', fillvalue=old[ 0, 0 ] )
    #     img /= np.mean( img )

    #     ot.set_domain( ScaledImage( [ 0, 0 ], shape, img ) )
    #     ot.adjust_weights()

    # for i in range( 8 ):
    #     tau = 1.0
    #     ot.set_positions( ( 1 - tau ) * ot.get_positions() + tau * ot.get_centroids() )
    #     ot.set_weights( np.ones( nb_diracs ) )
    #     ot.adjust_weights()

    # np.save( "patatoide_centroids", ot.get_centroids() )

    # plt.scatter( c[ :, 0 ], c[ :, 1 ], s = 2 )
    # plt.axis( 'equal' )
    # plt.axis( 'off' )
    # plt.savefig( f'pat_{ 100 + i }.png', dpi = 600, bbox_inches='tight' )
    # plt.show()
    # plt.cla()
    pts = np.load( "patatoide_centroids.npy" )
    
    ot = OptimalTransport()
    ot.set_domain( ScaledImage( [ 0, 0 ], shape, img ) )
    ot.set_positions( pts )
    ot.max_iter = 200000
    ot.verbosity = 2

    ot.adjust_weights()

    for i, t in enumerate( np.linspace( 0, np.pi, 40 ) ):
        # diracs = pts.copy()
        a = np.pi / 2 * np.sin( t )
        R = np.array( [ [ np.cos( a ), -np.sin( a ) ], [ np.sin( a ), np.cos( a ) ] ] )
        diracs = np.dot( pts - shape / 2, R.T ) + shape / 2

        ot.set_domain( ScaledImage( [ 0, 0 ], shape, img ) )
        ot.set_positions( diracs )
        ot.adjust_weights()

        if i % 1 == 0:
            n_diracs = diracs + np.array( [ 0, 0 ] )
            ctd = ot.get_centroids()
            for d, p in zip( n_diracs, ctd ):
                plt.plot( [ d[ 0 ], p[ 0 ] ], [ d[ 1 ], p[ 1 ] ], color = "blue", linewidth = .5 )
            plt.scatter( ctd[ :, 0 ], ctd[ :, 1 ], s = 4, color = "black" )
            plt.imshow( np.max( img ) - img, cmap = "gray" )
            # plt.axis( 'equal' )
            plt.axis( 'off' )
            ax = plt.gca()
            ax.set_xlim( [ -100, shape[ 0 ] + 100 ] )
            ax.set_ylim( [ -100, shape[ 1 ] + 100 ] )
            plt.savefig( f'pat_{ 100 + i }.png', dpi = 600, bbox_inches='tight' )
            # plt.show()
            plt.cla()

def wasserstein( a, b ):
    import usdot

    a = a.copy()
    b = b.copy()
    a -= min( a )
    b -= min( b )

    parms = usdot.OtParms()
    parms.verbosity = 2 

    # res = usdot.from_p1_grid( x, x, b, relative_mass_ratios = a )

    d = np.arange( len( a ) )
    res = usdot.from_p1_grid( d, d, b, relative_mass_ratios = a, ot_parms = parms )
    # usdot.plot( res )
    print( res.cost )
    return res.cost


def misfit():

    xs = np.linspace( -7, 7, 1000 )
    ys = np.sin( xs ) * np.exp( - xs**2 )

    err_l = []
    err_w = []
    ol = np.linspace( -4, 4, 101 )
    for o in ol:
        os = xs - o
        zs = np.sin( os ) * np.exp( - os**2 )
        err_l.append( np.sum( np.pow( zs - ys, 2 ) ) )
        err_w.append( wasserstein( zs, ys ) )

    for io, o in enumerate( ol ):
        if io:
            continue
        os = xs - o
        zs = np.sin( os ) * np.exp( - os**2 )
        fig = plt.figure(figsize=(5*2, 2*2))
        # fig.subplots_adjust(wspace=0.4, bottom=0.3)
        ax1 = fig.add_subplot( 121 )
        ax2 = fig.add_subplot( 122 )

        ax1.plot( xs, ys )
        ax1.plot( xs, zs )

        ax2.plot( ol, err_l, label = "L2" )
        ax2.plot( ol, err_w, label = "Wasserstein" )
        ax2.legend()
        
        ax2.scatter( [ ol[ io ] ], [ err_l[ io ] ], s = 6, c = "black" )
        ax2.scatter( [ ol[ io ] ], [ err_w[ io ] ], s = 6, c = "black" )

        # plt.savefig( f'pat_{ 100 + int( io / 10 ) }.png', dpi = 200, bbox_inches='tight' )
        plt.show()
        plt.cla()

def test_wasserstein():
    import usdot

    parms = usdot.OtParms()
    parms.verbosity = 2 

    res = usdot.from_p1_grid( [ 0.25, 0.75 ], [ 0, 1 ], [ 1, 1 ], ot_parms = parms )
    # 
    print( res.boundaries )
    print( res.cost**2 )
    print( res.cost )

# timings()
# lloyd_patat()

# registration()
# misfit()
test_wasserstein()
