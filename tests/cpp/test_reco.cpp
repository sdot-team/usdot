#include <cstdlib>
#include <tl/support/string/to_string.h>
#include <tl/support/containers/Vec.h>
#include <tl/support/operators/sum.h>
#include <tl/support/operators/sp.h>
#include <tl/support/P.h>

#include <usdot/utility/ThreadPool.h>
#include <usdot/utility/glot.h>
#include <fstream>

struct Transform {
    using TF = double;
    using Pt = Vec<TF,2>;
  
    Pt apply( const Pt &p ) const {
        using namespace std;
        return { 
            p[ 0 ] * cos( rotation ) - p[ 1 ] * sin( rotation ) + translation[ 0 ],
            p[ 0 ] * sin( rotation ) + p[ 1 ] * cos( rotation ) + translation[ 1 ]
        };
    }

    Pt apply_rev( Pt p ) const {
        using namespace std;
        p -= translation;
        return { 
            + p[ 0 ] * cos( rotation ) + p[ 1 ] * sin( rotation ),
            - p[ 0 ] * sin( rotation ) + p[ 1 ] * cos( rotation )
        };
    }

    Pt translation;
    TF rotation; // in radians
};

struct ProjectionData {
    using Ts = Vec<Transform>;
    using TF = double;
    using Pt = Vec<TF,2>;
 
    PI nb_angles() const {
        return transforms.size();
    }

    Pt rot90( const Pt &p ) const {
        return { -p[ 1 ], +p[ 0 ] };
    }

    void for_each_ray( PI num_angle, auto &&func ) const {
        for( PI num_pixel = 0; num_pixel < nb_pixels; ++num_pixel ) {
            Pt p0 = screen_l + ( screen_r - screen_l ) * ( num_pixel + TF( 0 ) ) / nb_pixels;
            Pt p1 = screen_l + ( screen_r - screen_l ) * ( num_pixel + TF( 1 ) ) / nb_pixels;
            Pt n0 = - rot90( p0 - emission_center );
            Pt n1 =   rot90( p1 - emission_center );
            if ( sp( n0, ( p1 - p0 ) / 2 ) > 0 ) {
                n0 = - n0;
                n1 = - n1;
            }
            func( num_pixel, n0, sp( n0, p0 ), n1, sp( n1, p1 ) );
        }
    }

    static Pt solve( TF m00, TF m01, TF m10, TF m11, TF v0, TF v1 ) {
        TF det = m00 * m11 - m01 * m10;
        if ( abs( det ) < 1e-10 )
            throw std::runtime_error( "Singular matrix in solve" );
        return {
            ( v0 * m11 - v1 * m01 ) / det,
            ( v1 * m00 - v0 * m10 ) / det
        };
    }

    TF pixel( const Pt &p ) const {
        Pt r = solve( 
            p[ 0 ] - emission_center[ 0 ], screen_l[ 0 ] - screen_r[ 0 ], 
            p[ 1 ] - emission_center[ 1 ], screen_l[ 1 ] - screen_r[ 1 ], 
            screen_l[ 0 ] - emission_center[ 0 ],
            screen_l[ 1 ] - emission_center[ 1 ]
        );
        return r[ 1 ] * nb_pixels;
    }

    void add_to_system( Vec<TF,6> &system, PI num_angle, TF pixel_pos ) {
        Pt t_emission_center = transforms[ num_angle ].apply_rev( emission_center );
        Pt t_screen_l = transforms[ num_angle ].apply_rev( screen_l );
        Pt t_screen_r = transforms[ num_angle ].apply_rev( screen_r );

        Pt prj = t_screen_l + ( t_screen_r - t_screen_l ) * pixel_pos / nb_pixels;
        Pt dir = rot90( prj - t_emission_center );
        TF spd = sp( dir, prj );
        system[ 0 ] += dir[ 0 ] * dir[ 0 ];
        system[ 1 ] += dir[ 0 ] * dir[ 1 ];
        system[ 2 ] += dir[ 1 ] * dir[ 0 ];
        system[ 3 ] += dir[ 1 ] * dir[ 1 ];
        system[ 4 ] += dir[ 0 ] * spd;
        system[ 5 ] += dir[ 1 ] * spd;
    }

    Pt  emission_center;
    PI  nb_pixels;
    Pt  screen_l;
    Pt  screen_r;

    Ts  transforms;
};

struct Polygon {
    using TF = double;
    using Pt = Vec<TF,2>;

    void cut( Pt normal, TF offset ) {
        if ( points.size() < 3 )
            return; // nothing to cut

        auto on_edge = [&]( const Pt &p0, const Pt &p1 ) {
            const TF s0 = sp( normal, p0 ) - offset;
            const TF s1 = sp( normal, p1 ) - offset;
            
            if ( s0 > 0 ) { // p0 is outside
                if ( s1 > 0 ) { // p1 is outside
                    // do nothing
                } else { // p1 is inside
                    const TF t = s0 / ( s0 - s1 );
                    new_points.push_back( p0 + t * ( p1 - p0 ) );
                    new_points.push_back( p1 );
                }
            } else { // p0 is inside
                if ( s1 > 0 ) { // p1 is outside
                    const TF t = s0 / ( s0 - s1 );
                    new_points.push_back( p0 + t * ( p1 - p0 ) );
                } else { // p1 is inside
                    new_points.push_back( p1 );
                }
            }
        };

        new_points.clear();
        for( PI i = 1; i < points.size(); ++i )
            on_edge( points[ i - 1 ], points[ i ] );
        on_edge( points.back(), points.front() );

        std::swap( points, new_points );
    }

    TF area() const {
        if ( points.size() < 3 )
            return 0;
        TF area = 0;
        for( PI i = 1; i < points.size(); ++i ) {
            area += points[ i - 1 ][ 0 ] * points[ i ][ 1 ];
            area -= points[ i - 1 ][ 1 ] * points[ i ][ 0 ];
        }
        area += points.back()[ 0 ] * points.front()[ 1 ];
        area -= points.back()[ 1 ] * points.front()[ 0 ];
        return area / 2;
    }

    Polygon transformed( const Transform &tr ) const {
        Polygon transformed_polygon;
        for( const Pt &p : points )
            transformed_polygon.points.push_back( tr.apply( p ) );
        return transformed_polygon;
    }

    Vec<Pt> new_points;
    Vec<Pt> points;
};

struct Sinogram {
    using TF = double;
    using Pt = Vec<TF,2>;
    using Values = Vec<Vec<TF>>;
    
    Sinogram( ProjectionData *pd ) : pd( pd ) {
        values.resize( pd->nb_angles(), FromSizeAndItemValue(), pd->nb_pixels, 0 );
    }

    void add_triangle( Pt p0, Pt p1, Pt p2, TF coeff = 1 ) {
        Polygon p_ref;
        p_ref.points << p0;
        p_ref.points << p1;
        p_ref.points << p2;

        for( PI n = 0; n < pd->nb_angles(); ++n ) {
            pd->for_each_ray( n, [&]( PI m, const Pt &n0, TF o0, const Pt &n1, TF o1 ) {
                Polygon p = p_ref.transformed( pd->transforms[ n ] );
                p.cut( n0, o0 );
                p.cut( n1, o1 );

                values[ n ][ m ] += coeff * p.area();
            } );
        }
    }

    void add_rectangle( Pt p0, Pt p1, TF coeff = 1 ) {
        add_triangle( p0, { p1[ 0 ], p0[ 1 ] }, p1, coeff );
        add_triangle( p0, p1, { p0[ 0 ], p1[ 1 ] }, coeff );
    }

    void plot( Str f = "glot.py" ) const {
        usdot::glot_stream( [&]( std::ofstream &fs ) {
            fs << "pyplot.imshow( " << to_string( values ) << " )\n";
        } );
    }

    Values          values;
    ProjectionData* pd;
};

struct Reconstructor {
    using TF = double;
    using Pt = Vec<TF,2>;

    Reconstructor( ProjectionData *pd, PI nb_points ) : pd( pd ) {
        for( PI i = 0; i < nb_points; ++i )
            pts << Pt{ TF( random() ) / RAND_MAX, TF( random() ) / RAND_MAX };
        // const TF L = std::sqrt( nb_points );
        // for( PI i = 0; i < L; ++i )
        //     for( PI j = 0; j < L; ++j )
        //         pts << Pt{ i / L, j / L };
    }

    static Vec<TF> barycenters( const Vec<TF> &dirac_positions, const Vec<TF> &density ) {
        // sorted_dirac_nums
        Vec<PI> sorted_dirac_nums( FromSizeAndFunctionOnIndex(), dirac_positions.size(), []( PI i ) { return i; } );
        std::sort( sorted_dirac_nums.begin(), sorted_dirac_nums.end(), [&]( PI a, PI b ) {
            return dirac_positions[ a ] < dirac_positions[ b ];
        } );

        //
        const TF coeff = dirac_positions.size() / sum( density );
        Vec<TF> res( FromSize(), dirac_positions.size() );
        TF acd = 0, acp = 1, b0 = 0, bx = 0;
        PI num_sorted_dirac = 0;
        for( PI i = 0; i <= density.size(); ++i ) {
            const TF inc = i < density.size() ? coeff * density[ i ] : std::numeric_limits<TF>::max() - acd;
            while ( num_sorted_dirac < dirac_positions.size() && acd + inc >= acp ) {
                const TF b1 = ( inc ? i + ( acp - acd ) / inc : i + 0.5 );
                const PI nd = sorted_dirac_nums[ num_sorted_dirac++ ];
                bx += inc * ( b1 - b0 ) * ( b0 + b1 ) / 2;
                res[ nd ] = bx;

                b0 = b1;
                bx = 0;

                acp += 1;
            }

            const TF b1 = i + 1;
            bx += inc * ( b1 - b0 ) * ( b0 + b1 ) / 2;
            b0 = b1;

            acd += inc;
        }

        return res;
    }

    void iteration( const Sinogram &sino, const TF random_pos = 0, const TF ratio = 0.95 ) {
        using Systems = Vec<Vec<TF,6>>;
        Vec<Systems> systems( FromSizeAndItemValue(), usdot::thread_pool.nb_threads(), FromSizeAndItemValue(), pts.size(), FromItemValue(), 0 );        
        usdot::thread_pool.execute( pd->nb_angles(), [&]( PI na, int num_thread ) {
            Vec<TF> projected_pts( FromSize(), pts.size() );
            for( PI i = 0; i < pts.size(); ++i )
                projected_pts[ i ] = pd->pixel( pd->transforms[ na ].apply( pts[ i ] ) );

            Vec<TF> bs = barycenters( projected_pts, sino.values[ na ] );
            for( PI nd = 0; nd < pts.size(); ++nd )
                pd->add_to_system( systems[ num_thread ][ nd ], na, ( 1 - ratio ) * projected_pts[ nd ] + ratio * bs[ nd ] );
        } );
        // for( PI na = 0; na < pd->nb_angles(); ++na ) {
        // }
        for( PI i = 1; i < usdot::thread_pool.nb_threads(); ++i )
            systems[ 0 ] += systems[ i ];

        for( PI nd = 0; nd < pts.size(); ++nd ) {
            pts[ nd ] = ProjectionData::solve( 
                systems[ 0 ][ nd ][ 0 ], systems[ 0 ][ nd ][ 1 ],
                systems[ 0 ][ nd ][ 2 ], systems[ 0 ][ nd ][ 3 ],
                systems[ 0 ][ nd ][ 4 ], systems[ 0 ][ nd ][ 5 ]
            ) + Pt{ random() * random_pos / RAND_MAX, random() * random_pos / RAND_MAX };
        }
    }

    void plot( PI n, Str f = "glot.py" ) const {
        auto xs = to_string( map_vec( pts, []( auto v ) { return v[ 0 ]; } ) );
        auto ys = to_string( map_vec( pts, []( auto v ) { return v[ 1 ]; } ) );
        std::ofstream fs( "glot.py" );
        fs << "from matplotlib import pyplot\n";
        fs << "pyplot.plot( " << xs << ", " << ys << ", '.' )\n";
        fs << "pyplot.axis( [ -1, 1, -1, 1 ] )\n";
        fs << "pyplot.axis( 'off' )\n";
        fs << "pyplot.savefig( 'pts_" << n << ".png', bbox_inches='tight', dpi=300 )\n";
        fs.close();
        system( "python glot.py" );
    }

    Vec<Pt>         pts;
    ProjectionData* pd;
};

int main() {
    using TF = double;
    using Pt = Vec<TF,2>;

    ProjectionData pd;
    pd.emission_center = { -100, 0 };
    pd.screen_l = { 1, -1 };
    pd.screen_r = { 1, +1 };
    pd.nb_pixels = 10000;
    for( PI n = 0, na = 100; n < na; ++n )
        pd.transforms.push_back_br( Pt{ 0, 0 }, 1 * M_PI * n / na );
        
    Sinogram sino( &pd );
    TF L = 0.30, l = 0.28;
    TF x = -0.35, y = 0.0;
    sino.add_rectangle( { x-L, y-L }, { x+L, y+L }, +1 );
    sino.add_rectangle( { x-l, y-l }, { x+l, y+l }, -1 );
    x = 0.35;
    sino.add_rectangle( { x-L, y-L }, { x+L, y+L }, +1 );
    sino.add_rectangle( { x-l, y-l }, { x+l, y+l }, -1 );
    // sino.plot();

    Reconstructor reco( &pd, 10000 );
    for( PI n = 0; n < 400; ++n ) {
        P( n );
        reco.iteration( sino, 1e-4, 0.95 ); // 1e-2 * ( n < 40 ) + 1e-3
        reco.plot( 1000 + n / 2 );
    }
    system( "convert pts_10* pts.gif" );
}

