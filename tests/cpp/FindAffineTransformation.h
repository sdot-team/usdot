#include "../../src/cpp/partial1D/Solver.h"

#include <tl/support/operators/norm_2.h>
#include <tl/support/operators/sp.h>
#include <chrono>

#include "AffineTransformation.h"

#include <eigen3/Eigen/src/SVD/JacobiSVD.h>
#include <eigen3/Eigen/Dense>

/** */
template<int dim>
struct FindAffineTransformation {
    using TF = double;
    using Pt = Vec<double,dim>;
    using Tr = AffineTransformation<dim>;

    using TM = Eigen::Matrix<TF,dim,dim>;
    using TV = Eigen::Matrix<TF,dim,1>;

    FindAffineTransformation() {
    }

    void display_points( auto &img, const Tr &glob_trans ) {
        unsigned char black[] = { 0, 0, 0 };
        unsigned char blue[] = { 0, 0, 255 };
        unsigned char red[] = { 255, 0, 0 };
        for( const Pt &pt : target ) {
            Pt fp = glob_trans.apply( pt );
            img.draw_point( fp[ 0 ], fp[ 1 ], black );
        }
        for( const Pt &pt : diracs ) {
            Pt fp = glob_trans.apply( dirac_trans.apply( pt ) );
            img.draw_circle( fp[ 0 ], fp[ 1 ], 2, red );
        }
        for( PI n = 0; n < new_diracs.size(); ++n ) {
            Pt op = glob_trans.apply( dirac_trans.apply( diracs[ n ] ) );
            Pt fp = glob_trans.apply( new_diracs[ n ] );
            // img.draw_circle( fp[ 0 ], fp[ 1 ], 2, blue );
            img.draw_line( op[ 0 ], op[ 1 ], fp[ 0 ], fp[ 1 ], blue, 1 );
        }
    }

    void display_diracs_vtk( VtkOutput &vo ) {
        for( const Pt &pt : diracs ) {
            VtkOutput::Pt p = dirac_trans.apply( pt );
            vo.add_point( p );
        }
    }

    void display_new_diracs_vtk( VtkOutput &vo ) {
        for( PI n = 0; n < diracs.size(); ++n ) {
            VtkOutput::Pt p = dirac_trans.apply( diracs[ n ] );
            VtkOutput::Pt q = new_diracs[ n ];
            vo.add_line( Vec<Pt>{ p, q } );
        }
    }

    void display_target_vtk( VtkOutput &vo ) {
        for( const Pt &pt : target ) {
            VtkOutput::Pt p = pt;
            vo.add_point( p );
        }
    }

    std::pair<Vec<TF>,Vec<TF>> make_density( Pt proj_dir ) {
        TF max_pos = std::numeric_limits<TF>::lowest();
        TF min_pos = std::numeric_limits<TF>::max();

        for( const Pt &p : target ) {
            TF v = sp( p, proj_dir );
            max_pos = max( max_pos, v );
            min_pos = min( min_pos, v );
        }

        Vec<TF> x{ FromSizeAndFunctionOnIndex(), nb_bins, [&]( PI ind ) { return min_pos + ( max_pos - min_pos ) * ind / ( nb_bins - 1 ); } };
        Vec<TF> y{ FromSizeAndItemValue(), nb_bins, 0 };
        for( const Pt &p : target ) {
            TF v = sp( p, proj_dir );
            
            TF pos = ( v - min_pos ) * ( nb_bins - 1 ) / ( max_pos - min_pos );
            PI ind = min( PI( pos ), nb_bins - 2 );
            // TF fra = pos - ind;

            if ( ind + 0 < y.size() ) y[ ind + 0 ] += 1; // ( 1 - fra );
            // if ( ind + 1  < y.size() ) y[ ind + 1 ] += fra;
        }
        // y.pop_back();

        return { x, y };
    }

    /// return translation to apply to each dirac AFTER application of `dirac_trans`
    Vec<TF> delta_for_dir( Pt proj_dir ) {
        // projected_density
        std::pair<Vec<TF>,Vec<TF>> projected_density = make_density( proj_dir );

        // projected_diracs
        Vec<TF> projected_diracs( FromReservationSize(), diracs.size() );
        for( PI n = 0; n < diracs.size(); ++n )
            projected_diracs << sp( new_diracs[ n ], proj_dir );

        auto t0 = std::chrono::high_resolution_clock::now();
        Solver s( projected_diracs, projected_density.first, projected_density.second, 1e-3, mass_ratio );
        s.solve();
        auto t1 = std::chrono::high_resolution_clock::now();
        PE( std::chrono::duration<double>{ t1 - t0 } );

        // optimal transport
        projected_density.second.pop_back();
        // auto t0 = std::chrono::high_resolution_clock::now();
        InitialSolution solver = make_initial_solution( projected_diracs, projected_density.first, projected_density.second, mass_ratio );
        // auto t1 = std::chrono::high_resolution_clock::now();
        // PE( std::chrono::duration<double>{ t1 - t0 } );

        // delta
        return solver.barycenters() - projected_diracs;
    }

    void compute_with_rand_dir() {
        // proj_dir
        TF u = rand() * 1.0 / RAND_MAX;
        TF v = rand() * 1.0 / RAND_MAX;
        TF theta = 2 * M_PI * u;
        TF phi = acos( 2 * v - 1 );
        TF x = sin( phi ) * cos( theta );
        TF y = sin( phi ) * sin( theta );
        TF z = cos( phi );

        Pt proj_dir{ x, y, z };

        // projected_density
        std::pair<Vec<TF>,Vec<TF>> projected_density = make_density( proj_dir );

        // projected_diracs
        Vec<TF> projected_diracs( FromReservationSize(), diracs.size() );
        for( PI n = 0; n < diracs.size(); ++n )
            projected_diracs << sp( new_diracs[ n ], proj_dir );

        // optimal transport
        InitialSolution solver = make_initial_solution( projected_diracs, projected_density.first, projected_density.second, mass_ratio );

        // delta
        Vec<TF> delta = solver.barycenters() - projected_diracs;

        //
    }

    Vec<Pt> get_dirs() {
        // TF a = num_dir * std::numbers::pi / nb_dirs;
        Vec<Pt> res;

        if ( dim == 2 ) {
            TF o = rand() * M_PI / RAND_MAX;
            for( PI n = 0; n < nb_dirs; ++n ) {
                TF a = o + n * M_PI / nb_dirs;
                res << Pt{ cos( a ), sin( a ) };
            }
            return res;
        }


        // for( PI n = 0; n < dim; ++n ) {
        //     Pt t;
        //     for( PI d = 0; d < dim; ++d )
        //         t[ d ] = ( d == n );
        //     res << t;
        // }

        for( PI n = 0; n < nb_dirs; ++n ) {
            TF u = rand() * 1.0 / RAND_MAX;
            TF v = rand() * 1.0 / RAND_MAX;
            TF theta = 2 * M_PI * u;
            TF phi = acos( 2 * v - 1 );
            TF x = sin( phi ) * cos( theta );
            TF y = sin( phi ) * sin( theta );
            TF z = cos( phi );

            res << Pt{ x, y, z };
        }

        return res;
    }

    void compute_new_diracs( PI nb_iters = 1 ) {
        new_diracs.resize( diracs.size() );
        for( PI n = 0; n < diracs.size(); ++n )
            new_diracs[ n ] = dirac_trans.apply( diracs[ n ] );

        for( PI num_iter = 0; num_iter < nb_iters; ++num_iter ) {
            Vec<Pt,dim> M;
            for( PI r = 0; r < dim; ++r )
                for( PI c = 0; c < dim; ++c )
                    M[ r ][ c ] = 0;

            Vec<Pt> V( FromSize(), diracs.size() );
            for( auto &v : V )
                for( PI r = 0; r < dim; ++r )
                    v[ r ] = 0;

            Vec<Pt> dirs = get_dirs();

            for( PI num_dir = 0; num_dir < dirs.size(); ++num_dir ) {
                Pt proj_dir = dirs[ num_dir ];

                for( PI r = 0; r < dim; ++r )
                    for( PI c = 0; c < dim; ++c )
                        M[ r ][ c ] += proj_dir[ r ] * proj_dir[ c ];

                Vec<TF> de = delta_for_dir( proj_dir );

                for( PI num_dirac = 0; num_dirac < diracs.size(); ++num_dirac )
                    for( PI r = 0; r < dim; ++r )
                        V[ num_dirac ][ r ] += proj_dir[ r ] * de[ num_dirac ];
            }

            // make matrix        
            TM eM;
            for( PI r = 0; r < dim; ++r )
                for( PI c = 0; c < dim; ++c )
                    eM.coeffRef( r, c ) = M[ r ][ c ];
            Eigen::FullPivLU<TM> lu( eM );

            // solve
            for( PI num_dirac = 0; num_dirac < diracs.size(); ++num_dirac ) {
                TV eV;
                for( PI r = 0; r < dim; ++r )
                    eV[ r ] = V[ num_dirac ][ r ];

                new_diracs[ num_dirac ] += 0.5 * Pt( lu.solve( eV ) );
            }
        }
    }

    Pt rotation_center() {
        using TM = Eigen::Matrix<TF,dim,dim>;
        using TV = Eigen::Matrix<TF,dim,1>;
        TM Mr = TM::Zero();
        TV Vr = TV::Zero();
        // min ( dir[ 0 ] * ( x[ 0 ] - odi[ 0 ] ) + dir[ 1 ] * ( x[ 1 ] - odi[ 1 ] ) )^2
        // dx0 -> dir[ 0 ] * ( dir[ 0 ] * ( x[ 0 ] - odi[ 0 ] ) + dir[ 1 ] * ( x[ 1 ] - odi[ 1 ] ) )
        // dx1 -> dir[ 1 ] * ( dir[ 0 ] * ( x[ 0 ] - odi[ 0 ] ) + dir[ 1 ] * ( x[ 1 ] - odi[ 1 ] ) )
        // Tr trans = Tr::rotation( 20 * M_PI / 180 ); //  * Tr::translation( { 1, 0 } )
        for( PI num_dirac = 0; num_dirac < diracs.size(); ++num_dirac ) {
            Pt o_d = dirac_trans.apply( diracs[ num_dirac ] );
            Pt n_d = new_diracs[ num_dirac ];
            Pt dir = n_d - o_d;
            Pt mid = ( n_d + o_d ) / 2;
            for( PI r = 0; r < dim; ++r )
                for( PI c = 0; c < dim; ++c )
                    Mr.coeffRef( r, c ) += dir[ r ] * dir[ c ];
            for( PI r = 0; r < dim; ++r )
                Vr.coeffRef( r ) += dir[ r ] * sp( dir, mid );
        }

        Eigen::FullPivLU<TM> lu( Mr );
        return lu.solve( Vr );
    }

    TF rotation_angle( const Pt &rotation_center ) {
        auto cross_prod = []( Pt a, Pt b ) {
            return a[ 0 ] * b[ 1 ] - a[ 1 ] * b[ 0 ];
        };

        TF rot_sum = 0;
        TF rot_coe = 0;
        for( PI num_dirac = 0; num_dirac < diracs.size(); ++num_dirac ) {
            Pt o = dirac_trans.apply( diracs[ num_dirac ] ) - rotation_center;
            Pt n = new_diracs[ num_dirac ] - rotation_center;
            TF no = norm_2( o ); 
            TF nn = norm_2( n ); 
            if ( nn == 0 || no == 0 )
                continue;
            o /= no;
            n /= nn;

            TF ang = asin( cross_prod( n, o ) ); 
            TF coe = no;

            rot_sum += coe * ang;
            rot_coe += coe;
        }

        if ( rot_coe )
             rot_sum /= rot_coe;

        return rot_sum;
    }

    Pt translation() {
        Pt translation( FromItemValue(), 0 );
        for( PI num_dirac = 0; num_dirac < diracs.size(); ++num_dirac )
            translation += new_diracs[ num_dirac ] - dirac_trans.apply( diracs[ num_dirac ] );
        translation /= diracs.size();
        return translation;
    }

    TF iteration_SVD() {
        // centers
        Pt center_old( FromItemValue(), 0 );
        Pt center_new( FromItemValue(), 0 );
        for( PI num_dirac = 0; num_dirac < diracs.size(); ++num_dirac ) {
            Pt old_pos = dirac_trans.apply( diracs[ num_dirac ] );
            Pt new_pos = new_diracs[ num_dirac ];
            center_old += old_pos;
            center_new += new_pos;
        }
        center_old /= diracs.size();
        center_new /= diracs.size();

        // covariance
        TM cov = TM::Zero();
        for( PI num_dirac = 0; num_dirac < diracs.size(); ++num_dirac ) {
            Pt old_pos = dirac_trans.apply( diracs[ num_dirac ] );
            Pt new_pos = new_diracs[ num_dirac ];
            for( PI r = 0; r < dim; ++r )
                for( PI c = 0; c < dim; ++c )
                    cov.coeffRef( r, c ) += ( new_pos[ r ] - center_new[ r ] ) * ( old_pos[ c ] - center_old[ c ] );
        }
        Eigen::JacobiSVD<TM> svd( cov, Eigen::ComputeFullU | Eigen::ComputeFullV );
        // std::cout << svd.singularValues() << std::endl;
        // std::cout << svd.matrixU() << std::endl;
        // std::cout << svd.matrixV() << std::endl;

        TM orth = svd.matrixU() * svd.matrixV().transpose();
        // std::cout << orth << std::endl;
        // std::cout << "det " << orth.determinant() << std::endl;

        TM rot = svd.matrixU() * svd.matrixV().transpose();
        // PE( ( rot - TM::Identity() ).norm() );
        double err = 0;
        err += pow( 1 - rot( 0 ), 2 );
        err += pow(     rot( 1 ), 2 );
        err += pow(     rot( 2 ), 2 );
        err += pow(     rot( 3 ), 2 );
        err += pow( 1 - rot( 4 ), 2 );
        err += pow(     rot( 5 ), 2 );
        err += pow(     rot( 6 ), 2 );
        err += pow(     rot( 7 ), 2 );
        err += pow( 1 - rot( 8 ), 2 );
        std::cout << std::sqrt( err ) << std::endl;

        dirac_trans *= Tr::translation( - center_old ) * Tr::linear( rot ) * Tr::translation( center_new );

        return std::sqrt( err );
    }

    void iteration_calc() {
        TF sa = 0;
        for( PI num_iter = 0; num_iter < 1; ++num_iter ) {
            // rotation
            Pt C = rotation_center();
            TF a = rotation_angle( C ); 
            dirac_trans *= Tr::translation( - C ) * Tr::rotation( a ) * Tr::translation( C );

            // translation
            Pt T = translation();
            dirac_trans *= Tr::translation( T );

            sa += a;
        }

        PE( sa );
    }

    void iteration_mat() {
        for( PI num_iter = 0; num_iter < 1; ++num_iter ) {
            using MT = Eigen::Matrix<TF,dim*(dim+1),dim*(dim+1)>;
            using VT = Eigen::Matrix<TF,dim*(dim+1),1>;
            MT M = MT::Zero();
            VT V = VT::Zero();

            // 2D case:
            //  ( U_00 * O_0 + U_01 * O_1 + U_02 - N_0 ) ^ 2 +
            //  ( U_10 * O_0 + U_11 * O_1 + U_12 - N_1 ) ^ 2
            // 
            for( PI num_dirac = 0; num_dirac < diracs.size(); ++num_dirac ) {
                Pt O = dirac_trans.apply( diracs[ num_dirac ] );
                Pt N = new_diracs[ num_dirac ];

                Vec<TF,dim+1> sf;
                for( PI r = 0; r < dim; ++r )
                    sf[ r ] = O[ r ];
                sf[ dim ] = 1;

                for( PI n = 0; n < dim; ++n ) {                    
                    for( PI r = 0; r <= dim; ++r ) {                    
                        for( PI c = 0; c <= dim; ++c )
                            M.coeffRef( r + n * ( dim + 1 ), c + n * ( dim + 1 ) ) += sf[ r ] * sf[ c ];
                        V.coeffRef( r + n * ( dim + 1 ) ) += sf[ r ] * N[ n ];
                    }
                }
            }

            // std::cout << M << std::endl;
            // std::cout << V << std::endl;
            //   1.03519 0.0231744 -1.03137
            //    0.1824 0.96904 -0.109215
    
            Eigen::FullPivLU<MT> lu( M );
            VT S = lu.solve( V );

            Tr ltrans;
            for( PI r = 0; r < dim; ++r )          
                for( PI c = 0; c <= dim; ++c )
                    ltrans.trans.coeffRef( r, c ) = S[ r * ( dim + 1 ) + c ];
            P( ltrans );

            dirac_trans *= ltrans;
        }
    }

    void iteration() {
        iteration_calc();
        // iteration_mat();

        // rotate new dirs
        // Tr ltrans = Tr::translation( - center ) * Tr::rotation( rot_sum ) * Tr::translation( center );
        // for( PI num_dirac = 0; num_dirac < diracs.size(); ++num_dirac )
        //     new_dirs[ num_dirac ] = dirac_trans.apply( diracs[ num_dirac ] ) + new_dirs[ num_dirac ] - ltrans.apply( dirac_trans.apply( diracs[ num_dirac ] ) );
        // ltrans *= Tr::translation( translation );
    }

    TF      mass_ratio = 1.0;
    PI      nb_bins = 100;
    PI      nb_dirs = 200;
    Vec<Pt> diracs; ///< 
    Vec<Pt> target; ///< large number of point

    Tr      dirac_trans; ///<
    Vec<Pt> new_diracs; ///<
};  