#include <tl/support/containers/Vec.h>
#include <tl/support/Displayer.h>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/LU>

/** */
template<class TF,int dim>
struct AffineTransformation {
    using Pt = Vec<TF,dim>;
    using TM = Eigen::Matrix<TF,dim+1,dim+1>;
    using TV = Eigen::Matrix<TF,dim+1,1>;

    AffineTransformation() {
        for( PI r = 0; r < dim + 1; ++r )
            for( PI c = 0; c < dim + 1; ++c )
                trans.coeffRef( r, c ) = ( r == c );
    }

    void display( Displayer &ds ) const {
        ds.start_array();
        for( PI r = 0; r < dim + 1; ++r )
            for( PI c = 0; c < dim + 1; ++c )
                ds << trans( r, c );
        ds.end_array();
    }

    static AffineTransformation rotation( TF angle ) {
        TF c = std::cos( angle );
        TF s = std::sin( angle );
        
        AffineTransformation res;
        res.trans.coeffRef( 0, 0 ) =  c; res.trans.coeffRef( 0, 1 ) =  s;
        res.trans.coeffRef( 1, 0 ) = -s; res.trans.coeffRef( 1, 1 ) =  c;
        return res;
    }

    static AffineTransformation translation( Pt tr ) {
        AffineTransformation res;
        for( PI d = 0; d < dim; ++d )
            res.trans.coeffRef( d, dim ) = tr[ d ];
        return res;
    }

    static AffineTransformation scale( TF s ) {
        AffineTransformation res;
        for( PI d = 0; d < dim; ++d )
            res.trans.coeffRef( d, d ) = s;
        return res;
    }

    static AffineTransformation linear( const auto &M ) {
        AffineTransformation res;
        for( PI r = 0; r < dim; ++r )
            for( PI c = 0; c < dim; ++c )
                res.trans.coeffRef( r, c ) = M( r, c );
        return res;
    }

    AffineTransformation operator*( const AffineTransformation &that ) const {
        AffineTransformation res;
        res.trans = that.trans * trans;
        return res;
    }

    AffineTransformation &operator*=( const AffineTransformation &that ) {
        *this = *this * that;
        return *this;
    }

    Pt apply( Pt p ) const {
        TV pv;
        for( PI r = 0; r < dim; ++r )
            pv[ r ] = p[ r ];
        pv[ dim ] = 1;

        TV resv = trans * pv;

        Pt res;
        for( PI r = 0; r < dim; ++r )
            res[ r ] = resv[ r ];
        return res;
    }

    TM trans;
};