#pragma once

#include "SumOfDensities.h"
#include <limits>

namespace usdot {

template<class TF>
class SumOfDensitiesDensityIterator : public DensityIterator<TF> {
public:
    virtual void barycenter( TF &pint, TF &area, TF l0, TF l1 ) const override {
        for( const auto &sd : iterators ) {
            if ( sd.second->x0 < l1 && sd.second->x1 > l0 ) {
                TF lpint = 0, larea = 0;
                sd.second->barycenter( lpint, larea, max( l0, sd.second->x0 ), min( l1, sd.second->x1 ) );
                pint += sd.first * lpint;
                area += sd.first * larea;
            }
        }
    }
    
    virtual TF integral( TF l0, TF l1 ) const override {
        TF res = 0;
        for( const auto &sd : iterators )
            if ( sd.second->x0 < l1 && sd.second->x1 > l0 )
                res += sd.first * sd.second->integral( max( l0, sd.second->x0 ), min( l1, sd.second->x1 ) );
        return res;
    }

    virtual TF value( TF pos ) const override {
        TF res = 0;
        for( const auto &sd : iterators )
            if ( sd.second->x0 <= pos && sd.second->x1 > pos )
                res += sd.first * sd.second->value( pos );
        return res;
    }

    virtual void display( Displayer &ds ) const override {
        DS_OBJECT( SumOfDensitiesDensityIterator, this->x0, this->x1, iterators );
    }
    
    virtual bool move_backward() override {
        using namespace std;

        // try find the previous "event" (i.e. the previous bnd of interval < this->x0)
        TF new_x0 = numeric_limits<TF>::lowest();
        for( const auto &sd : iterators ) {
            while ( sd.second->x0 >= this->x0 )
                if ( ! sd.second->move_backward() )
                    break;
            if ( sd.second->x0 < this->x0 && sd.second->x0 > new_x0 )
                new_x0 = sd.second->x0;
            if ( sd.second->x1 < this->x0 && sd.second->x1 > new_x0 )
                new_x0 = sd.second->x1;
        }

        // no new event -> impossible to move forward
        if ( new_x0 == this->x0 )
            return false;

        // else, use the new interval
        this->x1 = this->x0;
        this->x0 = new_x0;
        return true;
    }

    virtual bool move_forward() override {
        using namespace std;

        // try find the next "event" (i.e. the next bnd of interval > this->x1)
        TF new_x1 = numeric_limits<TF>::max();
        for( const auto &sd : iterators ) {
            while ( sd.second->x1 <= this->x1 )
                if ( ! sd.second->move_forward() )
                    break;
            if ( sd.second->x0 > this->x1 && sd.second->x0 < new_x1 )
                new_x1 = sd.second->x0;
            if ( sd.second->x1 > this->x1 && sd.second->x1 < new_x1 )
                new_x1 = sd.second->x1;
        }

        // no new event -> impossible to move forward
        if ( new_x1 == this->x1 )
            return false;

        // else, use the new interval
        this->x0 = this->x1;
        this->x1 = new_x1;
        return true;
    }

    Vec<std::pair<TF,RcPtr<DensityIterator<TF>>>> iterators; ///< beware: these iterators may be "valid" but without any intersection with *this
};

#define DTP template<class TF>
#define UTP SumOfDensities<TF>

DTP UTP::SumOfDensities( const SubDensities &sub_densities ) : sub_densities( sub_densities.filtered( []( const auto &sd ) { return sd.first; } ) ) {
}

DTP CdfApproximation<TF> UTP::cdf_approximation( TF epsilon ) const {
    if ( sub_densities.size() == 1 ) {
        CdfApproximation<TF> res = sub_densities[ 0 ].second->cdf_approximation( epsilon );
        for( TF &v : res.ys )
            v *= sub_densities[ 0 ].first;
        for( TF &v : res.zs )
            v *= sub_densities[ 0 ].first;
        return res;
    }
    throw std::runtime_error( "TODO (sum of CdfApproximations)" );
}

DTP RcPtr<Density<TF>> UTP::convoluted( RcPtr<Convolution<TF>> convolution ) const {
    SubDensities ns;
    for( const auto &sd : sub_densities )
        ns.push_back( sd.first, sd.second->convoluted( convolution ) );
    return new SumOfDensities<TF>( ns );
}

DTP RcPtr<DensityIterator<TF>> UTP::iterator() const {
    using namespace std;

    // try to make a SumOfDensitiesDensityIterator
    auto res = RcPtr<SumOfDensitiesDensityIterator<TF>>::New();

    // get valid iterators + x0 (first event)
    res->x0 = numeric_limits<TF>::max();
    for( const auto &sd : sub_densities ) {
        if ( sd.first == 0 )
            continue;
        if ( RcPtr<DensityIterator<TF>> iterator = sd.second->iterator() ) {
            res->iterators.push_back( sd.first, iterator );
            res->x0 = min( res->x0, iterator->x0 );
        }
    }

    // no valid sub-iterator ?
    if ( res->iterators.empty() )
        return {};

    // 
    if ( res->iterators.size() == 1 && res->iterators[ 0 ].first == 1 )
        return res->iterators[ 0 ].second;

    // get x1
    res->x1 = numeric_limits<TF>::max();
    for( const auto &sd : res->iterators ) {
        if ( sd.second->x0 > res->x0 )
            res->x1 = min( res->x1, sd.second->x0 );
        if ( sd.second->x1 > res->x0 )
           res->x1 = min( res->x1, sd.second->x1 );
    }

    return res;
}

DTP void UTP::display( Displayer &ds ) const {
    DS_OBJECT( SumOfDensities, sub_densities );
}

DTP TF UTP::min_x( TF eps ) const {
    using namespace std;
    TF res = sub_densities[ 0 ].second->min_x( eps );
    for( PI i = 1; i < sub_densities.size(); ++i )
        res = min( res, sub_densities[ i ].second->min_x( eps ) );
    return res;
}

DTP TF UTP::max_x( TF eps ) const {
    using namespace std;
    TF res = sub_densities[ 0 ].second->max_x( eps );
    for( PI i = 1; i < sub_densities.size(); ++i )
        res = max( res, sub_densities[ i ].second->max_x( eps ) );
    return res;
}

DTP TF UTP::mass() const {
    TF res = 0;
    for( const auto &sd : sub_densities )
        res += sd.first * sd.second->mass();
    return res;
}

} // namespace usdot

#undef DTP
#undef UTP
