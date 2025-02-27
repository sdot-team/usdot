#pragma once

#include <tl/support/containers/Vec.h>
#include "Density.h"

/**
*/
template<class TF>
class PiecewiseAffineDensity : public Density<TF> {
public:
    /**/         PiecewiseAffineDensity( const Vec<TF> &xs, const Vec<TF> &ys ); ///< xs.size() must be equal to ys.size() + 1
       
    virtual auto convoluted            ( TF std ) const -> RcPtr<Density<TF>> override;
    virtual auto iterator              () const -> RcPtr<DensityIterator<TF>> override;
    virtual void display               ( Displayer &ds ) const override;
    virtual TF   mass                  () const override;

    Vec<TF>      xs;
    Vec<TF>      ys;
};

#include "PiecewiseAffineDensity.cxx"
