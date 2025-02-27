#pragma once

#include <tl/support/containers/Vec.h>
#include "Density.h"

/**
*/
template<class TF>
class SumOfDensities : public Density<TF> {
public:
    using        SubDensities   = Vec<RcPtr<Density<TF>>>;
 
    /**/         SumOfDensities ( const SubDensities &sub_densities = {} );
        
    virtual auto convoluted     ( TF std ) const -> RcPtr<Density<TF>> override;
    virtual auto iterator       () const -> RcPtr<DensityIterator<TF>> override;
    virtual void display        ( Displayer &ds ) const override;
    virtual TF   mass           () const override;

    SubDensities sub_densities; ///<
};

#include "SumOfDensities.cxx"
