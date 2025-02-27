#pragma once

#include <tl/support/memory/RcPtr.h>

/**
  Sub-interval of a density where `integral` and `value` are "easy" to compute (i.e. there's no need to find in which sub-interval x or [x0,x1] belong)

  We assume that iterator.x0 < iterator.x1 (strictly)
*/
template<class TF>
class DensityIterator : public WithRefCount {
public:
    virtual     ~DensityIterator() {}
    
    virtual bool move_backward  () = 0; ///< try to move backward. return false if not possible (i.e. if x0 is at the beginning of the support)
    virtual bool move_forward   () = 0; ///< try to move forward. return false if not possible (i.e. if x1 is at the end of the support)

    virtual void barycenter     ( TF &pint, TF &area, TF x0, TF x1 ) const = 0;
    virtual TF   integral       ( TF x0, TF x1 ) const = 0;
    virtual void display        ( Displayer &ds ) const = 0;
    virtual TF   value          ( TF pos ) const = 0;
     
    TF           x0;            ///<
    TF           x1;            ///<
};
