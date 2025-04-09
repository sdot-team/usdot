#pragma once

#include "common_macros.h"
#include "common_types.h"
#include <utility>

namespace usdot {

/**
  Simple object pool


*/
class BumpPointerPool {
public:
    /* */       BumpPointerPool( const BumpPointerPool &that ) = delete;
    /* */       BumpPointerPool( BumpPointerPool &&that );
    /* */       BumpPointerPool();
    /* */      ~BumpPointerPool();

    void        operator=      ( const BumpPointerPool &that ) = delete;

    auto        allocate_max   ( PI min_size, PI max_size = 2048ul, PI alig = 8 ) -> std::pair<char *,PI>; ///< take end of the current buffer, up to `max_size`. If the current buffer does not contain `min_size`, get room from a new buffer.
    char*       allocate       ( PI size, PI alig );
    char*       allocate       ( PI size );

    T_TA T*     create         ( A &&...args );

    void        clear          ();
    void        free           ();

    static auto include_path   () { return "tl/support/BumpPointerPool.h"; }
    static auto type_name      () { return "TL_NAMESPACE::BumpPointerPool"; }

private:
    struct      Frame          { Frame *prev_frame; char *ending_ptr; char content[ 8 ]; };
    struct      Item           { virtual ~Item() {} Item *prev; };
    union       Exof           { char *cp; PI vp; };

    T_T struct  Inst : Item    { template<class... Args> Inst( Args &&...args ) : object{ std::forward<Args>( args )... } {} virtual ~Inst() {} T object; };

    Exof        current_ptr;   ///<
    char*       ending_ptr;    ///<
    Frame*      last_frame;    ///<
    Item*       last_item;     ///<
};

} // namespace usdot

#include "BumpPointerPool.cxx" // IWYU pragma: export
