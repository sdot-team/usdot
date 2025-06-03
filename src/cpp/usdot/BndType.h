#pragma once

#include "utility/common_types.h"

namespace usdot {

enum class BndType : PI8 {
    // Disk: 0
    // Void: 1
    // Eps: 2
    VoidToDisk = 1 + 4 * 0,
    DiskToDisk = 0 + 4 * 0,
    DiskToVoid = 0 + 4 * 1,
    DiskToEps  = 0 + 4 * 2,
    EpsToDisk  = 2 + 4 * 0,
    EpsToEps   = 2 + 4 * 2,
};

} // namespace usdot

#ifdef TL_DISPLAYER_IS_DEFINED
inline void display( Displayer &ds, const usdot::BndType &value ) {
    switch ( value ) {
    case usdot::BndType::VoidToDisk: ds << "VoidToDisk"; return;
    case usdot::BndType::DiskToDisk: ds << "DiskToDisk"; return;
    case usdot::BndType::DiskToVoid: ds << "DiskToVoid"; return;
    case usdot::BndType::DiskToEps : ds << "DiskToEps" ; return;
    case usdot::BndType::EpsToDisk : ds << "EpsToDisk" ; return;
    case usdot::BndType::EpsToEps  : ds << "EpsToEps"  ; return;
    }
}
#endif
