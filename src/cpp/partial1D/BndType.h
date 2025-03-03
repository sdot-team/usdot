#pragma once

#include <tl/support/Displayer.h>

namespace usdot {

enum class BndType {
    Density,
    Cell,
    Ball,
    Inf,
};

} // namespace usdot

inline void display( Displayer &ds, const usdot::BndType &value ) {
    switch ( value ) {
    case usdot::BndType::Density: ds << "Density"; return;
    case usdot::BndType::Cell   : ds << "Cell"   ; return;
    case usdot::BndType::Ball   : ds << "Ball"   ; return;
    case usdot::BndType::Inf    : ds << "Inf"    ; return;
    }
}
