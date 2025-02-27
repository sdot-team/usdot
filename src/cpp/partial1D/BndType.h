#pragma once

#include <tl/support/Displayer.h>

enum class BndType {
    Density,
    Cell,
    Ball,
    Inf,
};

inline void display( Displayer &ds, const BndType &value ) {
    switch ( value ) {
    case BndType::Density: ds << "Density"; return;
    case BndType::Cell   : ds << "Cell"   ; return;
    case BndType::Ball   : ds << "Ball"   ; return;
    case BndType::Inf    : ds << "Inf"    ; return;
    }
}
