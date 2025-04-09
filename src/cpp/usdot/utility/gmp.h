#pragma once

#include <boost/multiprecision/cpp_bin_float.hpp>
#include <tl/support/Displayer.h>

T_T void display( Displayer &ds, const boost::multiprecision::number<T> &n ) { ds << n.str(); }
