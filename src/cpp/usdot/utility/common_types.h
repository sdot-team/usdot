#pragma once

//#include <string_view>
#include <cstdint>
#include <string>

namespace usdot {

// number types
using PI64    = std::conditional< sizeof( unsigned long ) == 8, unsigned long, unsigned long long >::type;
using PI32    = std::uint32_t;
using PI32    = std::uint32_t;
using PI16    = std::uint16_t;
using PI8     = std::uint8_t;

using SI64    = std::conditional< sizeof( signed long ) == 8, signed long, signed long long >::type;
using SI32    = std::int32_t;
using SI16    = std::int16_t;
using SI8     = std::int8_t;

using Bool    = bool;

#if __cplusplus >= 202002L
using Byte    = std::byte;
#else
enum class byte : unsigned char {};
#endif

using PI      = std::conditional< sizeof( void * ) == 8, PI64, PI32 >::type;
using SI      = std::conditional< sizeof( void * ) == 8, SI64, SI32 >::type;

using FP80    = long double;
using FP64    = double;
using FP32    = float;

// common std types
//using StrView = std::string_view;
using Str     = std::string;

} // namespace usdot
