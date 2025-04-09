#pragma once

// shortcuts for template<...>
#define     T_AUV                         template<class... A,class U,class V>
#define     T_Tij                         template<class T,int i,int j>
#define     T_UVi                         template<class U,class V,int i>
#define     T_ijs                         template<int i,int... j>
#define     T_TY                          template<template<typename...> class Y>
#define     T_TA                          template<class T,class... A>
#define     T_TI                          template<class T,std::size_t i>
#define     T_Ti                          template<class T,int i>
#define     T_iT                          template<int i,class T>
#define     T_UV                          template<class U,class V>
#define     T_iA                          template<int i,class... A>
#define     T_ij                          template<int i,int j>
#define     T_Ss                          template<CtStringValue... S>
#define     T_is                          template<int... i>
#define     T_S                           template<CtStringValue S>
#define     T_A                           template<class... A>
#define     T_T                           template<class T>
#define     T_U                           template<class U>
#define     T_I                           template<std::size_t i>
#define     T_i                           template<int i>

//
#define     SCPI                          static constexpr PI

// type handling
#define     VALUE_IN_DECAYED_TYPE_OF( v ) std::decay_t<decltype( v )>::value
#define     CT_DECAYED_TYPE_OF( v )       TL_NAMESPACE::CtType<std::decay_t<decltype( v )>>()
#define     STORAGE_TYPE_OF( v )          typename StorageTypeFor<decltype( v )>::value
#define     DECAYED_TYPE_OF( v )          std::decay_t<decltype( v )>
#define     IS_BASE_OF( A, V )            std::is_base_of_v<A,std::decay_t<V>>
#define     FORWARD( v )                  std::forward<decltype( v )>( v )

//
#define     NUA                           [[no_unique_address]]

//
#define     TL_OBJECT( NAME, INCL, ... )  static auto type_name() { return #NAME; } void for_each_attribute( auto &&f ) const { auto fea = [&f]( StrView names, const auto &...args ) { ( f( read_arg_name( names ), args ), ... ); }; fea( #__VA_ARGS__, ##__VA_ARGS__ ); }
