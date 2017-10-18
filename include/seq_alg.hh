#ifndef IVANP_SEQ_ALG_HH
#define IVANP_SEQ_ALG_HH

#include <utility>
#include "meta.hh"

namespace ivanp { namespace seq {

template <typename... SS> struct join;

template <typename T, size_t... I>
struct join<std::integer_sequence<T,I...>> {
  using type = std::integer_sequence<T,I...>;
};

template <typename T, size_t... I1, size_t... I2>
struct join<std::integer_sequence<T,I1...>,std::integer_sequence<T,I2...>> {
  using type = std::integer_sequence<T,I1...,I2...>;
};

template <typename S1, typename S2, typename... SS>
struct join<S1,S2,SS...>: join<typename join<S1,S2>::type,SS...> { };

template <typename... SS> using join_t = typename join<SS...>::type;

template <typename Seq> struct head;
template <typename T>
struct head<std::integer_sequence<T>>: maybe<nothing> { };
template <typename T, T Head, T... I>
struct head<std::integer_sequence<T,Head,I...>>
: maybe<just<std::integral_constant<T,Head>>> { };
template <typename Seq> using head_t = typename head<Seq>::type;

}}

#endif
