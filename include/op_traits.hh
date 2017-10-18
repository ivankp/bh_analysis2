#ifndef IVANP_OP_TRAITS
#define IVANP_OP_TRAITS

#include "meta.hh"

namespace ivanp {

template <typename, typename = void> // ++x
struct has_pre_increment : std::false_type { };
template <typename T>
struct has_pre_increment<T,
  void_t<decltype( ++std::declval<T&>() )>
> : std::true_type { };

template <typename, typename = void> // x++
struct has_post_increment : std::false_type { };
template <typename T>
struct has_post_increment<T,
  void_t<decltype( std::declval<T&>()++ )>
> : std::true_type { };

template <typename, typename = void> // --x
struct has_pre_decrement : std::false_type { };
template <typename T>
struct has_pre_decrement<T,
  void_t<decltype( --std::declval<T&>() )>
> : std::true_type { };

template <typename, typename = void> // x--
struct has_post_decrement : std::false_type { };
template <typename T>
struct has_post_decrement<T,
  void_t<decltype( std::declval<T&>()-- )>
> : std::true_type { };

template <typename, typename, typename = void> // x += rhs
struct has_plus_eq : std::false_type { };
template <typename T, typename R>
struct has_plus_eq<T,R,
  void_t<decltype( std::declval<T&>()+=std::declval<R>() )>
> : std::true_type { };

template <typename, typename, typename = void> // x -= rhs
struct has_minus_eq : std::false_type { };
template <typename T, typename R>
struct has_minus_eq<T,R,
  void_t<decltype( std::declval<T&>()-=std::declval<R>() )>
> : std::true_type { };

}

#endif
