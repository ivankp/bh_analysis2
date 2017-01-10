// Written by Ivan Pogrebnyak

#ifndef IVANP_EXPRESSION_TRAITS_HH
#define IVANP_EXPRESSION_TRAITS_HH

#include <type_traits>

// DEFINE_TRAIT from:
// https://www.youtube.com/watch?v=CZi6QqZSbFg
// https://github.com/CppCon/CppCon2015/blob/master/Presentations/
// C++ Metaprogramming/C++ Metaprogramming - Fedor Pikus - CppCon 2015.pdf

#define DEFINE_UNARY_TRAIT(NAME, EXPR) \
template <typename T> struct NAME { \
  typedef char yes; \
  typedef char no[2]; \
  template <typename U> static auto f(U&& x) -> decltype(EXPR, NAME::yes()); \
  template <typename U> static no&  f(...); \
  enum { value = sizeof(NAME::f<T>(std::declval<T>())) \
              == sizeof(NAME::yes) }; \
};

#define DEFINE_BINARY_TRAIT(NAME, EXPR) \
template <typename T1, typename T2> struct NAME { \
  typedef char yes; \
  typedef char no[2]; \
  template <typename U1, typename U2> \
  static auto f(U1&& x1, U2&& x2) -> decltype(EXPR, NAME::yes()); \
  template <typename U1, typename U2> static no& f(...); \
  enum { value = sizeof(NAME::f<T1,T2>(std::declval<T1>(),std::declval<T2>())) \
              == sizeof(NAME::yes) }; \
};

#define DEFINE_VARIADIC_TRAIT(NAME, EXPR) \
template <typename T, typename... TT> struct NAME { \
  typedef char yes; \
  typedef char no[2]; \
  template <typename U, typename... UU> \
  static auto f(U&& x, UU&&... xx) -> decltype(EXPR, NAME::yes()); \
  template <typename U, typename... UU> static no& f(...); \
  enum { value = sizeof(NAME::f<T,TT...>(std::declval<T>(),std::declval<TT>()...)) \
              == sizeof(NAME::yes) }; \
};

namespace ivanp {

DEFINE_UNARY_TRAIT(has_op_pre_increment,  ++x)
DEFINE_UNARY_TRAIT(has_op_post_increment, x++)
DEFINE_UNARY_TRAIT(has_op_pre_decrement,  --x)
DEFINE_UNARY_TRAIT(has_op_post_decrement, x--)

DEFINE_BINARY_TRAIT(has_op_plus_eq, x1+=x2)
DEFINE_BINARY_TRAIT(has_op_minus_eq, x1-=x2)

DEFINE_VARIADIC_TRAIT(is_callable, x(xx...))
DEFINE_VARIADIC_TRAIT(is_constructible, T(xx...))

} // end namespace

#endif
