#ifndef IVANP_TYPE_TRAITS_HH
#define IVANP_TYPE_TRAITS_HH

#include "meta.hh"

namespace ivanp {

template <typename T> struct rm_rref { using type = T; };
template <typename T> struct rm_rref<T&&> { using type = T; };
template <typename T> using rm_rref_t = typename rm_rref<T>::type;

// Expression traits ================================================
template <typename T, typename... Args> // T(Args...)
class is_callable {
  template <typename, typename = void> struct impl: std::false_type { };
  template <typename U> struct impl<U,
    void_t<decltype( std::declval<U&>()(std::declval<Args>()...) )>
  > : std::true_type { };
public:
  static constexpr bool value = impl<T>::value;
};

template <typename T, typename... Args> // T(Args...)
class is_constructible {
  template <typename, typename = void> struct impl: std::false_type { };
  template <typename U> struct impl<U,
    void_t<decltype( U(std::declval<Args>()...) )>
  > : std::true_type { };
public:
  static constexpr bool value = impl<T>::value;
};

template <typename T, typename Arg> // T = Arg
class is_assignable {
  template <typename, typename = void> struct impl: std::false_type { };
  template <typename U> struct impl<U,
    void_t<decltype( std::declval<U&>() = std::declval<Arg>() )>
  > : std::true_type { };
public:
  static constexpr bool value = impl<T>::value;
};

// value_type =======================================================
namespace detail {
template <typename, typename = void>
struct maybe_value_type { using type = nothing; };
template <typename T>
struct maybe_value_type<T, void_t<typename T::value_type> > {
  using type = just<typename T::value_type>;
};
}
template <typename T>
using m_value_type = typename detail::maybe_value_type<T>::type;

// Remove const from tuple elements =================================
template <typename T> struct rm_const_elements {
  using type = T;
};
template <typename T1, typename T2>
struct rm_const_elements<std::pair<T1,T2>> {
  using type = std::pair<std::remove_const_t<T1>,std::remove_const_t<T2>>;
};
template <typename... Ts>
struct rm_const_elements<std::tuple<Ts...>> {
  using type = std::tuple<std::remove_const_t<Ts>...>;
};
template <typename T, size_t N>
struct rm_const_elements<std::array<T,N>> {
  using type = std::array<std::remove_const_t<T>,N>;
};
template <typename T>
using rm_const_elements_t = typename rm_const_elements<T>::type;

// Detect STD types =================================================
template <typename T> struct is_tuple: std::false_type { };
template <typename... T> struct is_tuple<std::tuple<T...>>: std::true_type { };

template <typename T> struct is_std_array: std::false_type { };
template <typename T, size_t N>
struct is_std_array<std::array<T,N>>: std::true_type { };

template <typename T> struct is_std_vector: std::false_type { };
template <typename T, typename Alloc>
struct is_std_vector<std::vector<T,Alloc>>: std::true_type { };

} // end namespace ivanp

#include "old_type_traits.hh"
#include "op_traits.hh"

#endif

