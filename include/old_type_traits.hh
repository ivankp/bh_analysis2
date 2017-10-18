#ifndef IVANP_OLD_TYPE_TRAITS_HH
#define IVANP_OLD_TYPE_TRAITS_HH

#include <type_traits>
#include <tuple>

namespace ivanp {

// ******************************************************************

template <typename Old, typename New> using type_to_type = New;
template <size_t I, typename T> using index_to_type = T;

template <typename T>
struct add_const_to_ref { using type = T; };
template <typename T>
struct add_const_to_ref<T&> { using type = const T&; };
template <typename T>
using add_const_to_ref_t = typename add_const_to_ref<T>::type;

// make_subtuple
template <typename Tuple, size_t... I>
inline auto make_subtuple(Tuple&& tup, std::index_sequence<I...>) {
  return std::tuple<
    rm_rref_t<std::tuple_element_t<I,Tuple>>...
  >{ std::get<I>(tup)... };
}

// forward_subtuple
template <typename Tuple, size_t... I>
inline auto forward_subtuple(Tuple&& tup, std::index_sequence<I...>) {
  return std::forward_as_tuple( std::get<I>(tup)... );
}

// tie_some
template <typename... T, size_t... I>
inline auto tie_some(T&... args, std::index_sequence<I...>) {
  auto&& tmp = std::tie(args...);
  return std::tie( std::get<I>(tmp)... );
}

// forward_some_as_tuple
template <typename... T, size_t... I>
inline auto forward_some_as_tuple(T&&... args, std::index_sequence<I...>) {
  auto&& tmp = std::forward_as_tuple(std::forward<T>(args)...);
  return std::forward_as_tuple( std::get<I>(tmp)... );
}

// subtuple
template <typename Tup, typename Elems> struct subtuple;
template <typename... T, size_t... I>
struct subtuple<std::tuple<T...>,std::index_sequence<I...>> {
  using type = std::tuple<std::tuple_element_t<I,std::tuple<T...>>...>;
};
template <typename Tup, typename Elems>
using subtuple_t = typename subtuple<Tup,Elems>::type;

// pack_is_tuple
template <typename...> struct pack_is_tuple : std::false_type { };
template <typename T> struct pack_is_tuple<T>
: std::integral_constant<bool, is_tuple<T>::value> { };

// tuple_of_same
template <typename T, size_t N>
class tuple_of_same {
  template <typename Seq> struct impl { };
  template <size_t... I> struct impl<std::index_sequence<I...>> {
    using type = std::tuple<index_to_type<I,T>...>;
  };
public:
  using type = typename impl<std::make_index_sequence<N>>::type;
};
template <typename T, size_t N>
using tuple_of_same_t = typename tuple_of_same<T,N>::type;

} // end namespace ivanp

#endif
