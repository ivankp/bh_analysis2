#ifndef IVANP_TYPE_TRAITS_HH
#define IVANP_TYPE_TRAITS_HH

#include <type_traits>

namespace ivanp {

// ******************************************************************

template <typename... T> struct make_void { typedef void type; };
template <typename... T> using void_t = typename make_void<T...>::type;

template <typename Old, typename New> using type_to_type = New;
template <size_t I, typename T> using index_to_type = T;

template <typename T>
struct add_const_to_ref { using type = T; };
template <typename T>
struct add_const_to_ref<T&> { using type = const T&; };
template <typename T>
using add_const_to_ref_t = typename add_const_to_ref<T>::type;

template <typename T> struct remove_rref      { using type = T; };
template <typename T> struct remove_rref<T&&> { using type = T; };
template <typename T> using remove_rref_t = typename remove_rref<T>::type;

// trait composition ************************************************

template<class...> struct conjunction: std::true_type { };
template<class B1> struct conjunction<B1>
: std::integral_constant<bool,bool(B1::value)> { };
template<class B1, class... Bn> struct conjunction<B1, Bn...> 
: std::conditional_t<bool(B1::value),
  conjunction<Bn...>, conjunction<B1> > {};

template<class...> struct disjunction: std::false_type { };
template<class B1> struct disjunction<B1>
: std::integral_constant<bool,bool(B1::value)> { };
template<class B1, class... Bn> struct disjunction<B1, Bn...> 
: std::conditional_t<bool(B1::value),
  disjunction<B1>, disjunction<Bn...>>  { };

// IS ***************************************************************

#ifdef _GLIBCXX_ARRAY
template <typename> struct is_std_array: std::false_type { };
template <typename T, size_t N>
struct is_std_array<std::array<T,N>>: std::true_type { };
#endif

#ifdef _GLIBCXX_VECTOR
template <typename> struct is_std_vector: std::false_type { };
template <typename T, typename Alloc>
struct is_std_vector<std::vector<T,Alloc>>: std::true_type { };
#endif

// Tuple ************************************************************

#ifdef _GLIBCXX_TUPLE
template <size_t... I, typename... T>
auto make_ref_subtuple(std::tuple<T...>& t, std::index_sequence<I...>) {
  return std::tie( std::get<I>(t)... );
}

template <typename Tup, typename Elems> struct subtuple;
template <typename... T, size_t... I>
struct subtuple<std::tuple<T...>,std::index_sequence<I...>> {
  using type = std::tuple<std::tuple_element_t<I,std::tuple<T...>>...>;
};
template <typename Tup, typename Elems>
using subtuple_t = typename subtuple<Tup,Elems>::type;

namespace detail {
template <typename T, typename S> struct tuple_of_same_impl { };
template <typename T, size_t... I>
struct tuple_of_same_impl<T,std::index_sequence<I...>> {
  template <typename U, size_t J> using identity = U;
  using type = std::tuple<identity<T,I>...>;
};
}

template <typename T, size_t N>
struct tuple_of_same
: detail::tuple_of_same_impl<T,std::make_index_sequence<N>> { };

template <typename T, size_t N>
using tuple_of_same_t = typename tuple_of_same<T,N>::type;
#endif

// Expression traits ************************************************

// void_t technique from Walter Brown
// https://www.youtube.com/watch?v=Am2is2QCvxY
// https://www.youtube.com/watch?v=a0FliKwcwXE

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

template <typename, typename, typename = void> // x += r
struct has_plus_eq : std::false_type { };
template <typename T, typename R>
struct has_plus_eq<T,R,
  void_t<decltype( std::declval<T&>()+=std::declval<R>() )>
> : std::true_type { };

template <typename, typename, typename = void> // x -= r
struct has_minus_eq : std::false_type { };
template <typename T, typename R>
struct has_minus_eq<T,R,
  void_t<decltype( std::declval<T&>()-=std::declval<R>() )>
> : std::true_type { };

template <typename T, typename... Args> // x(args...)
class is_callable {
  template <typename, typename = void>
  struct impl: std::false_type { };
  template <typename U>
  struct impl<U,
    void_t<decltype( std::declval<U&>()(std::declval<Args>()...) )>
  > : std::true_type { };
public:
  static constexpr bool value = impl<T>::value;
};

template <typename T, typename... Args> // x(args...)
class is_constructible {
  template <typename, typename = void>
  struct impl: std::false_type { };
  template <typename U>
  struct impl<U,
    void_t<decltype( U(std::declval<Args>()...) )>
  > : std::true_type { };
public:
  static constexpr bool value = impl<T>::value;
};

// ******************************************************************

} // end namespace ivanp

#endif
