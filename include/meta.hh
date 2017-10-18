#ifndef IVANP_META_HH
#define IVANP_META_HH

#include <type_traits>

namespace ivanp {

// void_t
template <typename... T> struct make_void { typedef void type; };
template <typename... T> using void_t = typename make_void<T...>::type;

// curry
template <template<typename,typename> typename Pred, typename T1>
struct bind_first {
  template <typename T2> using type = Pred<T1,T2>;
};
template <template<typename,typename...> typename Pred, typename... T2>
struct bind_tail {
  template <typename T1> using type = Pred<T1,T2...>;
};

// Logical ==========================================================
// bool const
template <bool B> using bool_constant = std::integral_constant<bool, B>;

// AND
template <class...> struct conjunction: std::true_type { };
template <class B1> struct conjunction<B1>: B1 { };
template <class B1, class... Bs>
struct conjunction<B1, Bs...> 
  : std::conditional_t<bool(B1::value), conjunction<Bs...>, B1> { };
// OR
template <class...> struct disjunction: std::false_type { };
template <class B1> struct disjunction<B1>: B1 { };
template <class B1, class... Bs>
struct disjunction<B1, Bs...> 
  : std::conditional_t<bool(B1::value), B1, disjunction<Bs...>> { };
// NOT
template <class B> struct negation: bool_constant<!bool(B::value)> { };

// Only last
template <typename...> struct only_last: std::false_type { };
template <typename B1> struct only_last<B1>: bool_constant<B1::value> { };
template <typename B1, typename... Bs>
struct only_last<B1, Bs...>
  : bool_constant<!B1::value && only_last<Bs...>::value> { };

// Maybe ============================================================
struct nothing { };
template <typename T> struct just { using type = T; };

template <typename T> struct is_nothing: std::false_type { };
template <> struct is_nothing<nothing>: std::true_type { };

template <typename T> struct is_just: std::false_type { };
template <typename T> struct is_just<just<T>>: std::true_type { };

template <typename T> struct maybe {
  static_assert(is_just<T>::value || is_nothing<T>::value);
  using type = T;
};

// Maybe value ------------------------------------------------------
template <typename T> struct just_value {
  constexpr const T& operator*() const noexcept { return value; }
  constexpr T& operator*() noexcept { return value; }
private:
  T value;
};
template <> struct just_value<void> {
  constexpr void operator*() const noexcept { }
  constexpr void operator*() noexcept { }
};
template <typename T> struct is_just_value: std::false_type { };
template <typename T> struct is_just_value<just_value<T>>: std::true_type { };

template <typename T> struct maybe_value : public T {
  static_assert(is_just_value<T>::value || is_nothing<T>::value);
};

// maybe enable -----------------------------------------------------
template <typename M, typename T = void>
using enable_if_just_t = std::enable_if_t<is_just<M>::value,T>;
template <typename M, typename T = void>
using enable_if_nothing_t = std::enable_if_t<is_nothing<M>::value,T>;

// maybe is ---------------------------------------------------------
template <template<typename> typename Pred, typename M, bool N = false>
struct maybe_is;
template <template<typename> typename Pred, bool N>
struct maybe_is<Pred,nothing,N>: bool_constant<N> { };
template <template<typename> typename Pred, typename T, bool N>
struct maybe_is<Pred,just<T>,N>: bool_constant<Pred<T>::value> { };

#ifdef __cpp_variable_templates
template <template<typename> typename Pred, typename M, bool N = false>
constexpr bool maybe_is_v = maybe_is<Pred,M,N>::value;
#endif

// First in pack ====================================================
template <typename...> struct first: maybe<nothing> { };
template <typename T, typename... Other>
struct first<T,Other...>: maybe<just<T>> { };
template <typename... Args> using first_t = typename first<Args...>::type;

// Find first =======================================================
// first pack element matching predicate
template <template<typename> typename Pred, typename... Args>
struct find_first: maybe<nothing> { };

template <template<typename> typename Pred, typename Arg1, typename... Args>
class find_first<Pred,Arg1,Args...> {
  template <typename, typename = void>
  struct impl: find_first<Pred,Args...> { };
  template <typename Arg>
  struct impl<Arg,std::enable_if_t<Pred<Arg>::value>>: maybe<just<Arg>> { };
public:
  using type = typename impl<Arg1>::type;
};

template <template<typename> typename Pred, typename... Args>
using find_first_t = typename find_first<Pred,Args...>::type;

// Find first index =================================================
namespace detail {
template <template<typename> typename Pred, size_t I, typename... Args>
struct first_index_impl: maybe<nothing> { };

template <template<typename> typename Pred, size_t I,
          typename Arg1, typename... Args>
class first_index_impl<Pred,I,Arg1,Args...> {
  template <typename, typename = void>
  struct impl: first_index_impl<Pred,I+1,Args...> { };
  template <typename Arg>
  struct impl<Arg,std::enable_if_t<Pred<Arg>::value>>
  : maybe<just<std::integral_constant<size_t,I>>> { };
public:
  using type = typename impl<Arg1>::type;
};
}

template <template<typename> typename Pred, typename... Args>
using first_index = typename detail::first_index_impl<Pred,0,Args...>::type;

// Detect ===========================================================
// http://en.cppreference.com/w/cpp/experimental/is_detected
namespace detail {
template <class Default, class AlwaysVoid,
          template<class...> class Op, class... Args>
struct detector {
  using value_t = std::false_type;
  using type = Default;
};

template <class Default, template<class...> class Op, class... Args>
struct detector<Default, void_t<Op<Args...>>, Op, Args...> {
  using value_t = std::true_type;
  using type = Op<Args...>;
};
} // end namespace detail

template <template<class...> class Op, class... Args>
using is_detected = typename detail::detector<nothing, void, Op, Args...>::value_t;

template <template<class...> class Op, class... Args>
using detected_t = typename detail::detector<nothing, void, Op, Args...>::type;

template <class Default, template<class...> class Op, class... Args>
using detected_or = detail::detector<Default, void, Op, Args...>;

template <class Default, template<class...> class Op, class... Args>
using detected_or_t = typename detected_or<Default, Op, Args...>::type;

template <class Expected, template<class...> class Op, class... Args>
using is_detected_exact = std::is_same<Expected, detected_t<Op, Args...>>;

template <class To, template<class...> class Op, class... Args>
using is_detected_convertible = std::is_convertible<detected_t<Op, Args...>, To>;

#ifdef __cpp_variable_templates
template <template<class...> class Op, class... Args>
constexpr bool is_detected_v = is_detected<Op, Args...>::value;

template <class Expected, template<class...> class Op, class... Args>
constexpr bool is_detected_exact_v =
  is_detected_exact<Expected, Op, Args...>::value;

template <class To, template<class...> class Op, class... Args>
constexpr bool is_detected_convertible_v =
  is_detected_convertible<To, Op, Args...>::value;
#endif

// Enable case ======================================================
namespace detail {
template <template<size_t,class...> class, typename, typename...>
struct is_case_impl;
template <template<size_t,class...> class Pred, size_t... I, typename... Args>
struct is_case_impl<Pred,std::index_sequence<I...>,Args...>
  : only_last<Pred<I,Args...>...> { };
} // end namespace detail

template <template<size_t,class...> class Pred, size_t I, typename... Args>
using is_case = detail::is_case_impl<Pred, std::make_index_sequence<I+1>, Args...>;

template <template<size_t,class...> class Pred, size_t I, typename... Args>
using enable_case = std::enable_if_t<is_case<Pred,I,Args...>::value>;

// template <template<typename...> typename... Preds>
// struct pred_conjunction {
//   template <typename... Args> struct pred: conjunction<Preds<Args...>...> { };
// };
//
// template <template<typename...> typename Pred>
// struct pred_negation {
//   template <typename... Args> struct pred: negation<Pred<Args...>> { };
// };

} // end namespace ivanp

#endif
