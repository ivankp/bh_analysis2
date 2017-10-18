#ifndef IVANP_OPT_INIT_HH
#define IVANP_OPT_INIT_HH

#include "program_options/opt_parser.hh"

namespace ivanp { namespace po {
namespace _ {

template <typename... Args> class opt_init_base {
  std::tuple<Args...> args;
  template <typename T, size_t... I>
  inline void construct(T& x, std::index_sequence<I...>) {
    x = { std::get<I>(args)... };
  }
  template <typename T>
  using direct = bool_constant<sizeof...(Args)==1 && maybe_is_v<
      bind_first<is_assignable,T>::template type, first_t<Args...>
    >>;
public:
  template <typename... TT>
  opt_init_base(TT&&... xx): args(std::forward<TT>(xx)...) { }
  template <typename... TT>
  opt_init_base(const std::tuple<TT...>& tup): args(tup) { }
  template <typename... TT>
  opt_init_base(std::tuple<TT...>&& tup): args(std::move(tup)) { }
  template <typename T>
  inline std::enable_if_t<direct<T>::value> construct(T& x) {
    x = std::get<0>(args);
  }
  template <typename T>
  inline std::enable_if_t<!direct<T>::value> construct(T& x) {
    construct(x,std::index_sequence_for<Args...>{});
  }
};
template <> struct opt_init_base<> {
  template <typename T> inline void construct(T& x) const { x = { }; }
};

template <typename... Args> struct switch_init: opt_init_base<Args...> {
  using opt_init_base<Args...>::opt_init_base;
};
template <typename... Args> struct default_init: opt_init_base<Args...> {
  using opt_init_base<Args...>::opt_init_base;
};

// detect -----------------------------------------------------------
template <typename T> struct is_switch_init : std::false_type { };
template <typename... T>
struct is_switch_init<switch_init<T...>> : std::true_type { };

template <typename T> struct is_default_init : std::false_type { };
template <typename... T>
struct is_default_init<default_init<T...>> : std::true_type { };

} // end _::

// Factories --------------------------------------------------------
template <typename... Args>
inline _::switch_init<std::decay_t<Args>...> switch_init(Args&&... args) {
  return { std::forward<Args>(args)... };
}
template <typename... Args>
inline _::switch_init<Args...> switch_init(const std::tuple<Args...>& args) {
  return { args };
}
template <typename... Args>
inline _::switch_init<Args...> switch_init(std::tuple<Args...>&& args) {
  return { std::move(args) };
}

template <typename... Args>
inline _::default_init<std::decay_t<Args>...> default_init(Args&&... args) {
  return { std::forward<Args>(args)... };
}
template <typename... Args>
inline _::default_init<Args...> default_init(const std::tuple<Args...>& args) {
  return { args };
}
template <typename... Args>
inline _::default_init<Args...> default_init(std::tuple<Args...>&& args) {
  return { std::move(args) };
}

}} // po::

#endif
