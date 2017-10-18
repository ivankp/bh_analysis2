#ifndef IVANP_MAYBE_VALID_HH
#define IVANP_MAYBE_VALID_HH

#include "meta.hh"

namespace ivanp {

// based on hana::sfinae
template <typename F>
class maybe_valid_wrap {
  F f;

  template <typename... Args>
  using ret_t = decltype(f(std::declval<Args&&>()...));

public:
  template <typename... Args>
  using is_valid = is_detected<ret_t,Args...>;

  template <typename T>
  constexpr maybe_valid_wrap(T&& x): f(std::forward<T>(x)) { }

  template <typename... Args>
  constexpr std::enable_if_t<
    !std::is_void<ret_t<Args...>>::value, just_value<ret_t<Args...>> >
  operator()(Args&&... args) noexcept(noexcept(f(std::forward<Args>(args)...)))
  { return f(std::forward<Args>(args)...); }

  template <typename... Args>
  constexpr std::enable_if_t<
    std::is_void<ret_t<Args...>>::value, just_value<void> >
  operator()(Args&&... args) noexcept(noexcept(f(std::forward<Args>(args)...)))
  { f(std::forward<Args>(args)...); return { }; }

  template <typename... Args>
  constexpr std::enable_if_t<!is_valid<Args...>::value,nothing>
  operator()(Args&&...) noexcept { return {}; }
};

template <typename F>
constexpr maybe_valid_wrap<F> maybe_valid(F&& f)
noexcept(noexcept(maybe_valid_wrap<F>{std::forward<F>(f)}))
{ return { std::forward<F>(f) }; }

template <typename... Args>
struct is_valid_pred {
  template <typename F>
  using type = typename F::template is_valid<Args...>;
};

template <typename... Fs>
class first_valid_wrap {
  std::tuple<maybe_valid_wrap<Fs>...> fs;

  template <typename... Args>
  using index = first_index<
    is_valid_pred<Args...>::template type, maybe_valid_wrap<Fs>...  >;

public:
  template <typename... _Fs>
  first_valid_wrap(_Fs&&... fs): fs(std::forward<_Fs>(fs)...) { }

  template <typename... Args, typename Index = index<Args...>,
            enable_if_just_t<Index,size_t> I = Index::type::value>
  constexpr decltype(auto) operator()(Args&&... args)
  noexcept(noexcept(std::get<I>(fs)(std::forward<Args>(args)...)))
  { return std::get<I>(fs)(std::forward<Args>(args)...); }

  template <typename... Args, typename Index = index<Args...>>
  constexpr enable_if_nothing_t<Index,nothing>
  operator()(Args&&... args) noexcept { return { }; }
};

template <typename... Fs>
constexpr first_valid_wrap<Fs...> first_valid(Fs&&... fs)
noexcept(noexcept(first_valid_wrap<Fs...>{std::forward<Fs>(fs)...}))
{ return { std::forward<Fs>(fs)... }; }

#ifdef SFINAE_EXPR
#error macro named 'SFINAE_EXPR' already defined
#endif
#define SFINAE_EXPR(EXPR, ...) \
    [](__VA_ARGS__) noexcept(noexcept(EXPR))->void_t<decltype(EXPR)> { EXPR; }

}

#endif
