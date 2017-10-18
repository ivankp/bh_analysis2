#ifndef IVANP_OPT_DEF_HH
#define IVANP_OPT_DEF_HH

#include "program_options/fwd/opt_def.hh"

namespace ivanp { namespace po {

struct no_var_t { } no_var;
template <typename T> struct is_no_var : std::false_type { };
template <> struct is_no_var<no_var_t> : std::true_type { };

// Opt def props ---------------------------------------------------

namespace _ {

struct name { std::string name; };
struct named { named(const name&) noexcept { } };
template <typename T> struct is_named : std::false_type { };
template <> struct is_named<named> : std::true_type { };
template <> struct is_named<name> : std::true_type { };

struct multi { };
template <typename T> struct is_multi : std::false_type { };
template <> struct is_multi<multi> : std::true_type { };

struct pos { };
template <typename T> struct is_pos : std::false_type { };
template <> struct is_pos<pos> : std::true_type { };

struct npos { unsigned npos; };
template <typename T> struct is_npos : std::false_type { };
template <> struct is_npos<npos> : std::true_type { };

struct req { };
template <typename T> struct is_req : std::false_type { };
template <> struct is_req<req> : std::true_type { };

template <typename T, typename F>
struct parser {
  F f;
  template <typename G> constexpr parser(G&& g): f(std::forward<G>(g)) { }
  inline auto operator()(const char* str, T& x)
  noexcept(noexcept(f(str,x))) { return f(str,x); }
};
template <typename T, typename F>
struct parser<T,F&> {
  F& f;
  constexpr parser(F& g): f(g) { }
  inline auto operator()(const char* str, T& x)
  noexcept(noexcept(f(str,x))) { return f(str,x); }
};

template <typename T> struct is_parser {
  template <typename F>
  using type = bool_constant<is_callable<F,const char*,T&>::value>;
};
template <typename T> struct is_parser<const T> {
  template <typename F>
  using type = std::false_type;
};

} // end namespace _

template <typename... Args>
inline _::name name(Args&&... args) { return {{std::forward<Args>(args)...}}; }
constexpr _::multi multi() noexcept { return {}; }
constexpr _::pos   pos() noexcept { return {}; }
constexpr _::npos  pos(unsigned n) noexcept { return { n }; }
constexpr _::req   req() noexcept { return {}; }

namespace detail {

// Option definition objects ----------------------------------------
// These provide common interface for invoking argument parsers
// and assigning new values to recepients via pointers
// These are created as a result of calling program_options::operator()

template <typename T, typename... Props>
class opt_def_impl final: public opt_def, Props... {
  T *x; // recepient of parsed value

public:
  using type = std::decay_t<T>;

#define OPT_PROP_TYPE(NAME) \
  using NAME##_t = find_first_t<_::is_##NAME,Props...>;

  OPT_PROP_TYPE(switch_init)
  OPT_PROP_TYPE(default_init)
  OPT_PROP_TYPE(pos)
  OPT_PROP_TYPE(npos)
  OPT_PROP_TYPE(req)
  OPT_PROP_TYPE(multi)
  OPT_PROP_TYPE(named)

#undef OPT_PROP_TYPE

  using parser_t = find_first_t<_::is_parser<T>::template type,Props...>;

  static constexpr bool _is_pos =
    is_just<pos_t>::value || is_just<npos_t>::value;
  static constexpr bool _is_switch = std::is_same<type,bool>::value;

  static_assert( !(_is_pos && _is_switch),
    ASSERT_MSG("definition of positional switch option"));

private:
  // parse ----------------------------------------------------------
  template <typename P = parser_t> inline std::enable_if_t<
    is_just<P>::value && !_is_switch>
  parse_impl(const char* arg) { parser_t::type::operator()(arg,*x); }
  template <typename P = parser_t> inline std::enable_if_t<
    is_nothing<P>::value && !_is_switch>
  parse_impl(const char* arg) {
    static_assert(!is_no_var<T>::value,
      ASSERT_MSG("po::no_var variable requires a user-define parser"));
    ivanp::po::arg_parser(arg,*x);
  }
  template <bool S = _is_switch> static inline std::enable_if_t<S>
  parse_impl(const char* arg) noexcept { }

  // switch ---------------------------------------------------------
  template <typename U = switch_init_t> inline enable_if_just_t<U>
  as_switch_impl() { U::type::construct(*x); }
  template <typename U = switch_init_t> inline std::enable_if_t<
    is_nothing<U>::value && std::is_same<type,bool>::value>
  as_switch_impl() const { (*x) = true; }
  template <typename U = switch_init_t> [[noreturn]]
  inline std::enable_if_t<
    is_nothing<U>::value && !std::is_same<type,bool>::value>
  as_switch_impl() const { throw error(name + " without value"); }

  // default --------------------------------------------------------
  template <typename U = default_init_t> inline enable_if_just_t<U>
  default_init_impl() { U::type::construct(*x); }
  template <typename U = default_init_t> static inline enable_if_nothing_t<U>
  default_init_impl() noexcept { }

public:
  // ----------------------------------------------------------------
  template <typename... M>
  opt_def_impl(T* x, std::string&& descr, M&&... m)
  : opt_def(std::move(descr)), Props(std::forward<M>(m))..., x(x) { };

  inline void parse(const char* arg) { parse_impl(arg); ++count; }
  inline void as_switch() { as_switch_impl(); ++count; }
  inline void default_init() { default_init_impl(); }

  inline bool is_switch() const noexcept { return _is_switch; }
  inline bool is_switch_init() const noexcept {
    return is_just<switch_init_t>::value;
  }

  inline bool is_multi() const noexcept {
    return disjunction< // allow `std::vector`s to automatically be multi
        is_just<multi_t>, is_std_vector<T>
      >::value;
  }
  inline bool is_pos() const noexcept { return _is_pos; }
  inline bool is_pos_end() const noexcept { return is_just<pos_t>::value; }
  inline bool is_req() const noexcept { return is_just<req_t>::value; }
  inline bool is_signed() const noexcept {
    return std::is_signed<type>::value;
  }
  inline bool is_named() const noexcept { return is_just<named_t>::value; }
};

// Factory ----------------------------------------------------------
template <typename T, typename Prop, typename = void>
struct prop_type_subst { using type = std::decay_t<Prop>; };
template <typename T, typename Prop>
struct prop_type_subst<T, Prop,
  std::enable_if_t< _::is_parser<T>::template type<Prop>::value >
> { using type = _::parser<T,rm_rref_t<Prop>>; };
template <typename T, typename Prop>
struct prop_type_subst<T, Prop,
  std::enable_if_t< _::is_named<std::decay_t<Prop>>::value >
> { using type = _::named; };

template <typename T, typename Tuple, size_t... I>
inline auto make_opt_def(
  T* x, std::string&& descr, Tuple&& tup, std::index_sequence<I...>
) {
  return new opt_def_impl<T,
    typename prop_type_subst<T,
      std::tuple_element_t<I,std::decay_t<Tuple>>
    >::type... >
  ( x, std::move(descr), std::get<I>(std::forward<Tuple>(tup))... );
}

} // end namespace detail

}}

#endif
