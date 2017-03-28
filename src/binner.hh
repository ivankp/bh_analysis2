// Written by Ivan Pogrebnyak

#ifndef IVANP_BINNER_HH
#define IVANP_BINNER_HH

#include <tuple>
#include <array>
#include <vector>
#include <iterator>
#include <stdexcept>

#include "axis.hh"
#include "default_bin_filler.hh"
#include "utility.hh"
#include "sequence_traits.hh"

namespace ivanp {

template <typename Axis, bool Uf=true, bool Of=true, bool Ex=false>
struct axis_spec {
  using axis  = Axis;
  using under = std::integral_constant<bool,Uf>;
  using over  = std::integral_constant<bool,Of>;
  using nover = std::integral_constant<axis_size_type,Uf+Of>;
  using excep = std::integral_constant<bool,Ex>;
};

template <typename Bin,
          typename AxesSpecs = std::tuple<axis_spec<uniform_axis<double>>>,
          typename Container = std::vector<Bin>,
          typename Filler = bin_filler<Bin>>
class binner;

template <typename Bin, typename... Ax, typename Container, typename Filler>
class binner<Bin,std::tuple<Ax...>,Container,Filler> {
  static_assert(sizeof...(Ax)>0,"");
public:
  using bin_type = Bin;
  using axes_specs = std::tuple<Ax...>;
  using axes_tuple = std::tuple<typename Ax::axis...>;
  template <unsigned I>
  using axis_spec = std::tuple_element_t<I,axes_specs>;
  template <unsigned I>
  using axis_type = std::decay_t<typename axis_spec<I>::axis>;
  template <unsigned I>
  using edge_type = typename axis_type<I>::edge_type;
  using container_type = Container;
  using filler_type = Filler;
  using value_type = typename container_type::value_type;
  using size_type = ivanp::axis_size_type;
  static constexpr unsigned naxes = sizeof...(Ax);
  using index_array_type = std::array<size_type,naxes>;
  using index_array_cref = const index_array_type&;

  using named_ptr_type = named_ptr<binner>;
  static std::vector<named_ptr_type> all;

private:
  axes_tuple _axes;
  container_type _bins;

  template <size_t I, typename... Args, size_t... A>
  constexpr axis_type<I> _make_axis(
    std::index_sequence<A...>,
    std::tuple<Args...> args
  ) { return axis_type<I>(std::get<A>(args)...); }

  template <size_t I, typename T, size_t N, size_t... A>
  constexpr axis_type<I> _make_axis(
    std::index_sequence<A...>,
    const std::array<T,N>& args
  ) { return axis_type<I>({std::get<A>(args)...}); }

public:
  template <unsigned I=0>
  constexpr const auto& axis() const noexcept { return std::get<I>(_axes); }
  constexpr const auto& axes() const noexcept { return _axes; }

  template <unsigned I=0>
  constexpr size_type nbins() const noexcept {
    static_assert( I < naxes, "Axis index out of range" );
    return axis<I>().nbins() + axis_spec<I>::nover::value;
  }

  template <unsigned I> // including I
  constexpr std::enable_if_t<I!=0,size_type> nbins_left() const noexcept {
    return nbins<I>() * nbins_left<I-1>();
  }
  template <unsigned I> // including I
  constexpr std::enable_if_t<I==0,size_type> nbins_left() const noexcept {
    return nbins<I>();
  }
  template <unsigned I> // excluding I
  constexpr std::enable_if_t<I!=0,size_type> nbins_before() const noexcept {
    return nbins_left<I-1>();
  }
  template <unsigned I> // excluding I
  constexpr std::enable_if_t<I==0,size_type> nbins_before() const noexcept {
    return 1;
  }

  template <unsigned I> // including I
  constexpr std::enable_if_t<I!=naxes-1,size_type> nbins_right() const noexcept {
    return nbins<I>() * nbins_right<I+1>();
  }
  template <unsigned I> // including I
  constexpr std::enable_if_t<I==naxes-1,size_type> nbins_right() const noexcept {
    return nbins<I>();
  }
  template <unsigned I> // excluding I
  constexpr std::enable_if_t<(I<naxes-1),size_type> nbins_after() const noexcept {
    return nbins_right<I+1>();
  }
  template <unsigned I> // excluding I
  constexpr std::enable_if_t<(I>=naxes-1),size_type> nbins_after() const noexcept {
    return 1;
  }

  constexpr size_type nbins_total() const noexcept {
    return nbins_left<naxes-1>();
  }

private:
  template <size_t I, typename A = axis_spec<I>>
  constexpr std::enable_if_t<A::under::value,bool>
  guard_under(size_type bin) const noexcept { return false; }
  template <size_t I, typename A = axis_spec<I>>
  constexpr std::enable_if_t<!A::under::value && !A::excep::value,bool>
  guard_under(size_type bin) const noexcept { return (bin == 0); }
  template <size_t I, typename A = axis_spec<I>>
  constexpr std::enable_if_t<!A::under::value && A::excep::value,bool>
  guard_under(size_type bin) const {
    if (bin == 0) throw std::out_of_range("axis underflow");
    return false;
  }

  template <size_t I, typename A = axis_spec<I>>
  constexpr std::enable_if_t<A::over::value,bool>
  guard_over(size_type bin) const noexcept { return false; }
  template <size_t I, typename A = axis_spec<I>>
  constexpr std::enable_if_t<!A::over::value && !A::excep::value,bool>
  guard_over(size_type bin) const noexcept {
    return (bin > axis<I>().nbins());
  }
  template <size_t I, typename A = axis_spec<I>>
  constexpr std::enable_if_t<!A::over::value && A::excep::value,bool>
  guard_over(size_type bin) const {
    if (bin > axis<I>().nbins()) throw std::out_of_range("axis overflow");
    return false;
  }

  // find_bin_impl --------------------------------------------------
  template <size_t I=0, typename T, typename... Args>
  inline size_type find_bin_impl(const T& x, const Args&... args) const {
    const size_type bin = axis<I>().find_bin(x);
    if (guard_under<I>(bin)) return size_type(-1);
    if (guard_over <I>(bin)) return size_type(-1);
    const size_type bin1 = find_bin_impl<I+1>(args...);
    if (bin1 == size_type(-1)) return size_type(-1);
    return bin - !axis_spec<I>::under::value
      + (axis<I>().nbins() + axis_spec<I>::nover::value) * bin1;
  }
  template <size_t I=0, typename T>
  inline size_type find_bin_impl(const T& x) const {
    const size_type bin = axis<I>().find_bin(x);
    if (guard_under<I>(bin)) return size_type(-1);
    if (guard_over <I>(bin)) return size_type(-1);
    return bin - !axis_spec<I>::under::value;
  }
  // ----------------------------------------------------------------
  template <typename... T, size_t... I>
  constexpr size_type find_bin_tuple(
    const std::tuple<T...>& t, std::index_sequence<I...>)
  const { return find_bin_impl(std::get<I>(t)...); }

  template <typename... T, size_t... I>
  inline size_type fill_bin_tuple(size_type bin,
    const std::tuple<T...>& t, std::index_sequence<I...>
  ) {
    filler_type()(_bins[bin], std::get<I>(t)...);
    return bin;
  }

  template <typename T, typename... TT>
  constexpr size_type index_impl(T i, TT... ii) const noexcept {
    return i + (axis<naxes-sizeof...(TT)-1>().nbins()
             + axis_spec<naxes-sizeof...(TT)-1>::nover::value
             - !axis_spec<naxes-sizeof...(TT)-1>::under::value)
      * index_impl(ii...);
  }
  constexpr size_type index_impl(size_type i) const noexcept { return i; }
  template <size_t... I>
  constexpr size_type index_impl(index_array_cref ia,std::index_sequence<I...>)
  const noexcept { return index_impl(std::get<I>(ia)...); }

public:
  binner() = default;
  ~binner() = default;

  template <typename C=container_type,
            std::enable_if_t<is_std_vector<C>::value>* = nullptr>
  binner(typename Ax::axis... axes)
  : _axes{std::forward<typename Ax::axis>(axes)...}, _bins(nbins_total()) { }
  template <typename C=container_type,
            std::enable_if_t<is_std_array<C>::value>* = nullptr>
  binner(typename Ax::axis... axes)
  : _axes{std::forward<typename Ax::axis>(axes)...}, _bins{} { }

  template <typename Name, typename C=container_type,
            std::enable_if_t<is_std_vector<C>::value>* = nullptr>
  binner(Name&& name, typename Ax::axis... axes)
  : _axes{std::forward<typename Ax::axis>(axes)...}, _bins(nbins_total()) {
    all.emplace_back(this,std::forward<Name>(name));
  }
  template <typename Name, typename C=container_type,
            std::enable_if_t<is_std_array<C>::value>* = nullptr>
  binner(Name&& name, typename Ax::axis... axes)
  : _axes{std::forward<typename Ax::axis>(axes)...}, _bins{} {
    all.emplace_back(this,std::forward<Name>(name));
  }

  template <typename C>
  binner(std::tuple<typename Ax::axis...>&& axes, C&& bins)
  : _axes(axes), _bins(std::forward<C>(bins)) { }

  binner(const binner& o): _axes(o._axes), _bins(o._bins) { }
  binner(binner&& o): _axes(std::move(o._axes)), _bins(std::move(o._bins)) { }
  binner& operator=(const binner& rhs) {
    _axes = rhs._axes;
    _bins = rhs._bins;
    return *this;
  }
  binner& operator=(binner&& rhs) {
    _axes = std::move(rhs._axes);
    _bins = std::move(rhs._bins);
    return *this;
  }
  binner(const std::string& name, const binner& o)
  : _axes(o._axes), _bins(o._bins) {
    all.emplace_back(this,name);
  }

  binner& operator+=(const binner& rhs) {
    if (nbins_total() != rhs.nbins_total()) throw std::length_error(
      "binner::operator+=: nbins_total does not match");
    for (size_type i=nbins_total(); i!=0;) {
      --i;
      _bins[i] += rhs._bins[i];
    }
    return *this;
  }

  inline const container_type& bins() const noexcept { return _bins; }
  inline container_type& bins() noexcept { return _bins; }

  constexpr size_type index(type_to_type<Ax,size_type>... ii) const noexcept {
    return index_impl(ii...);
  }
  constexpr size_type index(index_array_cref ii) const noexcept {
    return index_impl(ii,std::make_index_sequence<naxes>());
  }

  inline const value_type& bin(type_to_type<Ax,size_type>... ii) const {
    return _bins[index_impl(ii...)];
  }
  inline const value_type& bin(index_array_cref ii) const {
    return _bins[index_impl(ii,std::make_index_sequence<naxes>())];
  }

  // find bin -------------------------------------------------------
  template <typename... Args>
  inline size_type find_bin(const Args&... args) const {
    static_assert(sizeof...(Args)==naxes,"");
    return find_bin_impl(args...);
  }
  template <typename... Args>
  inline size_type find_bin(const std::tuple<Args...>& args) const {
    static_assert(sizeof...(Args)==naxes,"");
    return find_bin_tuple(args,std::make_index_sequence<naxes>());
  }

  // fill bin -------------------------------------------------------
  template <typename... Args>
  inline size_type fill_bin(size_type bin, Args&&... args) {
    // NOT safe for out of range indices
    filler_type()(_bins[bin], std::forward<Args>(args)...);
    return bin;
  }
  template <typename... Args>
  inline size_type fill_bin(index_array_cref ia, Args&&... args) {
    // NOT safe for out of range indices
    return fill_bin(_bins[index(ia)], std::forward<Args>(args)...);
  }
  template <typename... Args>
  inline size_type fill_bin_check(size_type bin, Args&&... args) {
    if (bin == size_type(-1)) return bin;
    return fill_bin(bin,std::forward<Args>(args)...);
  }

  template <typename... Args>
  inline std::enable_if_t<sizeof...(Args)==naxes,size_type>
  fill(const Args&... args) {
    const size_type bin = find_bin_impl(args...);
    if (bin == size_type(-1)) return size_type(-1);
    return fill_bin(bin);
  }
  template <typename... Args>
  inline std::enable_if_t<(sizeof...(Args)>naxes),size_type>
  fill(const Args&... args) {
    auto tup = std::forward_as_tuple(args...);
    const size_type bin =
      find_bin_tuple(tup,std::make_index_sequence<naxes>());
    if (bin == size_type(-1)) return size_type(-1);
    return fill_bin_tuple( bin,
      tup, ivanp::seq::make_index_range<naxes,sizeof...(Args)>()
    );
  }
  template <typename... Args>
  inline std::enable_if_t<(sizeof...(Args)>=naxes),size_type>
  operator()(const Args&... args) {
    return fill(args...);
  }

  // Algorithms -----------------------------------------------------
  template <unsigned I=0> void integrate_right() {
    const size_type nb = nbins_before<I>();
    const size_type n  = nbins<I>()-1;
    const size_type na = nbins_after<I>();

    auto begin = _bins.begin();
    for (size_type a=na; a; --a) {
      for (size_type b=nb; b;) { --b;
        auto prev = begin;
        for (size_type i=n; i; --i) {
          auto it = std::next(prev,nb);
          *it += *prev;
          prev = it;
        }
        if (b) ++begin;
        else begin = std::next(prev);
      }
    }
  }
  template <unsigned I=0> void integrate_left() {
    const size_type nb = nbins_before<I>();
    const size_type n  = nbins<I>()-1;
    const size_type na = nbins_after<I>();

    auto begin = _bins.rbegin();
    for (size_type a=na; a; --a) {
      for (size_type b=nb; b;) { --b;
        auto prev = begin;
        for (size_type i=n; i; --i) {
          auto it = std::next(prev,nb);
          *it += *prev;
          prev = it;
        }
        if (b) ++begin;
        else begin = std::next(prev);
      }
    }
  }
};

template <typename B, typename... A, typename C, typename F>
std::vector<typename binner<B,std::tuple<A...>,C,F>::named_ptr_type>
binner<B,std::tuple<A...>,C,F>::all;

// ##################################################################

// Metafunctions

template <typename Binner, typename... AxesSpecs> struct add_axes;
template <typename Bin, typename... Ax, typename Container, typename Filler,
          typename... AxesSpecs>
struct add_axes<binner<Bin,std::tuple<Ax...>,Container,Filler>,AxesSpecs...> {
  using type = binner<Bin,std::tuple<Ax...,AxesSpecs...>,Container,Filler>;
};
template <typename Binner, typename... AxesSpecs>
using add_axes_t = typename add_axes<Binner,AxesSpecs...>::type;

} // end namespace ivanp

#endif
