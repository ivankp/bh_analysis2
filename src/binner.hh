// Written by Ivan Pogrebnyak

#ifndef IVANP_BINNER_HH
#define IVANP_BINNER_HH

#include <tuple>
#include <array>

#include "axis.hh"
#include "default_bin_filler.hh"
#include "utility.hh"

namespace ivanp {

template <typename Axis, bool Uf=true, bool Of=true>
struct axis_spec {
  using axis  = Axis;
  using under = std::integral_constant<bool,Uf>;
  using over  = std::integral_constant<bool,Of>;
  using nover = std::integral_constant<axis_size_type,Uf+Of>;
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
  using axis_type = typename axis_spec<I>::axis;
  template <unsigned I>
  using edge_type = typename std::decay_t<axis_type<I>>::edge_type;
  using container_type = Container;
  using filler_type = Filler;
  using value_type = typename container_type::value_type;
  using size_type = ivanp::axis_size_type;
  static constexpr unsigned naxes = sizeof...(Ax);
  using index_array_type = std::array<size_type,naxes>;
  using index_array_cref = const index_array_type&;

  static std::vector<named<binner>> all;

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

  template <size_t I>
  constexpr std::enable_if_t<I!=0,size_type> nbins_total_impl() const noexcept {
    return ( axis<I>().nbins() + axis_spec<I>::nover::value
      ) * nbins_total_impl<I-1>();
  }
  template <size_t I>
  constexpr std::enable_if_t<I==0,size_type> nbins_total_impl() const noexcept {
    return ( axis<0>().nbins() + axis_spec<0>::nover::value );
  }

  template <size_t I=0, typename T, typename... Args>
  inline size_type find_bin_impl(const T& x, const Args&... args) const {
    return axis<I>().find_bin(x)
      + (axis<I>().nbins() - !axis_spec<I>::under::value)
      * find_bin_impl<I+1>(args...);
  }
  template <size_t I=0, typename T>
  inline size_type find_bin_impl(const T& x) const {
    return axis<I>().find_bin(x) - !axis_spec<I>::under::value;
  }
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
             - !axis_spec<naxes-sizeof...(TT)-1>::under::value)
      * index_impl(ii...);
  }
  constexpr size_type index_impl(size_type i) const noexcept { return i; }
  template <size_t... I>
  constexpr size_type index_impl(index_array_cref ia, std::index_sequence<I...>)
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

  binner(const binner& o): _axes(o._axes), _bins(o._bins) { }
  binner(binner&& o): _axes(std::move(o._axes)), _bins(std::move(o._bins)) { }
  binner& operator=(const binner& o) {
    _axes = o._axes;
    _bins = o._bins;
    return *this;
  }
  binner& operator=(binner&& o) {
    _axes = std::move(o._axes);
    _bins = std::move(o._bins);
    return *this;
  }

  template <unsigned I=0>
  constexpr const axis_type<I>& axis() const noexcept {
    return std::get<I>(_axes);
  }

  constexpr size_type nbins_total() const noexcept {
    return nbins_total_impl<naxes-1>();
  }

  inline const container_type& bins() const noexcept { return _bins; }
  inline container_type& bins() noexcept { return _bins; }

  constexpr size_type index(replace_t<size_type,Ax>... ii) const noexcept {
    return index_impl(ii...);
  }
  constexpr size_type index(index_array_cref bin) const noexcept {
    return index_impl(bin,std::make_index_sequence<naxes>());
  }

  inline const value_type& bin(replace_t<size_type,Ax>... ii) const {
    return _bins[index_impl(ii...)];
  }

  template <typename... Args>
  inline size_type fill_bin(size_type bin, Args&&... args) {
    filler_type()(_bins[bin], std::forward<Args>(args)...);
    return bin;
  }
  template <typename... Args>
  inline size_type fill_bin(index_array_cref ia, Args&&... args) {
    const auto bin = index(ia);
    filler_type()(_bins[bin], std::forward<Args>(args)...);
    return bin;
  }

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

  template <typename... Args>
  inline std::enable_if_t<sizeof...(Args)==naxes,size_type>
  fill(const Args&... args) {
    return fill_bin(find_bin_impl(args...));
  }
  template <typename... Args>
  inline std::enable_if_t<(sizeof...(Args)>naxes),size_type>
  fill(const Args&... args) {
    auto tup = std::forward_as_tuple(args...);
    return fill_bin_tuple(
      find_bin_tuple(tup,std::make_index_sequence<naxes>()),
      tup, index_sequence_tail<naxes,sizeof...(Args)>()
    );
  }
  template <typename... Args>
  inline std::enable_if_t<(sizeof...(Args)>=naxes),size_type>
  operator()(const Args&... args) {
    return fill(args...);
  }
};

template <typename B, typename... A, typename C, typename F>
std::vector<named<binner<B,std::tuple<A...>,C,F>>>
binner<B,std::tuple<A...>,C,F>::all;

} // end namespace ivanp

#endif
