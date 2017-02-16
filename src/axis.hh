// Written by Ivan Pogrebnyak

#ifndef IVANP_AXIS_HH
#define IVANP_AXIS_HH

#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>
#include <memory>
#include <limits>

#include "type_traits.hh"
#include "exception.hh"

namespace ivanp {

using axis_size_type = unsigned;

template <typename T>
class edge_proxy {
public:
  enum state : char { minf = -1, ok = 0, pinf = 1 };
  using type = T;

private:
  state s;
  T x;

public:
  edge_proxy(state s): s{s}, x{
    s==pinf ?  std::numeric_limits<T>::infinity() :
    s==minf ? -std::numeric_limits<T>::infinity() : T{ }
  } { }

  template <typename U>
  edge_proxy& operator=(U&& x) {
    s = ok;
    this->x = std::forward<U>(x);
    return *this;
  }

  inline operator bool() const noexcept { return s == ok; }
  template <typename U = T>
  inline operator std::enable_if_t<!std::is_same<U,bool>::value,U>()
    const noexcept { return x; }
  inline T get() const noexcept { return x; }

  friend std::ostream& operator<<(std::ostream& os, const edge_proxy& p) {
    switch (p.s) {
      case   ok: os << p.x;  break;
      case pinf: os <<  "∞"; break;
      case minf: os << "-∞"; break;
    }
    return os;
  }
};


// INFO =============================================================

/*
 * Same convention as in ROOT TH1:
 * bin = 0;       underflow bin
 * bin = 1;       first bin with low-edge xlow INCLUDED
 * bin = nbins;   last bin with upper-edge xup EXCLUDED
 * bin = nbins+1; overflow bin
 */

// Abstract Axis ====================================================

template <typename EdgeType>
class abstract_axis {
public:
  using size_type = axis_size_type;
  using edge_type = EdgeType;
  using edge_ptype = edge_proxy<edge_type>;

  virtual size_type nbins () const = 0;
  virtual size_type nedges() const = 0;

  virtual size_type vfind_bin(edge_type x) const = 0;
  inline  size_type  find_bin(edge_type x) const { return vfind_bin(x); }

  virtual edge_type edge(size_type i) const = 0;
  virtual edge_type min() const = 0;
  virtual edge_type max() const = 0;
  virtual edge_ptype lower(size_type i) const = 0;
  virtual edge_ptype upper(size_type i) const = 0;

  virtual bool is_uniform() const = 0;
};

// Blank axis base ==================================================

struct axis_base { };

// Container Axis ===================================================

template <typename Container, bool Inherit=false>
class container_axis final: public std::conditional_t<Inherit,
  abstract_axis<typename std::decay_t<Container>::value_type>,
  axis_base>
{
public:
  using base_type = std::conditional_t<Inherit,
    abstract_axis<typename std::decay_t<Container>::value_type>,
    axis_base>;
  using container_type = Container;
  using edge_type = typename std::decay_t<container_type>::value_type;
  using edge_ptype = edge_proxy<edge_type>;
  using size_type = ivanp::axis_size_type;

private:
  container_type _edges;

public:
  container_axis() = default;
  ~container_axis() = default;

  container_axis(const container_type& edges): _edges(edges) { }
  template <typename C=container_type,
            std::enable_if_t<!std::is_reference<C>::value>* = nullptr>
  container_axis(container_type&& edges): _edges(std::move(edges)) { }

  container_axis(const container_axis& axis): _edges(axis._edges) { }
  template <typename C=container_type,
            std::enable_if_t<!std::is_reference<C>::value>* = nullptr>
  container_axis(container_axis&& axis): _edges(std::move(axis._edges)) { }

  template <typename T, typename C=container_type,
            std::enable_if_t<!is_std_array<C>::value>* = nullptr>
  container_axis(std::initializer_list<T> edges): _edges(edges) { }
  template <typename T, typename C=container_type,
            std::enable_if_t<is_std_array<C>::value>* = nullptr>
  container_axis(std::initializer_list<T> edges) {
    std::copy(edges.begin(),edges.end(),_edges.begin());
  }

  container_axis& operator=(const container_type& edges) {
    _edges = edges;
    return *this;
  }
  template <typename C=container_type,
            std::enable_if_t<!std::is_reference<C>::value>* = nullptr>
  container_axis& operator=(container_type&& edges) {
    _edges = std::move(edges);
    return *this;
  }

  container_axis& operator=(const container_axis& axis) {
    _edges = axis._edges;
    return *this;
  }
  template <typename C=container_type,
            std::enable_if_t<!std::is_reference<C>::value>* = nullptr>
  container_axis& operator=(container_axis&& axis) {
    _edges = std::move(axis._edges);
    return *this;
  }

  inline size_type nedges() const noexcept(noexcept(_edges.size()))
  { return _edges.size(); }
  inline size_type nbins () const noexcept(noexcept(_edges.size()))
  { return _edges.size()-1; }

  inline edge_type edge(size_type i) const noexcept { return _edges[i]; }

  inline edge_type min() const noexcept(noexcept(_edges.front()))
  { return _edges.front(); }
  inline edge_type max() const noexcept(noexcept(_edges.back()))
  { return _edges.back(); }

  inline edge_ptype lower(size_type i) const noexcept {
    edge_ptype proxy(
      i==0 ? edge_ptype::minf :
      i>nedges()+1 ? edge_ptype::pinf : edge_ptype::ok );
    if (proxy) proxy = edge(i-1);
    return proxy;
  }
  inline edge_ptype upper(size_type i) const noexcept {
    edge_ptype proxy( i>=nedges() ? edge_ptype::pinf : edge_ptype::ok );
    if (proxy) proxy = edge(i);
    return proxy;
  }

  template <typename T>
  size_type find_bin(const T& x) const noexcept {
    return std::distance(
      _edges.begin(), std::upper_bound(_edges.begin(), _edges.end(), x)
    );
  }
  inline size_type vfind_bin(edge_type x) const { return find_bin(x); }

  template <typename T>
  inline size_type operator[](const T& x) const noexcept {
    return find_bin(x);
  }

  inline const container_type& edges() const { return _edges; }

  constexpr bool is_uniform() const noexcept { return false; }

};

// Uniform Axis =====================================================

template <typename EdgeType, bool Inherit=false>
class uniform_axis final: public std::conditional_t<Inherit,
  abstract_axis<EdgeType>, axis_base>
{
public:
  using base_type = std::conditional_t<Inherit,
    abstract_axis<EdgeType>, axis_base>;
  using edge_type = EdgeType;
  using edge_ptype = edge_proxy<edge_type>;
  using size_type = ivanp::axis_size_type;

private:
  size_type _nbins;
  edge_type _min, _max;

public:
  uniform_axis() = default;
  ~uniform_axis() = default;
  uniform_axis(size_type nbins, edge_type min, edge_type max)
  : _nbins(nbins), _min(std::min(min,max)), _max(std::max(min,max)) { }
  uniform_axis(const uniform_axis& axis)
  : _nbins(axis._nbins), _min(axis._min), _max(axis._max) { }
  uniform_axis& operator=(const uniform_axis& axis) {
    _nbins = axis._nbins;
    _min = axis._min;
    _max = axis._max;
    return *this;
  }

  inline size_type nbins () const noexcept { return _nbins; }
  inline size_type nedges() const noexcept { return _nbins+1; }

  inline edge_type edge(size_type i) const noexcept {
    const auto width = (_max - _min)/_nbins;
    return _min + i*width;
  }

  inline edge_type min() const noexcept { return _min; }
  inline edge_type max() const noexcept { return _max; }

  inline edge_ptype lower(size_type i) const noexcept {
    edge_ptype proxy(
      i==0 ? edge_ptype::minf :
      i>nedges()+1 ? edge_ptype::pinf : edge_ptype::ok );
    if (proxy) proxy = edge(i-1);
    return proxy;
  }
  inline edge_ptype upper(size_type i) const noexcept {
    edge_ptype proxy( i>=nedges() ? edge_ptype::pinf : edge_ptype::ok );
    if (proxy) proxy = edge(i);
    return proxy;
  }

  template <typename T>
  size_type find_bin(const T& x) const noexcept {
    if (x < _min) return 0;
    if (!(x < _max)) return _nbins+1;
    return _nbins*(x-_min)/(_max-_min) + 1;
  }

  inline size_type vfind_bin(edge_type x) const noexcept
  { return find_bin(x); }

  template <typename T>
  inline size_type operator[](const T& x) const noexcept
  { return find_bin(x); }

  constexpr bool is_uniform() const noexcept { return true; }

};

// Index Axis =======================================================

template <typename EdgeType = ivanp::axis_size_type, bool Inherit=false>
class index_axis final: public std::conditional_t<Inherit,
  abstract_axis<ivanp::axis_size_type>, axis_base>
{
  static_assert(std::is_integral<EdgeType>::value,"");
public:
  using base_type = std::conditional_t<Inherit,
    abstract_axis<ivanp::axis_size_type>, axis_base>;
  using edge_type = EdgeType;
  using edge_ptype = edge_proxy<edge_type>;
  using size_type = ivanp::axis_size_type;

private:
  edge_type _min, _max;

public:
  index_axis() = default;
  ~index_axis() = default;
  constexpr index_axis(edge_type min, edge_type max)
  : _min(std::min(min,max)), _max(std::max(min,max)) { }
  constexpr index_axis(const index_axis& axis)
  : _min(axis._min), _max(axis._max) { }
  constexpr index_axis& operator=(const index_axis& axis) {
    _min = axis._min;
    _max = axis._max;
    return *this;
  }

  constexpr size_type nbins () const noexcept { return _max-_min; }
  constexpr size_type nedges() const noexcept { return _max-_min+1; }

  constexpr edge_type edge(size_type i) const noexcept { return _min + i; }

  constexpr edge_type min() const noexcept { return _min; }
  constexpr edge_type max() const noexcept { return _max; }

  inline edge_ptype lower(size_type i) const noexcept {
    edge_ptype proxy(
      i==0 ? edge_ptype::minf :
      i>nedges()+1 ? edge_ptype::pinf : edge_ptype::ok );
    if (proxy) proxy = edge(i-1);
    return proxy;
  }
  inline edge_ptype upper(size_type i) const noexcept {
    edge_ptype proxy( i>=nedges() ? edge_ptype::pinf : edge_ptype::ok );
    if (proxy) proxy = edge(i);
    return proxy;
  }

  template <typename T>
  constexpr size_type find_bin(const T& x) const noexcept {
    if (x < _min) return 0;
    if (!(x < _max)) return _max-_min+1;
    return x-_min+1;
  }

  constexpr size_type vfind_bin(edge_type x) const noexcept
  { return find_bin(x); }

  template <typename T>
  constexpr size_type operator[](const T& x) const noexcept
  { return find_bin(x); }

  constexpr bool is_uniform() const noexcept { return true; }

};

// Indirect Axis ====================================================

template <typename EdgeType, typename Ref = const abstract_axis<EdgeType>*,
          bool Inherit = false>
class ref_axis final: public std::conditional_t<Inherit,
  abstract_axis<EdgeType>, axis_base>
{
public:
  using base_type = std::conditional_t<Inherit,
    abstract_axis<EdgeType>, axis_base>;
  using edge_type = EdgeType;
  using edge_ptype = edge_proxy<edge_type>;
  using size_type = ivanp::axis_size_type;
  using axis_ref  = Ref;

private:
  axis_ref _ref;

public:
  ref_axis() = default;
  ~ref_axis() = default;

  ref_axis(abstract_axis<EdgeType>* ptr): _ref(ptr) { }

  ref_axis(axis_ref ref): _ref(ref) { }
  ref_axis& operator=(axis_ref ref) {
    _ref = ref;
    return *this;
  }
  ref_axis(const ref_axis& axis): _ref(axis._ref) { }
  ref_axis(ref_axis&& axis): _ref(std::move(axis._ref)) { }
  ref_axis& operator=(const ref_axis& axis) {
    _ref = axis._ref;
    return *this;
  }
  ref_axis& operator=(ref_axis&& axis) {
    _ref = std::move(axis._ref);
    return *this;
  }

  inline size_type nedges() const { return _ref->nedges(); }
  inline size_type nbins () const { return _ref->nbins (); }

  inline size_type vfind_bin (edge_type x) const { return _ref->vfind_bin(x); }
  inline size_type  find_bin (edge_type x) const { return _ref->vfind_bin(x); }
  inline size_type operator[](edge_type x) const { return _ref->vfind_bin(x); }

  inline edge_type edge(size_type i) const { return _ref->edge(i); }
  inline edge_type min() const { return _ref->min(); }
  inline edge_type max() const { return _ref->max(); }
  inline edge_ptype lower(size_type i) const { return _ref->lower(i); }
  inline edge_ptype upper(size_type i) const { return _ref->upper(i); }

  inline bool is_uniform() const noexcept { return _ref->is_uniform(); }

};

// Factory functions ================================================

template <typename A, typename B, typename EdgeType = std::common_type_t<A,B>>
inline decltype(auto) make_axis(axis_size_type nbins, A min, B max) {
  return uniform_axis<EdgeType>(nbins,min,max);
}

template <typename T, size_t N>
inline decltype(auto) make_axis(const std::array<T,N>& edges) {
  return container_axis<std::array<T,N>>(edges);
}

template <typename EdgeType>
inline decltype(auto) make_unique_axis(const abstract_axis<EdgeType>* axis) {
  return ref_axis<EdgeType,std::unique_ptr<abstract_axis<EdgeType>>>(axis);
}

template <typename EdgeType>
inline decltype(auto) make_shared_axis(const abstract_axis<EdgeType>* axis) {
  return ref_axis<EdgeType,std::shared_ptr<abstract_axis<EdgeType>>>(axis);
}

template <typename A, typename B, typename EdgeType = std::common_type_t<A,B>>
inline decltype(auto) make_unique_axis(axis_size_type nbins, A min, B max) {
  return ref_axis<EdgeType,std::unique_ptr<abstract_axis<EdgeType>>>(
    std::make_unique<uniform_axis<EdgeType>>(nbins,min,max) );
}

template <typename A, typename B, typename EdgeType = std::common_type_t<A,B>>
inline decltype(auto) make_shared_axis(axis_size_type nbins, A min, B max) {
  return ref_axis<EdgeType,std::shared_ptr<abstract_axis<EdgeType>>>(
    std::make_shared<uniform_axis<EdgeType>>(nbins,min,max) );
}

template <typename T, size_t N>
inline decltype(auto) make_unique_axis(const std::array<T,N>& edges) {
  return ref_axis<T,std::unique_ptr<abstract_axis<T>>>(
    std::make_unique<uniform_axis<T>>(edges) );
}

template <typename T, size_t N>
inline decltype(auto) make_shared_axis(const std::array<T,N>& edges) {
  return ref_axis<T,std::shared_ptr<abstract_axis<T>>>(
    std::make_shared<uniform_axis<T>>(edges) );
}

// Constexpr Axis ===================================================

template <typename EdgeType, bool Inherit=false>
class const_axis final: public std::conditional_t<Inherit,
  abstract_axis<EdgeType>, axis_base>
{
public:
  using base_type = std::conditional_t<Inherit,
    abstract_axis<EdgeType>, axis_base>;
  using edge_type = EdgeType;
  using edge_ptype = edge_proxy<edge_type>;
  using size_type = ivanp::axis_size_type;

private:
  const edge_type* _edges;
  size_type _ne;

public:
  template <size_type N>
  constexpr const_axis(const edge_type(&a)[N]): _edges(a), _ne(N - 1) {}

  constexpr size_type nedges() const noexcept { return _ne+1; }
  constexpr size_type nbins () const noexcept { return _ne; }

  constexpr edge_type edge(size_type i) const noexcept { return _edges[i]; }

  constexpr edge_type min() const noexcept { return _edges[0]; }
  constexpr edge_type max() const noexcept { return _edges[_ne]; }

  inline edge_ptype lower(size_type i) const noexcept {
    edge_ptype proxy(
      i==0 ? edge_ptype::minf :
      i>nedges()+1 ? edge_ptype::pinf : edge_ptype::ok );
    if (proxy) proxy = edge(i-1);
    return proxy;
  }
  inline edge_ptype upper(size_type i) const noexcept {
    edge_ptype proxy( i>=nedges() ? edge_ptype::pinf : edge_ptype::ok );
    if (proxy) proxy = edge(i);
    return proxy;
  }

  constexpr size_type find_bin(edge_type x) const noexcept {
    size_type i = 0, j = 0, count = _ne, step = 0;

    if (!(x < _edges[_ne])) i = _ne + 1;
    else if (!(x < _edges[0])) while (count > 0) {
      step = count / 2;
      j = step + i;
      if (!(x < _edges[j])) {
        i = j + 1;
        count -= step + 1;
      } else count = step;
    }
    return i;
  }
  constexpr size_type operator[](edge_type x) const noexcept
  { return find_bin(x); }
  inline size_type vfind_bin(edge_type x) const noexcept
  { return find_bin(x); }

  constexpr bool is_uniform() const noexcept { return false; }

};

// ==================================================================

template <typename T, typename Axis>
auto vector_of_edges(const Axis& axis) {
  const auto n = axis.nedges();
  std::vector<T> edges;
  edges.reserve(n);
  for (typename Axis::size_type i=0; i<n; ++i)
    edges.emplace_back(axis.edge(i));
  return edges;
}

} // end namespace ivanp

#endif
