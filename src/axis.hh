// Written by Ivan Pogrebnyak

#ifndef IVANP_AXIS_HH
#define IVANP_AXIS_HH

#include <algorithm>
#include <cmath>
#include <utility>
#include <stdexcept>
#include <sstream>
#include <memory>

#include "type_traits.hh"

#ifdef _CN
#error In __FILE__: cannot define macro _CN
#endif
#define _CN const noexcept
#ifdef _CF
#error In __FILE__: cannot define macro _CF
#endif
#define _CF const final
#ifdef _CNF
#error In __FILE__: cannot define macro _CNF
#endif
#define _CNF const noexcept final
#ifdef _CNFN
#error In __FILE__: cannot define macro _CNFN
#endif
#define _CNFN(expr) const noexcept(noexcept(expr)) final
#ifdef _CNN
#error In __FILE__: cannot define macro _CNN
#endif
#define _CNN(expr) const noexcept(noexcept(expr))

namespace ivanp {

// INFO =============================================================

/*
 * Same convention as in ROOT TH1:
 * bin = 0;       underflow bin
 * bin = 1;       first bin with low-edge xlow INCLUDED
 * bin = nbins;   last bin with upper-edge xup EXCLUDED
 * bin = nbins+1; overflow bin
 */

// Abstract Axis ====================================================

using axis_size_type = unsigned;

template <typename EdgeType>
class abstract_axis {
public:
  using size_type = axis_size_type;
  using edge_type = EdgeType;
  using edge_cref = const_ref_if_not_scalar_t<edge_type>;

  virtual size_type nbins () const = 0;
  virtual size_type nedges() const = 0;

  virtual size_type vfind_bin(edge_cref x) const = 0;

  virtual edge_cref edge(size_type i) const = 0;
  virtual edge_cref min() const = 0;
  virtual edge_cref max() const = 0;
  virtual edge_cref lower(size_type bin) const = 0;
  virtual edge_cref upper(size_type bin) const = 0;

  inline bool check_edge(size_type i) const { return (i < nedges()); }
  inline void check_edge_throw(size_type i) const {
    if (!check_edge(i)) {
      std::ostringstream ss("Axis edge index ");
      ss << i << " >= " << nedges();
      throw std::range_error(ss.str());
    }
  }
  inline bool check_bin(size_type bin) const { return (bin < nbins()); }
  inline void check_bin_throw(size_type bin) const {
    if (!check_bin(bin)) {
      std::ostringstream ss("Axis bin index ");
      ss << bin << " >= " << nbins();
      throw std::range_error(ss.str());
    }
  }
};

// Container Axis ===================================================

template <typename Container>
class container_axis: public abstract_axis<typename Container::value_type> {
public:
  using base_type = abstract_axis<typename Container::value_type>;
  using size_type = typename base_type::size_type;
  using edge_type = typename base_type::edge_type;
  using edge_cref = typename base_type::edge_cref;
  using container_type = Container;

private:
  container_type _edges;

public:
  container_axis() = default;
  ~container_axis() = default;

  container_axis(const container_type& edges): _edges(edges) { }
  container_axis(container_type&& edges): _edges(std::move(edges)) { }
  container_axis(const container_axis& axis): _edges(axis._edges) { }
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
  container_axis& operator=(container_type&& edges) {
    _edges = std::move(edges);
    return *this;
  }
  container_axis& operator=(const container_axis& axis) {
    _edges = axis._edges;
    return *this;
  }
  container_axis& operator=(container_axis&& axis) {
    _edges = std::move(axis._edges);
    return *this;
  }

  inline size_type nedges() _CNFN(_edges.size()) { return _edges.size(); }
  inline size_type nbins () _CNFN(_edges.size()) { return _edges.size()-1; }

  inline edge_cref edge(size_type i) _CNF { return _edges[i]; }

  inline edge_cref min() _CNN(_edges.front()) { return _edges.front(); }
  inline edge_cref max() _CNN(_edges.back()) { return _edges.back(); }

  inline edge_cref lower(size_type bin) _CNF { return _edges[bin-1]; }
  inline edge_cref upper(size_type bin) _CNF { return _edges[bin]; }

  template <typename T>
  size_type find_bin(const T& x) const noexcept {
    return std::distance(
      _edges.begin(), std::upper_bound(_edges.begin(), _edges.end(), x)
    );
  }
  inline size_type vfind_bin(edge_cref x) _CF { return find_bin(x); }

  template <typename T>
  inline size_type operator[](const T& x) const noexcept {
    return find_bin(x);
  }

  inline const container_type& edges() const { return _edges; }

};

// Uniform Axis =====================================================

template <typename EdgeType>
class uniform_axis: public abstract_axis<EdgeType> {
public:
  using base_type = abstract_axis<EdgeType>;
  using size_type = typename base_type::size_type;
  using edge_type = typename base_type::edge_type;
  using edge_cref = typename base_type::edge_cref;

private:
  size_type _nbins;
  edge_type _min, _max;

public:
  uniform_axis() = default;
  ~uniform_axis() = default;
  uniform_axis(size_type nbins, edge_cref min, edge_cref max)
  : _nbins(nbins), _min(std::min(min,max)), _max(std::max(min,max)) { }
  uniform_axis(const uniform_axis& axis)
  : _nbins(axis._nbins), _min(axis._min), _max(axis._max) { }
  uniform_axis& operator=(const uniform_axis& axis) {
    _nbins = axis._nbins;
    _min = axis._min;
    _max = axis._max;
    return *this;
  }

  inline size_type nbins () _CNF { return _nbins; }
  inline size_type nedges() _CNF { return _nbins+1; }

  inline edge_cref edge(size_type i) _CNF {
    const auto width = (_max - _min)/_nbins;
    return _min + i*width;
  }

  inline edge_cref min() const noexcept { return _min; }
  inline edge_cref max() const noexcept { return _max; }

  inline edge_cref lower(size_type bin) _CNF { return edge(bin-1); }
  inline edge_cref upper(size_type bin) _CNF { return edge(bin); }

  template <typename T>
  size_type find_bin(const T& x) const noexcept {
    if (x < _min) return 0;
    if (!(x < _max)) return _nbins+1;
    return _nbins*(x-_min)/(_max-_min) + 1;
  }

  inline size_type vfind_bin(edge_cref x) _CNF { return find_bin(x); }

  template <typename T>
  inline size_type operator[](const T& x) const noexcept {
    return find_bin(x);
  }

  // exact edge calculation
  edge_type exact_edge(size_type i) const noexcept {
    auto x = edge(i);
    const auto j = find_bin(x);
    const auto i1 = i + 1;

    if (j < i1) {
      for ( ; find_bin(x = std::nextafter(x,x+1)) < i1; );
    } else {
      auto x2 = x;
      for ( ; find_bin(x2 = std::nextafter(x2,x2-1)) == i1; ) x = x2;
    }
    return x;
  }
  inline edge_cref exact_lower(size_type bin) const noexcept {
    return exact_edge(bin-1);
  }
  inline edge_cref exact_upper(size_type bin) const noexcept {
    return exact_edge(bin);
  }

};

// Indirect Axis ====================================================

template <typename EdgeType, typename Ref = const abstract_axis<EdgeType>*>
class ref_axis: public abstract_axis<EdgeType> {
public:
  using base_type = abstract_axis<EdgeType>;
  using size_type = typename base_type::size_type;
  using edge_type = typename base_type::edge_type;
  using edge_cref = typename base_type::edge_cref;
  using axis_ref  = Ref;

private:
  axis_ref _ref;

public:
  ref_axis() = default;
  ~ref_axis() = default;

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

  inline size_type nedges() _CF { return _ref->nedges(); }
  inline size_type nbins () _CF { return _ref->nbins (); }

  inline size_type vfind_bin (edge_cref x) _CF   { return _ref->vfind_bin(x); }
  inline size_type  find_bin (edge_cref x) const { return _ref->vfind_bin(x); }
  inline size_type operator[](edge_cref x) const { return _ref->vfind_bin(x); }

  inline edge_cref edge(size_type i) _CF { return _ref->edge(i); }
  inline edge_cref min() _CF { return _ref->min(); }
  inline edge_cref max() _CF { return _ref->max(); }
  inline edge_cref lower(size_type bin) _CF { return _ref->lower(bin); }
  inline edge_cref upper(size_type bin) _CF { return _ref->upper(bin); }

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

template <typename EdgeType>
class const_axis: public abstract_axis<EdgeType> {
public:
  using base_type = abstract_axis<EdgeType>;
  using size_type = typename base_type::size_type;
  using edge_type = typename base_type::edge_type;
  using edge_cref = typename base_type::edge_cref;

private:
  const edge_type* _edges;
  size_type _ne;

public:
  template <size_type N>
  constexpr const_axis(const edge_type(&a)[N]): _edges(a), _ne(N - 1) {}

  constexpr size_type nedges() _CNF { return _ne+1; }
  constexpr size_type nbins () _CNF { return _ne; }

  constexpr edge_cref edge(size_type i) _CNF { return _edges[i]; }

  constexpr edge_cref min() _CNF { return _edges[0]; }
  constexpr edge_cref max() _CNF { return _edges[_ne]; }

  constexpr edge_cref lower(size_type bin) _CNF { return _edges[bin-1]; }
  constexpr edge_cref upper(size_type bin) _CNF { return _edges[bin]; }

  constexpr size_type find_bin(edge_cref x) _CN {
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
  constexpr size_type operator[](edge_cref x) _CN { return find_bin(x); }
  inline size_type vfind_bin(edge_cref x) _CNF { return find_bin(x); }
};

// ==================================================================
} // end namespace ivanp

#undef _CF
#undef _CN
#undef _CNF
#undef _CNN
#undef _CNFN

#endif
