#ifndef CATEGORY_BIN_HH
#define CATEGORY_BIN_HH

#include "enum_traits.hh"

template <typename Bin, typename E, typename... Es>
struct category_bin {
  static unsigned _id;
  using category = E;

  template <typename...>
  struct next_bin { using type = Bin; };
  template <typename _E, typename... _Es>
  struct next_bin<_E,_Es...> { using type = category_bin<Bin,_E,_Es...>; };

  using bin_type = typename next_bin<Es...>::type;

  std::array< bin_type, enum_traits<E>::size > bins;

  inline category_bin& operator()(double w) noexcept {
    std::get<0>(bins)(w);
    if (_id) bins[_id](w);
    return *this;
  }
  inline category_bin& operator+=(const category_bin& b) noexcept {
    for (auto i=bins.size(); i; ) --i, bins[i] += b.bins[i];
    return *this;
  }
  inline bin_type& operator->() noexcept { return bins[_id]; }
  inline const bin_type& operator->() const noexcept { return bins[_id]; }

  template <typename _E>
  inline static std::enable_if_t<std::is_same<_E,E>::value,unsigned>& id() {
    return _id;
  }
  template <typename _E>
  inline static std::enable_if_t<!std::is_same<_E,E>::value,unsigned>& id() {
    return bin_type::template id<_E>();
  }
};

template <typename Bin, typename E, typename... Es>
unsigned category_bin<Bin,E,Es...>::_id = 0;

#endif
