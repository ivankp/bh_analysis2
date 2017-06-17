#ifndef IVANP_UTILITY_HH
#define IVANP_UTILITY_HH

#include <string>

namespace ivanp {

template <typename T> constexpr T zero(T) { return 0; }

template <typename T> const T& as_const(T& x) { return x; }

template <typename T>
struct named_ptr {
  using type = T;
  type *p;
  std::string name;

  named_ptr(): p(nullptr), name() { }
  named_ptr(const named_ptr& n) = default;
  named_ptr(named_ptr&& n) = default;
  template <typename Name>
  named_ptr(T* ptr, Name&& name): p(ptr), name(std::forward<Name>(name)) { }

  inline type& operator*() const noexcept { return *p; }
  inline type* operator->() const noexcept { return p; }
};

#ifdef _GLIBCXX_VECTOR
template <typename T>
auto reserve(size_t n) {
  std::vector<T> _v;
  _v.reserve(n);
  return _v;
}
#endif

#ifdef _GLIBCXX_ARRAY
template <typename T, typename... Args>
std::array<T,sizeof...(Args)> make_array(Args&&... args) {
  return { std::forward<Args>(args)... };
}
#endif

}

#endif
