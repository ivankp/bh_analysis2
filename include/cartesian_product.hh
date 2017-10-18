#ifndef IVANP_CARTESIAN_PRODUCT_HH
#define IVANP_CARTESIAN_PRODUCT_HH

#include <utility>
#include <tuple>

namespace ivanp {

template <typename B>
constexpr bool increment_ranges(const B&) { return true; }
template <typename T, typename B>
bool increment_ranges(const B& begins, std::pair<T,T>& r) {
  ++r.first;
  if (r.first == r.second) return true;
  return false;
}
template <typename T, typename... TT, typename B>
bool increment_ranges(const B& begins, std::pair<T,T>& r, std::pair<TT,TT>&... rr) {
  ++r.first;
  if (r.first == r.second) {
    r.first = std::get<std::tuple_size<B>::value-sizeof...(rr)-1>(begins);
    return increment_ranges(begins,rr...);
  }
  return false;
}

template <typename F, typename... R>
void cartesian_product(F&& f, std::pair<R,R>... ranges) {
  const auto begins = std::make_tuple(ranges.first...);
  for (;;) {
    f(ranges.first...); // apply
    if (increment_ranges(begins, ranges...)) break;
  }
}

}

#endif
