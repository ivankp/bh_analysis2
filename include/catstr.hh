// Developed by Ivan Pogrebnyak, MSU

#ifndef IVANP_CATSTR_HH
#define IVANP_CATSTR_HH

#include <string>
#include <sstream>

namespace ivanp {

namespace detail {

template <typename S>
inline void cat_impl(S& s) { }

template <typename S, typename T>
inline void cat_impl(S& s, const T& x) { s << x; }

template <typename S, typename T, typename... TT>
inline void cat_impl(S& s, const T& x, const TT&... xx) {
  s << x;
  cat_impl(s,xx...);
}

}

template <typename... TT>
inline std::string cat(const TT&... xx) {
  std::stringstream ss;
  detail::cat_impl(ss,xx...);
  return ss.str();
}

}

#endif
