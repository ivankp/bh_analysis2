// Developed by Ivan Pogrebnyak, MSU

#ifndef IVANP_CATSTR_HH
#define IVANP_CATSTR_HH

#include <string>
#include <sstream>
#include <utility>

namespace ivanp {

namespace detail {

template <typename S, typename T>
inline void cat_impl(S& s, const T& t) { s << t; }

template <typename S, typename T, typename... TT>
inline void cat_impl(S& s, const T& t, const TT&... tt) {
  s << t;
  cat_impl(s,tt...);
}

}

template <typename... TT>
inline std::string cat(const TT&... tt) {
  std::stringstream ss;
  detail::cat_impl(ss,tt...);
  return ss.str();
}

}

#endif
