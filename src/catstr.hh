// Developed by Ivan Pogrebnyak, MSU

#ifndef IVANP_CATSTR_HH
#define IVANP_CATSTR_HH

#include <string>
#include <sstream>
#include <utility>

namespace ivanp {

namespace detail {

template<typename T>
inline void cat_impl(std::stringstream& ss, T&& t) {
  ss << std::forward<T>(t);
}

template<typename T, typename... TT>
inline void cat_impl(std::stringstream& ss, T&& t, TT&&... tt) {
  ss << t;
  cat_impl(ss,std::forward<TT>(tt)...);
}

}

template<typename... TT>
inline std::string cat(TT&&... tt) {
  std::stringstream ss;
  ivanp::detail::cat_impl(ss,std::forward<TT>(tt)...);
  return ss.str();
}

}

#endif
