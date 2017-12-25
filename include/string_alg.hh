// Written by Ivan Pogrebnyak

#ifndef IVANP_STRING_ALGORITHMS_HH
#define IVANP_STRING_ALGORITHMS_HH

#include <string>
#include <array>

namespace ivanp {

template <typename Str, unsigned N>
inline bool starts_with(const Str& str, const char(&prefix)[N]) {
  for (unsigned i=0; i<N-1; ++i)
    if (!str[i] || str[i]!=prefix[i]) return false;
  return true;
}

template <size_t N, typename Str>
std::array<Str,N+1> lsplit(const Str& str, typename Str::value_type delim) {
  std::array<Str,N+1> arr;
  auto l = 0;
  size_t i = 0;
  for (; i!=N; ++i) {
    const auto r = str.find(delim,l);
    if (!(r+1)) break;
    arr[i] = str.substr(l,r-l);
    l = r+1;
  }
  arr[i] = str.substr(l);
  return arr;
}

template <size_t N, typename Str>
std::array<Str,N+1> rsplit(const Str& str, typename Str::value_type delim) {
  std::array<Str,N+1> arr;
  auto r = Str::npos-1;
  size_t i = N;
  for (; i!=0; --i) {
    const auto l = str.rfind(delim,r);
    if (!(l+1)) break;
    arr[i] = str.substr(l+1,r-l);
    r = l-1;
  }
  arr[i] = str.substr(0,r+1);
  return arr;
}

}

#endif
