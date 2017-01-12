// Written by Ivan Pogrebnyak

#ifndef IVANP_MATH_HH
#define IVANP_MATH_HH

#include <cmath>

namespace ivanp {

template <typename T> [[ gnu::const ]]
constexpr T sq(T x) noexcept { return x*x; }
template <typename T, typename... TT> [[ gnu::const ]]
constexpr T sq(T x, TT... xx) noexcept { return sq(x)+sq(xx...); }

// return absolute value of phi separation
template <typename T>
[[ gnu::const ]] inline T dphi(T phi1, T phi2) noexcept {
  T _dphi = phi1 - phi2;
  if (__builtin_expect(_dphi < 0.,0)) _dphi = -_dphi;
  return ( __builtin_expect(_dphi > M_PI,0) ? M_2_PI-_dphi : _dphi );
}

}

#endif
