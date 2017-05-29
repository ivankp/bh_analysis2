// Written by Ivan Pogrebnyak

#ifndef IVANP_MATH_HH
#define IVANP_MATH_HH

#include <cmath>

namespace ivanp { namespace math {

template <typename T> [[ gnu::const ]]
constexpr T sq(T x) noexcept { return x*x; }
template <typename T, typename... TT> [[ gnu::const ]]
constexpr T sq(T x, TT... xx) noexcept { return sq(x)+sq(xx...); }

constexpr auto prod() noexcept { return 1; }
template <typename T> [[ gnu::const ]]
constexpr T prod(T x) noexcept { return x; }
template <typename T, typename... TT> [[ gnu::const ]]
constexpr T prod(T x, TT... xx) noexcept { return x*prod(xx...); }

// return absolute value of phi separation
template <typename T> [[ gnu::const ]]
inline T dphi(T phi1, T phi2) noexcept {
  static constexpr T twopi = M_PI*2;

  T _dphi = phi1 - phi2;
  if (__builtin_expect(_dphi < 0.,0)) _dphi = -_dphi;
  return ( __builtin_expect(_dphi > M_PI,0) ? twopi-_dphi : _dphi );
}

template <typename T> [[ gnu::const ]]
inline T deltaR(T eta1, T eta2, T phi1, T phi2) noexcept {
  return std::sqrt(sq(eta1-eta2,dphi(phi1,phi2)));
}

template <typename T1, typename T2>
inline void smaller(T1& x, const T2& y) noexcept { if (y < x) x = y; }
template <typename T1, typename T2>
inline void larger (T1& x, const T2& y) noexcept { if (x < y) x = y; }

[[ gnu::const ]]
constexpr size_t nut(size_t n) noexcept { return n*(n+1)/2; }

}}

#endif
