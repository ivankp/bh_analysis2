// Written by Ivan Pogrebnyak

#ifndef IVANP_MATH_HH
#define IVANP_MATH_HH

#include <cmath>

namespace ivanp { namespace math {

template <typename T> [[ gnu::const ]]
constexpr auto sq(T x) noexcept { return x*x; }
template <typename T, typename... TT> [[ gnu::const ]]
constexpr auto sq(T x, TT... xx) noexcept { return sq(x)+sq(xx...); }

template <typename... TT> [[ gnu::const ]]
constexpr auto qadd(TT... xx) noexcept { return std::sqrt(sq(xx...)); }

constexpr auto prod() noexcept { return 1; }
template <typename T> [[ gnu::const ]]
constexpr T prod(T x) noexcept { return x; }
template <typename T, typename... TT> [[ gnu::const ]]
constexpr auto prod(T x, TT... xx) noexcept { return x*prod(xx...); }

constexpr auto sum() noexcept { return 0; }
template <typename T> [[ gnu::const ]]
constexpr T sum(T x) noexcept { return x; }
template <typename T, typename... TT> [[ gnu::const ]]
constexpr auto sum(T x, TT... xx) noexcept { return x+sum(xx...); }

constexpr bool eq() noexcept { return true; }
template <typename T> [[ gnu::const ]]
constexpr bool eq(const T&) noexcept { return true; }
template <typename T1, typename T>
constexpr bool eq(const T& x1, const T& x) noexcept { return x1==x; }
template <typename T1, typename T, typename... TT>
constexpr bool eq(const T& x1, const T& x, const TT&... xx) noexcept {
  return (x1==x) && eq(x1,xx...);
}

// return absolute value of phi separation
[[ gnu::const ]]
inline double dphi(double phi1, double phi2) noexcept {
  static constexpr double twopi = M_PI*2;
  double _dphi = phi1 - phi2;
  if (__builtin_expect(_dphi < 0.,0)) _dphi = -_dphi;
  return ( __builtin_expect(_dphi > M_PI,0) ? twopi-_dphi : _dphi );
}

[[ gnu::const ]]
inline double deltaR(double eta1, double eta2, double phi1, double phi2)
noexcept {
  return std::sqrt(sq(eta1-eta2,dphi(phi1,phi2)));
}

template <typename J>
inline double tau(const J& jet, double higgs_y) noexcept {
  return std::sqrt( sq(jet[3])-sq(jet[2]) )/( 2.*std::cosh(jet.rap()-higgs_y) );
}
template <typename P>
inline double pTt(const P& a, const P& b) noexcept {
  return std::abs( a[0]*b[1]-b[0]*a[1] )/( std::sqrt(sq(a[0]+b[0],a[1]+b[1]))*2 );
}

[[ gnu::const ]]
inline double dphi_signed(double phi1, double phi2, double rap1, double rap2)
noexcept {
  static constexpr double twopi = M_PI*2;
  double _dphi = phi1 - phi2;
  if (rap1 < rap2) _dphi = -_dphi;
  while (_dphi >= M_PI) _dphi -= twopi;
  while (_dphi < -M_PI) _dphi += twopi;
  return _dphi;
}

template <typename T1, typename T2>
inline void smaller(T1& x, const T2& y) noexcept { if (y < x) x = y; }
template <typename T1, typename T2>
inline void larger (T1& x, const T2& y) noexcept { if (x < y) x = y; }

[[ gnu::const ]]
constexpr size_t utn(size_t n) noexcept { return n*(n+1)/2; }

}}

#endif
