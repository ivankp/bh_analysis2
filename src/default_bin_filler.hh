// Written by Ivan Pogrebnyak

#ifndef IVANP_DEFAULT_BIN_FILLER_HH
#define IVANP_DEFAULT_BIN_FILLER_HH

#include "type_traits.hh"

namespace ivanp {

template <typename BinType> struct default_bin_filler {

  template <typename Bin = BinType>
  inline typename std::enable_if<
    has_pre_increment<Bin>::value
  >::type
  operator()(Bin& bin) noexcept(noexcept(++bin)) { ++bin; }

  template <typename Bin = BinType>
  inline typename std::enable_if<
    !has_pre_increment<Bin>::value &&
    has_post_increment<Bin>::value
  >::type
  operator()(Bin& bin) noexcept(noexcept(bin++)) { bin++; }

  template <typename Bin = BinType>
  inline typename std::enable_if<
    !has_pre_increment<Bin>::value &&
    !has_post_increment<Bin>::value &&
    is_callable<Bin>::value
  >::type
  operator()(Bin& bin) noexcept(noexcept(bin())) { bin(); }

  template <typename T, typename Bin = BinType>
  inline typename std::enable_if<
    has_plus_eq<Bin,T>::value
  >::type
  operator()(Bin& bin, T&& x) noexcept(noexcept(bin+=std::forward<T>(x)))
  { bin+=std::forward<T>(x); }

  template <typename T, typename Bin = BinType>
  inline typename std::enable_if<
    !has_plus_eq<Bin,T>::value &&
    is_callable<Bin,T>::value
  >::type
  operator()(Bin& bin, T&& x) noexcept(noexcept(bin(x))) { bin(x); }

  template <typename T1, typename... TT, typename Bin = BinType>
  inline typename std::enable_if<
    is_callable<Bin,T1,TT...>::value
  >::type
  operator()(Bin& bin, T1&& arg1, TT&&... args)
  noexcept(noexcept(bin(std::forward<T1>(arg1), std::forward<TT>(args)...)))
  { bin(std::forward<T1>(arg1), std::forward<TT>(args)...); }

};

} // end namespace ivanp

#endif
