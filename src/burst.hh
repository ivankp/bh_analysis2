#ifndef IVANP_BINNER_BURST_HH
#define IVANP_BINNER_BURST_HH

#include <tuple>
#include <array>
#include <utility>

#include "axis.hh"
#include "type_traits.hh"

namespace ivanp {

namespace detail {

template <typename... A, size_t... I>
inline std::array<axis_size_type,sizeof...(I)>
all_nbins(const std::tuple<A...>& axes, std::index_sequence<I...>) {
  return {(std::get<I>(axes).nbins())...};
}

template <typename... A, size_t... I>
inline auto make_ranges(
  const std::tuple<A...>& axes,
  std::index_sequence<I...>
) -> std::tuple<std::array<
      typename std::tuple_element_t<I,std::tuple<A...>>::edge_type,2>...>
{
  return { { std::get<I>(axes).min(), std::get<I>(axes).max() }... };
}

template <typename... A, size_t... I>
inline auto make_slices(
  const std::tuple<A...>& axes,
  const std::array<axis_size_type,sizeof...(I)>& ii,
  std::index_sequence<I...>
) -> std::tuple<std::array<
      typename std::tuple_element_t<I,std::tuple<A...>>::edge_type,2>...>
{
  return {
    { std::get<I>(axes).lower(std::get<I-sizeof...(A)+sizeof...(I)>(ii)+1),
      std::get<I>(axes).upper(std::get<I-sizeof...(A)+sizeof...(I)>(ii)+1)
    }... };
}

template <size_t I=0, size_t N>
inline std::enable_if_t<(I==N-1),std::pair<bool,bool>>
advance_cnt(
  const std::array<axis_size_type,N>& nn,
  std::array<axis_size_type,N>& ii
) noexcept {
  if ((++std::get<I>(ii)) == std::get<I>(nn)) {
    std::get<I>(ii) = 0;
    return {true,true};
  }
  return {false,false};
}
template <size_t I=0, size_t N>
inline std::enable_if_t<(I<N-1),std::pair<bool,bool>>
advance_cnt(
  const std::array<axis_size_type,N>& nn,
  std::array<axis_size_type,N>& ii
) noexcept {
  if ((++std::get<I>(ii)) == std::get<I>(nn)) {
    std::get<I>(ii) = 0;
    return {true,advance_cnt<I+1>(nn,ii).second};
  }
  return {false,false};
}

template <size_t I, typename T, size_t N, typename F>
inline std::enable_if_t<(I==1),T>
prod(const std::array<T,N>& a, F f) {
  return f(std::get<0>(a));
}
template <size_t I, typename T, size_t N, typename F>
inline std::enable_if_t<(I>1 && I<=N),T>
prod(const std::array<T,N>& a, F f) {
  return f(std::get<I-1>(a)) * prod<I-1>(a,f);
}

} // end namespace detail

template <typename, typename, typename> struct binner_slice { };
template <typename It, typename... Ranges, typename... Slices>
struct binner_slice<It,std::tuple<Ranges...>,std::tuple<Slices...>> {
  using ranges_t = std::tuple<std::array<Ranges,2>...>;
  using slices_t = std::tuple<std::array<Slices,2>...>;
  static constexpr size_t ranges_size = sizeof...(Ranges);
  static constexpr size_t slices_size = sizeof...(Slices);
  It begin, end;
  std::array<axis_size_type,ranges_size> nbins;
  ranges_t ranges;
  slices_t slices;
};

template <typename, typename, typename, typename> struct make_binner_slice { };
template <typename It, typename... A, size_t... R, size_t... S>
struct make_binner_slice<It,std::tuple<A...>,
  std::index_sequence<R...>, std::index_sequence<S...>
> {
  using tuple_type = std::tuple<A...>;
  using type = binner_slice<It,
    std::tuple<typename std::tuple_element_t<R,tuple_type>::edge_type...>,
    std::tuple<typename std::tuple_element_t<S,tuple_type>::edge_type...>
  >;
};

// template <typename... T> struct test_type;

template <size_t D, typename... A, typename It>
auto burst(const std::tuple<A...>& axes, It first, It last) {
  static_assert(0<D && D<sizeof...(A),"");
  static_assert(D==1||D==2,"for now, burst only works for D==1,2");
  using seq  = std::make_index_sequence<sizeof...(A)>;
  using seq1 = std::make_index_sequence<D>;
  using seq2 = index_sequence_tail<D,sizeof...(A)>;
  using slice_t = typename make_binner_slice<It,
    std::tuple<A...>, seq1, seq2 >::type;

  const auto ranges = detail::make_ranges(axes,seq1{});

  const auto nbins  = detail::all_nbins(axes,seq{});
  const auto nbins1 = detail::all_nbins(axes,seq1{});
  const auto nbins2 = detail::all_nbins(axes,seq2{});
  std::remove_const_t<decltype(nbins2)> cnt{};

  std::vector<slice_t> slices;
  slices.reserve(
    detail::prod<seq2::size()>(nbins2,[](auto n){ return n; }) );

  const auto len1 = detail::prod<D>(nbins,[](auto n){ return n+2; });
  It it = first + detail::prod<seq::size()-1>(nbins,[](auto n){ return n+2; });
  if (seq2::size()>1) it += len1;

  for ( std::pair<bool,bool> check{false,false}; !check.second;
        check = detail::advance_cnt(nbins2,cnt)
  ) {
    if (check.first) it += len1*2;
    slices.push_back({ it, it+len1, nbins1,
      ranges, detail::make_slices(axes,cnt,seq2{}) });
    it += len1;
  }

  return slices;
}

} // end namespace ivanp

#endif