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
inline std::enable_if_t<(I<N-1),bool> advance_cnt(
  const std::array<axis_size_type,N>& nn,
  std::array<axis_size_type,N>& ii
) noexcept {
  if ((++std::get<I>(ii)) == std::get<I>(nn)) {
    std::get<I>(ii) = 0;
    return advance_cnt<I+1>(nn,ii);
  }
  return true;
}
template <size_t I=0, size_t N>
inline std::enable_if_t<(I==N-1),bool> advance_cnt(
  const std::array<axis_size_type,N>& nn,
  std::array<axis_size_type,N>& ii
) noexcept {
  if ((++std::get<I>(ii)) == std::get<I>(nn)) {
    std::get<I>(ii) = 0;
    return false;
  }
  return true;
}

}

template <typename It, typename... Edges>
struct binner_slice {
  It begin, end;
  std::tuple<std::array<Edges,2>...> ranges;
};

template <typename It, typename, typename> struct make_binner_slice { };
template <typename It, typename... A, size_t... I>
struct make_binner_slice<It,std::tuple<A...>,std::index_sequence<I...>> {
  using tuple_type = std::tuple<A...>;
  using type = binner_slice<It,
    typename std::tuple_element_t<I,tuple_type>::edge_type... >;
};

// template <typename... T> struct test_type;

template <size_t D, typename... A, typename It>
decltype(auto) burst(const std::tuple<A...>& axes, It first, It last) {
  static_assert(0<D && D<sizeof...(A),"");
  // using seq = std::make_index_sequence_for<A...>;
  using seq = index_sequence_tail<D,sizeof...(A)>;
  using arr = std::array<axis_size_type,seq::size()>;
  // using slice_t = binner_slice<It,typename A::edge_type...>;
  using slice_t = typename make_binner_slice<It,std::tuple<A...>,seq>::type;

  // test_type<decltype( detail::all_nbins(axes,seq{}) )> a;
  // test_type<std::array<int,sizeof(slice_t)>> a;

  arr nbins1 = detail::all_nbins(axes,std::make_index_sequence<D>{});
  const arr nbins2 = detail::all_nbins(axes,seq{});
  arr cnt;

  size_t len1 = 1;
  for (auto& n : nbins1) len1 *= (n + 2);
  size_t len2 = 1;
  for (auto  n : nbins2) len2 *= n;

  std::vector<slice_t> slices;
  slices.reserve(len2);

  It it = first + len1;

  // test_type<decltype( detail::make_ranges(axes,cnt,seq{}) )> a;

  std::cout << "TEST" << std::endl;
  do {
    // for (auto i : cnt) std::cout << ' ' << i;
    // std::cout << std::endl;
    // auto ranges = detail::make_ranges(axes,cnt,seq{});
    // std::cout << std::get<0>(ranges)[0] << ' ';
    // std::cout << std::get<0>(ranges)[1] << std::endl;

    slice_t s{it,it+len1,
      detail::make_ranges(axes,cnt,seq{})};
    slices.emplace_back(s);

    // slices.emplace_back(first,first+len1,
    //   detail::make_ranges(axes,cnt,seq{}));
    it += len1;
  } while (detail::advance_cnt(nbins2,cnt));

  // return std::move(slices);
  return slices;
}

} // end namespace ivanp

#endif
