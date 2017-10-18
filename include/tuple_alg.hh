#ifndef IVANP_TUPLE_ALG_HH
#define IVANP_TUPLE_ALG_HH

#include <tuple>
#include <utility>

namespace ivanp {
namespace detail {

template <template<typename> typename Pred, typename Tuple>
class indices_of_impl {
  static constexpr size_t size = std::tuple_size<Tuple>::value;
  template <size_t I, size_t... II> struct impl {
    using type = std::conditional_t<
      Pred<std::tuple_element_t<I,Tuple>>::value,
      typename impl<I+1, II..., I>::type,
      typename impl<I+1, II...>::type >;
  };
  template <size_t... II> struct impl<size,II...> {
    using type = std::index_sequence<II...>;
  };
public:
  using type = typename impl<0>::type;
};

/*
template <template<typename> typename Pred, typename Tuple>
class first_index_of_impl {
  static constexpr size_t size = std::tuple_size<Tuple>::value;
  template <size_t I, int...> struct impl {
    using type = std::conditional_t<
      Pred<std::tuple_element_t<I,Tuple>>::value,
      std::index_sequence<I>,
      typename impl<I+1>::type >;
  };
  template <int..._> struct impl<size,_...> {
    using type = std::index_sequence<>;
  };
public:
  using type = typename impl<0>::type;
};
*/

} // end namespace detail

template <template<typename> typename Pred, typename Tuple>
using indices_of = typename detail::indices_of_impl<Pred,Tuple>::type;

// template <template<typename> typename Pred, typename Tuple>
// using first_index_of = typename detail::first_index_of_impl<Pred,Tuple>::type;

} // end namespace ivanp

#endif
