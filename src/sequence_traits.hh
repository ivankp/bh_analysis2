#ifndef IVANP_SEQUENCE_TRAITS_HH
#define IVANP_SEQUENCE_TRAITS_HH

namespace ivanp {
namespace seq {

// is_integer_sequence **********************************************

template <typename> struct is_integer_sequence: std::false_type { };
template <typename T, T... Ints>
struct is_integer_sequence<std::integer_sequence<T,Ints...>>: std::true_type { };

// increment ********************************************************

template <typename Seq, typename Seq::value_type Inc> struct increment;
template <typename T, T... I, T Inc>
struct increment<std::integer_sequence<T,I...>, Inc> {
  using type = std::integer_sequence<T,(I+Inc)...>;
};

// range ************************************************************

template <typename T, size_t A, size_t B>
using make_integer_range = typename
  increment<std::make_integer_sequence<T,B-A>,A>::type;
template <size_t A, size_t B>
using make_index_range = typename
  increment<std::make_index_sequence<B-A>,A>::type;

// subsequence ******************************************************

template <size_t N, typename Seq>
struct rsubseq {
  template <size_t _N, typename _Seq> struct impl;
  template <size_t _N, typename T, T I1, T... II>
  struct impl<_N,std::integer_sequence<T,I1,II...>> {
    using type = typename
      rsubseq<_N-1,std::integer_sequence<T,II...>>::type;
  };
  using type = typename impl<N,Seq>::type;
};
template <typename Seq> struct rsubseq<0,Seq> { using type = Seq; };
template <size_t N, typename Seq>
using rsubseq_t = typename rsubseq<N,Seq>::type;

template <size_t N, typename Seq>
struct lsubseq {
  template <size_t _N, typename _Seq> struct impl;
  template <size_t _N, typename T, T I1, T... II>
  struct impl<_N,std::integer_sequence<T,II...,I1>> {
    using type = typename
      lsubseq<_N-1,std::integer_sequence<T,II...>>::type;
  };
  using type = typename impl<N,Seq>::type;
};
template <typename Seq> struct lsubseq<0,Seq> { using type = Seq; };
template <size_t N, typename Seq>
using lsubseq_t = typename lsubseq<N,Seq>::type;

// sequence element *************************************************

template <typename Seq> struct front;
template <typename T, T I1, T... II>
struct front<std::integer_sequence<T,I1,II...>> {
  static constexpr T value = I1;
};

template <size_t I, typename Seq>
struct element {
  static constexpr auto value = front<rsubseq_t<I,Seq>>::value;
};

// inverse **********************************************************

template <typename S, typename I
  = std::make_integer_sequence<typename S::value_type, S::size()>
> struct inverse;
template <typename T, T... S, T... I>
class inverse<std::integer_sequence<T,S...>,std::integer_sequence<T,I...>> {
  typedef T arr[sizeof...(S)];

  constexpr static auto inv(T j) {
    arr aseq = { S... };
    arr ainv = { };
    for (T i=0; i<T(sizeof...(S)); ++i) ainv[aseq[i]] = i;
    return ainv[j];
  }

public:
  using type = std::integer_sequence<T,inv(I)...>;
};
template <typename S> using inverse_t = typename inverse<S>::type;

// ******************************************************************

} // end namespace seq
} // end namespace ivanp

#endif
