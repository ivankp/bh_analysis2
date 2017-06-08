#ifndef IVANP_COMPREHENSION_HH
#define IVANP_COMPREHENSION_HH

#ifdef _GLIBCXX_VECTOR
template <typename T, typename Pred>
auto operator|(const std::vector<T>& in, Pred f) {
  std::vector<std::decay_t<decltype(f(std::declval<const T&>()))>> out;
  out.reserve(in.size());
  for (const auto& x : in) out.emplace_back(f(x));
  return std::move(out);
}
#endif

#ifdef _GLIBCXX_ARRAY
template <typename T, size_t N, typename Pred, size_t... I>
inline auto array_transform(const std::array<T,N>& in, Pred f,
                            std::index_sequence<I...>)
-> std::array<decltype(f(std::declval<T>())),N> {
  return { f(std::get<I>(in))... };
}
template <typename T, size_t N, typename Pred>
inline auto operator|(const std::array<T,N>& in, Pred f) {
  return array_transform(in,f,std::make_index_sequence<N>{});
}
#endif

#ifdef _GLIBCXX_TUPLE
template <typename... T, typename Pred, size_t... I>
inline auto tuple_transform(const std::tuple<T...>& in, Pred f,
                            std::index_sequence<I...>
) {
  return std::make_tuple(f(std::get<I>(in))...);
}
template <typename... T, typename Pred>
inline auto operator|(const std::tuple<T...>& in, Pred f) {
  return tuple_transform(in,f,std::index_sequence_for<T...>{});
}
#endif

#endif
