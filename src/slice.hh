#ifndef IVANP_BINNER_SLICE_HH
#define IVANP_BINNER_SLICE_HH

#include "binner.hh"
#include <memory>
#include <functional>

namespace ivanp {

namespace detail { namespace slice {

template <size_t N>
class counter {
  using size_type = ivanp::axis_size_type;

  std::array<size_type,N> _ii, _nn;
  bool _next, _done;

public:
  template <typename M, size_t... I>
  counter(const M& m, std::index_sequence<I...>) noexcept
  : _ii{zero(I)...}, _nn{std::get<I>(m)...},
    _next(false), _done(!sizeof...(I)) { }

  counter& operator++() noexcept {
    _next = false;
    ++_ii[0];
    for (size_type i=0; i<N; ++i) {
      if (_ii[i] == _nn[i]) {
        _ii[i] = 0;
        _next = true;
        if (i==N-1) _done = true;
        else ++_ii[i+1];
      }
    }
    return *this;
  }

  inline operator bool() const noexcept { return _done; }
  inline bool next() const noexcept { return _next; }
  inline const std::array<size_type,N>& ii() const noexcept { return _ii; }
  size_type size() const noexcept {
    size_type s = 1;
    for (const auto& n : _nn) s *= n;
    return s;
  }
  inline void reset() noexcept { _done = false; }
};

template <size_t I=0, size_t N>
inline std::enable_if_t<I==N-1, ivanp::axis_size_type> index(
  const std::array<const ivanp::axis_size_type*,N>& ii,
  const std::array<ivanp::axis_size_type,N>& nn
) noexcept {
  return *std::get<I>(ii);
}
template <size_t I=0, size_t N>
inline std::enable_if_t<I!=N-1, ivanp::axis_size_type> index(
  const std::array<const ivanp::axis_size_type*,N>& ii,
  const std::array<ivanp::axis_size_type,N>& nn
) noexcept {
  return *std::get<I>(ii) + std::get<I>(nn) * index<I+1>(ii,nn);
}

template <typename T, size_t N1, size_t N2, size_t... I1, size_t... I2>
inline auto cat_cptr(
  const std::array<T,N1>& arr1,
  const std::array<T,N2>& arr2,
  std::index_sequence<I1...>,
  std::index_sequence<I2...>
) -> std::array<const T*, N1+N2> {
  static_assert(sizeof...(I1) == N1,"");
  static_assert(sizeof...(I2) == N2,"");
  return { &std::get<I1>(arr1)..., &std::get<I2>(arr2)... };
}

template <typename T, size_t N, size_t... I>
inline std::array<T,N> sort(std::array<T,N> arr, std::index_sequence<I...>) {
  static_assert(sizeof...(I) == N,"");
  return { std::get<I>(arr)... };
}

} // end namespace slice
} // end namespace detail

template <size_t D, typename Bin, typename Ax> struct binner_slice;
template <size_t D, typename Bin, typename... Ax>
class binner_slice<D, Bin, std::tuple<Ax...>> {
  // types ----------------------------------------------------------
  using head = std::make_index_sequence<D>;
  using tail = seq::make_index_range<D,sizeof...(Ax)>;

  using specs = std::tuple<Ax...>;
  using specs_head = subtuple_t<specs,head>;
  using specs_tail = subtuple_t<specs,tail>;

  using tail_indices = std::array<ivanp::axis_size_type,tail::size()>;
  template <size_t I>
  using under = typename std::tuple_element_t<I,specs>::under;

  template <typename> struct axes_crefs_from_specs;
  template <typename... T> struct axes_crefs_from_specs<std::tuple<T...>> {
    using type = std::tuple<const typename T::axis&...>;
  };
  using axes_crefs = typename axes_crefs_from_specs<specs>::type;

  template <typename> struct edges_from_specs;
  template <typename... T> struct edges_from_specs<std::tuple<T...>> {
    using type = std::tuple<std::array<typename T::axis::edge_ptype,2>...>;
  };

  template <typename> struct cref_axis_specs;
  template <typename... T> struct cref_axis_specs<std::tuple<T...>> {
    using type = std::tuple< axis_spec<
      const typename T::axis&,
      T::under::value, T::over::value, T::excep::value >... >;
  };

public:
  // binner with reference to axes and bins
  using binner_type = binner< Bin,
    typename cref_axis_specs<specs_head>::type,
    std::vector<std::reference_wrapper<const Bin>> >;

private:
  // members --------------------------------------------------------
  typename edges_from_specs<specs_tail>::type _ranges;
  binner_type _slice;

  // constructor impl -----------------------------------------------
  template <size_t... IH, size_t... IT>
  binner_slice( const axes_crefs& axes, const tail_indices& ii,
    typename binner_type::container_type&& bins,
    std::index_sequence<IH...>, std::index_sequence<IT...>
  ) : _ranges{ {
      std::get<IT>(axes).lower( std::get<IT-D>(ii) + !under<IT>::value ),
      std::get<IT>(axes).upper( std::get<IT-D>(ii) + !under<IT>::value )
    }... },
      _slice( std::tie(std::get<IH>(axes)...), std::move(bins) )
  { }

  // name proxy struct ----------------------------------------------
  template <typename... L> struct name_proxy {
    static_assert(sizeof...(L) <= tail::size(),
      "\033[33mcannot have more labels than dimensions\033[0m");
    const binner_slice *ptr;
    std::tuple<L...> labels;

    template <size_t I, typename = void> struct label_type;
    template <size_t I>
    struct label_type<I,std::enable_if_t<(I<sizeof...(L))>> {
      using type = std::tuple_element_t<I,decltype(labels)>;
    };

    template <size_t I>
    std::enable_if_t<!(I<sizeof...(L)),const char*>
    label() const noexcept { return ""; }
    template <size_t I>
    std::enable_if_t<(I<sizeof...(L)),const typename label_type<I>::type&>
    label() const noexcept { return std::get<I>(labels); }

    template <size_t I=0>
    inline std::enable_if_t<I==tail::size(), std::ostream>&
    impl(std::ostream& os) const { return os; }
    template <size_t I=0>
    inline std::enable_if_t<I!=tail::size(), std::ostream>&
    impl(std::ostream& os) const {
      os << '_' << label<I>()
         << '[' << std::get<I>(ptr->_ranges)[0]
         << ',' << std::get<I>(ptr->_ranges)[1] << ')';
      return impl<I+1>(os);
    }

    friend std::ostream& operator<<( std::ostream& os, const name_proxy& p)
    { return p.impl(os); }
  };
  // ----------------------------------------------------------------

public:
  inline binner_slice(
    const axes_crefs& axes, const tail_indices& ii,
    typename binner_type::container_type&& bins
  ) : binner_slice(axes,ii,std::move(bins),head{},tail{}) { }

  inline auto& operator* () const noexcept { return  _slice; }
  inline auto* operator->() const noexcept { return &_slice; }

  // name proxy factory function ------------------------------------
  template <typename... L>
  inline auto name(L&&... labels) const noexcept
  -> std::enable_if_t< !pack_is_tuple<L...>::value,
    name_proxy<add_const_to_ref_t<L>...>
  > {
    return { this, std::forward_as_tuple(std::forward<L>(labels)...) };
  }
  template <typename... L>
  inline auto name(const std::tuple<L...>& labels) const noexcept
  -> name_proxy<remove_rref_t<add_const_to_ref_t<L>>...> {
    return { this, labels };
  }
  // ----------------------------------------------------------------
};

template <size_t D=1, // slicing into chunks of D dimensions
          typename Bin,
          typename... Ax,
          typename Container,
          typename Filler,
          size_t... I
>
auto slice(
  const binner<Bin,std::tuple<Ax...>,Container,Filler>& hist,
  std::index_sequence<I...>
) {
  static_assert( sizeof...(I) == sizeof...(Ax),
    "\033[33mthe number of indices must match the number of axes\033[0m");
  static_assert( D > 0,
    "\033[33mmust leave at least 1 dimension unsliced\033[0m");
  static_assert( D < sizeof...(I),
    "\033[33mnumber of unsliced dimensions must be less than total\033[0m");
  using namespace ivanp::detail::slice;

  using specs = std::tuple<Ax...>;
  using head = std::make_index_sequence<D>;
  using tail = seq::make_index_range<D,sizeof...(I)>;
  using inv  = seq::inverse_t<std::index_sequence<I...>>;

  using reordered_specs = std::tuple<
    typename std::tuple_element_t<I,specs>...>;

  const auto& axes = hist.axes();
  const auto reordered_axes = std::tie(as_const(std::get<I>(axes))...);

  const auto reordered_nbins = make_array<ivanp::axis_size_type>(
    ( std::get<I>(axes).nbins()
      + std::tuple_element_t<I,specs>::nover::value )... );
  const auto nbins = sort(reordered_nbins, inv{});

  counter<tail::size()> tcnt(reordered_nbins, tail{});
  counter<head::size()> hcnt(reordered_nbins, head{});
  const auto n_head_bins = hcnt.size();

  const auto ii = sort(
    cat_cptr(hcnt.ii(), tcnt.ii(),
      head{}, std::make_index_sequence<tail::size()>{}),
    inv{});

  using slice_t = binner_slice<D, Bin, reordered_specs>;
  std::vector<slice_t> slices;
  slices.reserve(tcnt.size());

  const auto& hbins = hist.bins();

  for ( ; !tcnt; ++tcnt ) {
    typename slice_t::binner_type::container_type bins;
    bins.reserve(n_head_bins);
    for ( ; !hcnt; ++hcnt ) {
      bins.emplace_back(std::ref(hbins[index(ii,nbins)]));
    }
    hcnt.reset();

    slices.emplace_back( reordered_axes, tcnt.ii(), std::move(bins) );
  }

  return slices;
}

template <size_t D=1, // slicing into chunks of D dimensions
          typename Bin,
          typename... Ax,
          typename Container,
          typename Filler
>
inline auto slice( // overload for default axis order
  const binner<Bin,std::tuple<Ax...>,Container,Filler>& hist
) { return slice<D>(hist,std::index_sequence_for<Ax...>{}); }

} // end namespace ivanp

#endif
