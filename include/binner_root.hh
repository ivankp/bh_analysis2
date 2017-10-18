#ifndef IVANP_BINNER_ROOT_HH
#define IVANP_BINNER_ROOT_HH

#include <sstream>

#include "slice.hh"

namespace ivanp {
namespace root {

namespace detail {

template <size_t D> struct TH {
  static_assert( D==1 || D==2 || D==3, "ROOT has only TH1, TH2, TH3");
#ifndef ROOT_TH1
  static_assert( D!=1, "\033[33mTH1.h is no included\033[0m");
#endif
#ifndef ROOT_TH2
  static_assert( D!=2, "\033[33mTH2.h is no included\033[0m");
#endif
#ifndef ROOT_TH3
  static_assert( D!=3, "\033[33mTH3.h is no included\033[0m");
#endif
};
#ifdef ROOT_TH1
template <> struct TH<1> { using type = TH1D; };
#endif
#ifdef ROOT_TH2
template <> struct TH<2> { using type = TH2D; };
#endif
#ifdef ROOT_TH3
template <> struct TH<3> { using type = TH3D; };
#endif
template <size_t D> using TH_t = typename TH<D>::type;

}

template <typename T>
inline std::enable_if_t<std::is_integral<T>::value,double>
half_shift(T x) noexcept { return x - 0.5; }
template <typename T>
inline std::enable_if_t<!std::is_integral<T>::value,T>
half_shift(T x) noexcept { return x; }

#ifdef ROOT_TH1
template <typename A1>
TH1D* make_TH(const std::string& name, const std::tuple<A1>& axes) {
  const auto& a1 = std::get<0>(axes);
  if (a1.is_uniform()) {
    return new TH1D(name.c_str(),"",
      a1.nbins(),half_shift(a1.min()),half_shift(a1.max()));
  } else {
    return new TH1D(name.c_str(),"",
      a1.nbins(),vector_of_edges<double>(a1).data());
  }
}
#endif

#ifdef ROOT_TH2
template <typename A1, typename A2>
TH2D* make_TH(const std::string& name, const std::tuple<A1,A2>& axes) {
  const auto& a1 = std::get<0>(axes);
  const auto& a2 = std::get<1>(axes);
  if (a1.is_uniform()) {
    if (a2.is_uniform()) {
      return new TH2D(name.c_str(),"",
        a1.nbins(),a1.min(),a1.max(), a2.nbins(),a2.min(),a2.max());
    } else {
      return new TH2D(name.c_str(),"",
        a1.nbins(),a1.min(),a1.max(),
        a2.nbins(),vector_of_edges<double>(a2).data());
    }
  } else {
    const auto e1 = vector_of_edges<double>(a1);
    if (a2.is_uniform()) {
      return new TH2D(name.c_str(),"",
        a1.nbins(),e1.data(), a2.nbins(),a2.min(),a2.max());
    } else {
      return new TH2D(name.c_str(),"",
        a1.nbins(),e1.data(), a2.nbins(),vector_of_edges<double>(a2).data());
    }
  }
}
#endif

#ifdef ROOT_TH3
template <typename A1, typename A2, typename A3>
TH3D* make_TH(const std::string& name, const std::tuple<A1,A2,A3>& axes) {
  const auto& a1 = std::get<0>(axes);
  const auto& a2 = std::get<1>(axes);
  const auto& a3 = std::get<2>(axes);
  if (a1.is_uniform() && a2.is_uniform() && a3.is_uniform()) {
    return new TH3D(name.c_str(),"",
      a1.nbins(),a1.min(),a1.max(),
      a2.nbins(),a2.min(),a2.max(),
      a3.nbins(),a3.min(),a3.max());
  } else {
    return new TH3D(name.c_str(),"",
      a1.nbins(),vector_of_edges<double>(a1).data(),
      a2.nbins(),vector_of_edges<double>(a2).data(),
      a3.nbins(),vector_of_edges<double>(a3).data());
  }
}
#endif

template <typename Bin>
struct bin_converter {
  inline const auto& weight(const Bin& b) const noexcept { return b.w;  }
  inline const auto&  sumw2(const Bin& b) const noexcept { return b.w2; }
  inline const auto&    num(const Bin& b) const noexcept { return b.n;  }
};

namespace detail {

template <typename F, typename Bin>
class bin_converter_traits {
  template <typename, typename = void>
  struct _has_weight : std::false_type { };
  template <typename T> struct _has_weight<T,
    void_t<decltype( std::declval<T>().weight(std::declval<Bin>()) )>
  > : std::true_type { };

  template <typename, typename = void>
  struct _has_sumw2 : std::false_type { };
  template <typename T> struct _has_sumw2<T,
    void_t<decltype( std::declval<T>().sumw2(std::declval<Bin>()) )>
  > : std::true_type { };

  template <typename, typename = void>
  struct _has_num : std::false_type { };
  template <typename T> struct _has_num<T,
    void_t<decltype( std::declval<T>().num(std::declval<Bin>()) )>
  > : std::true_type { };

public:
  static constexpr bool has_weight = _has_weight<F>::value;
  static constexpr bool has_sumw2  = _has_sumw2 <F>::value;
  static constexpr bool has_num    = _has_num   <F>::value;
};

template <bool Use, bool Uf, typename Bins, typename F>
inline std::enable_if_t<!Use> set_weight(TH1* h, const Bins& bins, F get) { }
template <bool Use, bool Uf, typename Bins, typename F>
inline std::enable_if_t<Use> set_weight(TH1* h, const Bins& bins, F get) {
  Double_t *weight = dynamic_cast<TArrayD*>(h)->GetArray();
  size_t i = !Uf;
  for (const auto& bin : bins) weight[i] = get.weight(bin), ++i;
}

template <bool Use, bool Uf, typename Bins, typename F>
inline std::enable_if_t<!Use> set_sumw2(TH1* h, const Bins& bins, F get) { }
template <bool Use, bool Uf, typename Bins, typename F>
inline std::enable_if_t<Use> set_sumw2(TH1* h, const Bins& bins, F get) {
  h->Sumw2();
  Double_t *sumw2 = h->GetSumw2()->GetArray();
  size_t i = !Uf;
  for (const auto& bin : bins) sumw2[i] = get.sumw2(bin), ++i;
}

template <bool Use, typename Bins, typename F>
inline std::enable_if_t<!Use> set_num(TH1* h, const Bins& bins, F get) { }
template <bool Use, typename Bins, typename F>
inline std::enable_if_t<Use> set_num(TH1* h, const Bins& bins, F get) {
  Double_t n_total = 0;
  for (const auto& bin : bins) n_total += get.num(bin);
  h->SetEntries(n_total);
}

template <typename... Args>
class last_is_empty {
  template <typename Last, typename = void> struct impl : std::false_type { };
  template <typename Last>
  struct impl<Last, std::enable_if_t< std::is_empty<Last>::value > >
  : std::true_type { };
public:
  static constexpr bool value = impl<
    std::tuple_element_t<sizeof...(Args),std::tuple<int,Args...>>
  >::value;
};

} // end namespace detail

template <typename Bin, typename... Ax, typename Container, typename Filler,
          typename F = bin_converter<Bin> >
auto* to_root(
  const binner<Bin,std::tuple<Ax...>,Container,Filler>& hist,
  const std::string& name, F convert = F{}
) {
  auto *h = make_TH(name.c_str(),hist.axes());

  using traits = detail::bin_converter_traits<F,Bin>;
  using under  = typename std::tuple_element_t<0,std::tuple<Ax...>>::under;

  detail::set_weight<traits::has_weight,under::value>(h, hist.bins(), convert);
  detail::set_sumw2 <traits::has_sumw2, under::value>(h, hist.bins(), convert);
  detail::set_num   <traits::has_num>(h, hist.bins(), convert);

  return h;
};

template <size_t D, typename Bin, typename... Ax, typename... Args>
auto to_root(
  const binner_slice<D, Bin, std::tuple<Ax...>>& hist,
  const std::string& name, Args&&... args
) -> std::enable_if_t<
  detail::last_is_empty<Args...>::value, detail::TH_t<D>*
> {
  auto&& args_tup = std::forward_as_tuple(std::forward<Args>(args)...);
  std::stringstream ss;
  ss.precision(4);
  ss << name << hist.name(forward_subtuple(
    args_tup, std::make_index_sequence<sizeof...(Args)-1>{} ));
  return to_root(*hist,ss.str(),std::get<sizeof...(Args)-1>(args_tup));
};

template <size_t D, typename Bin, typename... Ax, typename... Args>
auto to_root(
  const binner_slice<D, Bin, std::tuple<Ax...>>& hist,
  const std::string& name, Args&&... args
) -> std::enable_if_t<
  !detail::last_is_empty<Args...>::value, detail::TH_t<D>*
> {
  std::stringstream ss;
  ss.precision(4);
  ss << name << hist.name(std::forward<Args>(args)...);
  return to_root(*hist,ss.str());
};

template <size_t D, typename Seq,
          typename Bin, typename... Ax, typename Container, typename Filler,
          typename... Args>
auto slice_to_root(
  const binner_slice<D, Bin, std::tuple<Ax...>>& hist,
  const std::string& name, Args&&... args
) {
  std::vector<detail::TH_t<D>*> hh;
  for (const auto& h : slice<D>(hist,Seq{}))
    hh.push_back( to_root(h,name,std::forward<Args>(args)...) );
  return hh;
}

template <size_t D=1,
          typename Bin, typename... Ax, typename Container, typename Filler,
          typename... Args>
auto slice_to_root(
  const binner<Bin,std::tuple<Ax...>,Container,Filler>& hist,
  const std::string& name, Args&&... args
) {
  std::vector<detail::TH_t<D>*> hh;
  for (const auto& h : slice<D>(hist))
    hh.push_back( to_root(h,name,std::forward<Args>(args)...) );
  return hh;
}

} // end namespace root
} // end namespace ivanp

#endif
