// Written by Ivan Pogrebnyak

#include <iostream>
#include <cstring>
#include <algorithm>
#include <memory>
#include <regex>

#include <boost/optional.hpp>

#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TH1.h>
#include <TH2.h>
#include <TLorentzVector.h>

#include <fastjet/ClusterSequence.hh>

#include "math.hh"
#include "timed_counter.hh"
#include "catstr.hh"
#include "exception.hh"
#include "float_or_double_reader.hh"
#include "binner_root.hh"
#include "re_axes.hh"
#include "parse_args.hh"
#include "string_alg.hh"
#include "cartesian_product.hh"

#define test(var) \
  std::cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << std::endl;

#define _STR(S) #S
#define STR(S) _STR(S)

// #define NLO

using std::cout;
using std::cerr;
using std::endl;

using boost::optional;

using ivanp::cat;
using ivanp::reserve;

namespace fj = fastjet;
using namespace ivanp::math;

struct Jet {
  TLorentzVector p;
  double pT, y, eta, phi, mass;
  template <typename P>
  Jet(const P& _p): p(_p[0],_p[1],_p[2],_p[3]),
    pT(p.Pt()), y(p.Rapidity()), eta(p.Eta()), phi(p.Phi()), mass(p.M())
  { }
};

struct dijet {
  TLorentzVector p;
  double dpT, dy, deta, dphi, mass;
  dijet(const Jet& j1, const Jet& j2)
  : p(j1.p+j2.p),
    dpT (std::abs(j1.pT -j2.pT )),
    dy  (std::abs(j1.y  -j2.y  )),
    deta(std::abs(j1.eta-j2.eta)),
    dphi(ivanp::math::dphi(j1.phi,j2.phi)),
    mass(p.M()) {}
};

enum isp_t { any_isp=0, gg, gq, qq };
isp_t get_isp(Int_t id1, Int_t id2) noexcept {
  bool g1 = (id1 == 21);
  bool g2 = (id2 == 21);
  if (g1 && g1) return gg;
  else if ((!g1) && (!g2)) return qq;
  else return gq;
}

struct hist_bin {
  static double weight;

#ifdef NLO
  static int current_id;
  int id = 0;
  double wtmp = 0;
#endif

  double w = 0, w2 = 0;
  size_t n = 0;

  inline hist_bin& operator++() noexcept {
#ifdef NLO
    if (id == current_id) wtmp += weight;
    else {
      id = current_id;
      w2 += wtmp*wtmp;
      wtmp = weight;
    }
#else
    w2 += weight*weight;
#endif
    w += weight;
    ++n;
    return *this;
  }
  inline hist_bin& operator+=(const hist_bin& b) noexcept {
#ifdef NLO
    wtmp += b.wtmp;
#endif
    w += b.w;
    w2 += b.w2;
    n += b.n;
    return *this;
  }
};
double hist_bin::weight;
#ifdef NLO
int hist_bin::current_id;
#endif

namespace ivanp { namespace root {
template <> struct bin_converter<hist_bin> {
  inline auto weight(const hist_bin& b) const noexcept {
#ifdef NLO
    return b.w + b.wtmp;
#else
    return b.w;
#endif
  }
  inline auto sumw2(const hist_bin& b) const noexcept {
#ifdef NLO
    return b.w2 + sq(b.wtmp);
#else
    return b.w2;
#endif
  }
  inline auto num(const hist_bin& b) const noexcept { return b.n; }
};
}}

template <typename... Axes>
using hist_t = ivanp::binner<hist_bin,
  std::tuple<ivanp::axis_spec<Axes>...>>;
template <typename T>
using hist = hist_t<ivanp::uniform_axis<T>>;

using re_axis = typename re_axes::axis_type;
template <bool... OF>
using re_hist = ivanp::binner<hist_bin, std::tuple<
  ivanp::axis_spec<re_axis,OF,OF>...>>;

template <bool... OF>
struct isp_split {
  std::vector<re_hist<OF...>> gg, gq, qq;
  void emplace(const std::string& name, const re_axis& axis) {
    gg.emplace_back("gg_"+name,axis);
    gq.emplace_back("gq_"+name,axis);
    qq.emplace_back("qq_"+name,axis);
  }
  template <typename F, typename... I>
  isp_split(F f, I... i) {
    const auto n = ivanp::math::prod(i...);
    gg.reserve(n);
    gq.reserve(n);
    qq.reserve(n);
    ivanp::cartesian_product(
      [&,this](auto... args){ f(*this,args...); }, // bind first
      std::make_pair(I{},i)... );
  }
};

int main(int argc, char* argv[]) {
  // parse program arguments ========================================
  std::vector<const char*> ntuples, weights;
  const char* output_file_name = nullptr;
  const char* bins_file = STR(PREFIX) "/config/Hjets_mtop.bins";
  const char* tree_name = "t3";
  fj::JetDefinition jet_def;
  unsigned need_njets = 0;
  bool w_arg = false;

  for (int i=1; i<argc; ++i) {
    using namespace parse;
    if (!std::strcmp(argv[i],"w")) { w_arg = true; continue; }
    if (!std::strcmp(argv[i],"bh")) { w_arg = false; continue; }
    if (output_root_file_name(argv[i],output_file_name)) continue;
    if (root_file_name(argv[i], w_arg ? weights : ntuples )) continue;
    if (num_jets(argv[i],need_njets)) continue;
    if (bins_file_name(argv[i],bins_file)) continue;
    if (jetdef(argv[i],jet_def)) continue;
    if (parse::tree_name(argv[i],tree_name)) continue;

    cerr << "\033[31mUnexpected argument:\033[0m " << argv[i] << endl;
    return 1;
  }
  if (!ntuples.size()) {
    cerr << "\033[31mNo input root files provided\033[0m" << endl;
    return 1;
  }
  if (!output_file_name) {
    cerr << "\033[31mNo output root file provided\033[0m" << endl;
    return 1;
  }
  if (!need_njets) {
    cerr << "\033[31mNumber of jets is not specified\033[0m" << endl;
    return 1;
  }
  if (jet_def.jet_algorithm() == fj::undefined_jet_algorithm)
    jet_def = fj::JetDefinition(fj::antikt_algorithm,0.4);
  // ================================================================

  const double jet_pt_cut = 30.;
  const double jet_eta_cut = 4.4;
#ifdef NLO
  const unsigned njmax = need_njets + 1;
#else
  const unsigned njmax = need_njets;
#endif

  cout << "\033[36mBinning\033[0m: " << bins_file << '\n' << endl;
  re_axes ra(bins_file); // read axes

  // Open input ntuples root file ===================================
  TChain chain(tree_name);
  cout << "\033[36mInput ntuples\033[0m:" << endl;
  for (const char* name : ntuples) {
    if (!chain.Add(name,0)) return 1;
    cout << "  " << name << endl;
  }
  cout << endl;

  // Set up branches for reading
  TTreeReader reader(&chain);

  TTreeReaderValue<Int_t> _nparticle(reader,"nparticle");
  TTreeReaderArray<Int_t> _kf(reader,"kf");
  TTreeReaderValue<Double_t> _weight(reader,"weight2");
#ifdef NLO
  TTreeReaderValue<Int_t> _id(reader,"id");
#endif

  float_or_double_array_reader _px(reader,"px");
  float_or_double_array_reader _py(reader,"py");
  float_or_double_array_reader _pz(reader,"pz");
  float_or_double_array_reader _E (reader,"E" );

  optional<TTreeReaderValue<Int_t>> _ncount;
  for ( auto bo : *reader.GetTree()->GetListOfBranches() ) {
    if (!strcmp(bo->GetName(),"ncount")) _ncount.emplace(reader,"ncount");
  }
  TTreeReaderValue<Int_t> _id1(reader,"id1");
  TTreeReaderValue<Int_t> _id2(reader,"id2");

  // Define histograms ==============================================
  hist<int> h_Njets({njmax+1u,0,int(njmax+1)});

#define a_(name) auto a_##name = ra[#name];
#define h_(name) re_hist<1> h_##name(#name,ra[#name]);

  a_(y) a_(phi)

  h_(HT) h_(H_pT) h_(H_y) h_(H_eta) h_(H_phi) h_(H_mass)

  auto h_jet_pT   = reserve<re_hist<1>>(njmax);
  auto h_jet_y    = reserve<re_hist<1>>(njmax);
  auto h_jet_eta  = reserve<re_hist<1>>(njmax);
  auto h_jet_phi  = reserve<re_hist<1>>(njmax);
  auto h_jet_mass = reserve<re_hist<1>>(njmax);

  for (unsigned i=0; i<njmax; ++i) {
    auto name = cat("jet",i+1,"_pT");
    h_jet_pT.emplace_back(name,ra[name]);
  }
  for (unsigned i=0; i<njmax; ++i) {
    auto name = cat("jet",i+1,"_y");
    h_jet_y.emplace_back(name,a_y);
  }
  for (unsigned i=0; i<njmax; ++i) {
    auto name = cat("jet",i+1,"_eta");
    h_jet_eta.emplace_back(name,a_y);
  }
  for (unsigned i=0; i<njmax; ++i) {
    auto name = cat("jet",i+1,"_phi");
    h_jet_phi.emplace_back(name,a_phi);
  }
  for (unsigned i=0; i<njmax; ++i) {
    auto name = cat("jet",i+1,"_mass");
    h_jet_mass.emplace_back(name,ra[name]);
  }

  h_(Hjets_mass) h_(Hjets_mass_close)

#define jjpT_h_(name) \
  optional<re_hist<1>> h_jjpT_##name; \
  if (need_njets >= 2) { \
    const char* _name = need_njets > 2 ? "jjpT_"#name : "jj_"#name; \
    h_jjpT_##name.emplace( _name, ra[_name] ); \
  }
#define jjfb_h_(name) \
  optional<re_hist<1>> h_jjfb_##name; \
  if (need_njets > 2) h_jjfb_##name.emplace("jjfb_"#name,ra["jjfb_"#name]);

  jjpT_h_(dpT )  jjfb_h_(dpT )
  jjpT_h_(dy  )  jjfb_h_(dy  )
  jjpT_h_(deta)  jjfb_h_(deta)
  jjpT_h_(dphi)  jjfb_h_(dphi)
  jjpT_h_(mass)  jjfb_h_(mass)

  auto h_H_pT_isp = isp_split<1>(
    [&](isp_split<1>& hh){ hh.emplace("H_pT",h_H_pT.axis()); });

  auto h_H_y_isp = isp_split<1>(
    [&](isp_split<1>& hh){ hh.emplace("H_y",h_H_y.axis()); });

  auto h_jet_pT_isp = isp_split<1>(
    [&](isp_split<1>& hh, unsigned j){
      hh.emplace( cat("jet",j+1,"_pT"), h_jet_pT[j].axis() );
    }, njmax);

  auto h_jet_y_isp = isp_split<1>(
    [&](isp_split<1>& hh, unsigned j){
      hh.emplace( cat("jet",j+1,"_y"), h_jet_y[j].axis() );
    }, njmax);

#define j_h_(name) \
  optional<re_hist<1>> h_##name; \
  if (need_njets >= 1) h_##name.emplace(#name,ra[#name]);
#define jj_h_(name) \
  optional<re_hist<1>> h_##name; \
  if (need_njets >= 2) h_##name.emplace(#name,ra[#name]);

  // histograms from HGamma analysis
  h_(hgam_pT_yy)
  h_(hgam_yAbs_yy)
  // h_(hgam_cosTS_yy)
  // h_(hgam_pTt_yy)

  h_(hgam_HT)

  j_h_(hgam_pT_j1)   jj_h_(hgam_pT_j2)
  j_h_(hgam_yAbs_j1) jj_h_(hgam_yAbs_j2)

  // h_(hgam_sumTau_yyj)
  // h_(hgam_maxTau_yyj)

  jj_h_(hgam_Dphi_j_j) //h_(hgam_Dphi_j_j_signed)
  jj_h_(hgam_Dy_j_j)
  jj_h_(hgam_m_jj)

  jj_h_(hgam_pT_yyjj)
  jj_h_(hgam_Dphi_yy_jj)

  // histograms for mtop study
  a_(_x) a_(_HT) a_(_maxdy) a_(_pT) a_(_x2)

  using h_x_HT_maxdy_type = ivanp::add_axes_t< re_hist<1,0>,
    ivanp::axis_spec<re_axis,false,true> >;

  std::array<std::vector<h_x_HT_maxdy_type>,4> h_x_HT;
  std::vector<std::vector< re_hist<1,0> >> h_p1pT_p2x(njmax+1);
  std::vector<re_hist<1,1,0>> h_xx_HT;
  optional<re_hist<1,1,0>> h_maxdy_maxdphi_HT;

  if (need_njets >= 2) {

    for (int isp=0; isp<4; ++isp) { // any, gg, gq, qq
      static char name[] = "ii_xj_HT_maxdy";
      name[0] = (isp < qq ? 'g' : 'q');
      name[1] = (isp < gq ? 'g' : 'q');
      h_x_HT[isp].reserve(njmax+1); // +1 for Higgs
      for (unsigned j=0; j<=njmax; ++j) {
        name[4] = (j ? char('0'+j) : 'H');
        h_x_HT[isp].emplace_back(isp ? name : name+3, a__x, a__HT, a__maxdy);
      }
    }

    // p1's pT in bins of p2's x
    for (unsigned i=0; i<=njmax; ++i) {
      h_p1pT_p2x[i].reserve(njmax+1);
      for (unsigned j=0; j<=njmax; ++j) {
        h_p1pT_p2x[i].emplace_back(
          cat(i ? cat("jet",i) : "H","_pT_x",j ? char('0'+j) : 'H'),
          a__pT,a__x2);
      }
    }

    // xi vs xj in bins of HT
    a_(_x3)
    h_xx_HT.reserve((njmax*(njmax+1))/2);
    for (unsigned i=0; i<njmax; ++i) {
      for (unsigned j=i+1; j<=njmax; ++j) {
        static char name[] = "xi_xj_HT";
        name[1] = (i ? char('0'+i) : 'H');
        name[4] = (j ? char('0'+j) : 'H');
        h_xx_HT.emplace_back(name,a__x3,a__x3,a__HT);
      }
    }

    // maxdy vs maxdphi in bins of HT
    a_(_maxdy2) a_(_maxdphi2)
    h_maxdy_maxdphi_HT.emplace("maxdy_maxdphi_HT",a__maxdy2,a__maxdphi2,a__HT);

  }

  // ================================================================

  std::vector<fj::PseudoJet> particles;
#ifdef NLO
  Int_t prev_id = -1;
#endif

  fastjet::ClusterSequence::print_banner(); // get it out of the way
  cout << jet_def.description() << endl;
  cout << "Expecting \033[36m" << need_njets
       << "\033[0m or more jets per event\n" << endl;

  size_t ncount = 0, num_events = 0, num_selected = 0;

  // LOOP ===========================================================
  using tc = ivanp::timed_counter<Long64_t>;
  for (tc ent(reader.GetEntries(true)); reader.Next(); ++ent) {
    hist_bin::weight = *_weight; // Read weight

#ifdef NLO
    // Keep track of multi-entry events -----------------------------
    hist_bin::current_id = *_id;
    if (prev_id != hist_bin::current_id) {
      prev_id = hist_bin::current_id;
      ncount += ( _ncount ? **_ncount : 1);
      ++num_events;
    }
    // --------------------------------------------------------------
#else
    ncount += ( _ncount ? **_ncount : 1);
    ++num_events;
#endif

    const size_t np = *_nparticle;
    particles.clear();
    particles.reserve(np);
    optional<TLorentzVector> higgs;

    // Read particles -----------------------------------------------
    for (size_t i=0; i<np; ++i) {
      if (_kf[i] == 25) {
        higgs.emplace(_px[i],_py[i],_pz[i],_E[i]);
      } else {
        particles.emplace_back(_px[i],_py[i],_pz[i],_E[i]);
      }
    }
    if (!higgs) cerr << "\033[31mNo Higgs in entry " << ent <<"\033[0m"<< endl;
    // --------------------------------------------------------------

    // Cluster jets -------------------------------------------------
    auto fj_jets = fj::ClusterSequence(particles,jet_def) // cluster
      .inclusive_jets(jet_pt_cut); // apply pT cut
    // apply eta cut
    for (auto it=fj_jets.begin(); it!=fj_jets.end(); ) {
      if (std::abs(it->eta()) > jet_eta_cut) fj_jets.erase(it);
      else ++it;
    }
    // sort by pT
    std::sort( fj_jets.begin(), fj_jets.end(),
      [](const fj::PseudoJet& a, const fj::PseudoJet& b){
        return ( a.pt() > b.pt() );
      });
    // resulting number of jets
    const unsigned njets = fj_jets.size();
    // --------------------------------------------------------------

    // Fill Njets histograms before cuts
    h_Njets.fill_bin(njets+1); // njets+1 because njets==0 is bin 1
    // --------------------------------------------------------------

    // Apply event cuts ---------------------------------------------
    if (njets < need_njets) continue; // at least needed number of jets
    // --------------------------------------------------------------

    // Define variables ---------------------------------------------
    const double H_pT = higgs->Pt();

    std::vector<Jet> jets; // compute jet variables
    jets.reserve(njets);
    for (const auto& jet : fj_jets) jets.emplace_back(jet);

    double HT = H_pT;
    TLorentzVector Hjets = *higgs;
    for (const auto& j : jets) {
      HT += j.pT;
      Hjets += j.p;
    }

    double H_y = higgs->Rapidity();
    double H_phi = higgs->Phi();

    const auto isp = get_isp(*_id1, *_id2);
    // --------------------------------------------------------------

    ++num_selected;

    // Fill histograms ----------------------------------------------
    h_HT(HT);

    const auto bin_H_pT = h_H_pT(H_pT);
    const auto bin_H_y  = h_H_y(H_y);
    h_H_eta(higgs->Eta());
    h_H_phi(H_phi);
    h_H_mass(higgs->M());

    const auto Hjets_mass = Hjets.M();
    h_Hjets_mass(Hjets_mass);
    h_Hjets_mass_close(Hjets_mass);

    h_hgam_pT_yy(H_pT);
    h_hgam_yAbs_yy(std::abs(H_y));
    h_hgam_HT(HT);

    switch (isp) {
      case gg: h_H_pT_isp.gg.front().fill_bin(bin_H_pT);
               h_H_y_isp .gg.front().fill_bin(bin_H_y);
               break;
      case gq: h_H_pT_isp.gq.front().fill_bin(bin_H_pT);
               h_H_y_isp .gq.front().fill_bin(bin_H_y);
               break;
      case qq: h_H_pT_isp.qq.front().fill_bin(bin_H_pT);
               h_H_y_isp .qq.front().fill_bin(bin_H_y);
               break;
      default: ;
    }

    // jet histograms ...............................................
    for (unsigned i=0, n=std::min(njets,njmax); i<n; ++i) {
      const auto bin_pT = h_jet_pT  [i](jets[i].pT);
      const auto bin_y  = h_jet_y   [i](jets[i].y );
      h_jet_eta [i](jets[i].eta );
      h_jet_phi [i](jets[i].phi );
      h_jet_mass[i](jets[i].mass);

      switch (isp) {
        case gg: h_jet_pT_isp.gg[i].fill_bin(bin_pT);
                 h_jet_y_isp .gg[i].fill_bin(bin_y);
                 break;
        case gq: h_jet_pT_isp.gq[i].fill_bin(bin_pT);
                 h_jet_y_isp .gq[i].fill_bin(bin_y);
                 break;
        case qq: h_jet_pT_isp.qq[i].fill_bin(bin_pT);
                 h_jet_y_isp .qq[i].fill_bin(bin_y);
                 break;
        default: ;
      }
    }
    // ..............................................................

    if ( need_njets >= 2 ) {
      // find maximum rapidity separation in the event ................
      double max_dy = std::abs(jets[0].y-H_y);
      for (unsigned i=1; i<njets; ++i) {
        for (unsigned j=0; j<i; ++j) {
          larger(max_dy,std::abs(jets[i].y-jets[j].y));
        }
        larger(max_dy,std::abs(jets[i].y-H_y));
      }
      // ..............................................................

      // find maximum phi separation in the event .....................
      double max_dphi = dphi(jets[0].phi,H_phi);
      for (unsigned i=1; i<njets; ++i) {
        for (unsigned j=0; j<i; ++j) {
          larger(max_dphi,dphi(jets[i].phi,jets[j].phi));
        }
        larger(max_dphi,dphi(jets[i].phi,H_phi));
      }
      // ..............................................................

      (*h_maxdy_maxdphi_HT)(max_dy,max_dphi,HT);

      // HT fractions .................................................
      std::vector<double> frac_HT(njmax+1);
      auto x_HT_bin = h_x_HT[any_isp][0]( frac_HT[0] = H_pT/HT, HT, max_dy );
      if (x_HT_bin+1) h_x_HT[isp][0].fill_bin( x_HT_bin );
      for (unsigned i=njets; i!=0; --i) {
        x_HT_bin = h_x_HT[any_isp][i]( frac_HT[i] = jets[i-1].pT / HT, HT, max_dy );
        if (x_HT_bin+1) h_x_HT[isp][i].fill_bin( x_HT_bin );
      }
      // ..............................................................

      std::vector<double> pT(njets+1);
      pT[0] = H_pT;
      for (unsigned i=0; i<njets; ++i) pT[i+1] = jets[i].pT;

      for (unsigned i=0; i<=njets; ++i)
        for (unsigned j=0; j<=njets; ++j)
          h_p1pT_p2x[i][j]( pT[i], frac_HT[j] );

      for (unsigned i=0, k=0; i<njets; ++i)
        for (unsigned j=i+1; j<=njets; ++j)
          h_xx_HT[k++]( frac_HT[i], frac_HT[j], HT );
    }

    (*h_hgam_pT_j1)(jets[0].pT);
    (*h_hgam_yAbs_j1)(std::abs(jets[0].y));

    if (njets<2) continue; // 222222222222222222222222222222222222222

    // Jet pair with highest pT .....................................
    const dijet jjpT(jets[0],jets[1]);

    (*h_jjpT_dpT )(jjpT.dpT );
    (*h_jjpT_dy  )(jjpT.dy  );
    (*h_jjpT_deta)(jjpT.deta);
    (*h_jjpT_dphi)(jjpT.dphi);
    (*h_jjpT_mass)(jjpT.mass);
    // ..............................................................

    // HGam .........................................................
    (*h_hgam_pT_j2)(jets[1].pT);
    (*h_hgam_yAbs_j2)(std::abs(jets[1].y));

    (*h_hgam_Dphi_j_j)(jjpT.dphi);
    (*h_hgam_Dy_j_j)(jjpT.dy);
    (*h_hgam_m_jj)(jjpT.mass);

    const auto Hjj = *higgs + jjpT.p;

    (*h_hgam_pT_yyjj)(Hjj.Pt());
    (*h_hgam_Dphi_yy_jj)( M_PI - dphi(H_phi,jjpT.p.Phi()) );
    // --------------------------------------------------------------

    if (njets<3) continue; // 333333333333333333333333333333333333333

    // Two most forward-backward jets ...............................
    // These are the same as the two harders jets for njets == 2
    // so these histograms are duplicates for njets < 3
    unsigned jb = 0, jf = 0;
    if (njets==2) { jf = 1;
    } else {
      for (unsigned j=1; j<njets; ++j) {
        if (jets[j].y < jets[jb].y) jb = j;
        if (jets[j].y > jets[jf].y) jf = j;
      }
    }
    const dijet& jjfb = ( (jb==0 && jf==1) || (jb==1 && jf==0)
      ? jjpT : dijet(jets[jb],jets[jf]) // possibly reuse jjpT
    );

    (*h_jjfb_dpT )(jjfb.dpT );
    (*h_jjfb_dy  )(jjfb.dy  );
    (*h_jjfb_deta)(jjfb.deta);
    (*h_jjfb_dphi)(jjfb.dphi);
    (*h_jjfb_mass)(jjfb.mass);
    // ..............................................................

  } // END EVENT LOOP
  // ================================================================

  cout << "Selected entries: " << num_selected << endl;
  cout << "Processed events: " << num_events << endl;
  cout << "ncount: " << ncount << '\n' << endl;

  // Open output root file for histograms
  auto fout = std::make_unique<TFile>(output_file_name,"recreate");
  if (fout->IsZombie()) return 1;

  // write root historgrams
  auto* dir = fout->mkdir(cat(_weight.GetBranchName(),"_Jet",
      jet_def.jet_algorithm() == fj::antikt_algorithm ? "AntiKt"
    : jet_def.jet_algorithm() == fj::kt_algorithm ? "Kt"
    : jet_def.jet_algorithm() == fj::cambridge_algorithm ? "CA"
    : "", std::setprecision(2), jet_def.R()*10.
  ).c_str());
  cout << dir->GetName() << endl;
  dir->cd();

  using ivanp::root::to_root;
  using ivanp::root::slice_to_root;

  auto* h_Njets_excl = to_root(h_Njets,"jets_N_excl");
  h_Njets.integrate_left();
  auto* h_Njets_incl = to_root(h_Njets,"jets_N_incl");
  h_Njets_incl->SetEntries( h_Njets_excl->GetEntries() );

  for (auto& h : re_hist<1>::all) to_root(*h,h.name);

  for (auto& h : re_hist<1,0>::all) {
    const auto vars = ivanp::rsplit<1>(h.name,'_');
    slice_to_root(*h,vars[0],vars[1]);
  }

  std::regex _maxdy_right("(_maxdy)\\[(.+),(.+)\\)");
  std::string _maxdy_inf("_maxdy∞");
  for (auto& h : h_x_HT_maxdy_type::all) {
    h->integrate_right<2>();
    const auto vars = ivanp::rsplit<2>(h.name,'_');
    for ( auto* h : slice_to_root(*h,vars[0],vars[1],vars[2]) ) {
      std::string name = std::regex_replace(h->GetName(),_maxdy_right,"$1$3");
      // const auto inf = name.rfind(_maxdy_inf);
      // if (inf!=std::string::npos) name.erase(inf,_maxdy_inf.size());
      const auto len = _maxdy_inf.size();
      const auto pos = name.size()-len;
      if (!name.compare(pos,len,_maxdy_inf)) name.erase(pos);
      h->SetName( name.c_str() );
    }
  }

  for (auto& h : re_hist<1,1,0>::all) {
    const auto vars = ivanp::rsplit<1>(h.name,'_');
    slice_to_root<2>(*h,vars[0],vars[1]);
  }

  fout->cd();
  TH1D *h_N = new TH1D("N","N",1,0,1);
  h_N->SetBinContent(1,ncount);
  h_N->SetEntries(num_selected);
  fout->Write();
  cout << "\n\033[32mOutput\033[0m: " << fout->GetName() << endl;

  return 0;
}

