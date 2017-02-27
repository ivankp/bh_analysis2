// Written by Ivan Pogrebnyak

#include <iostream>
#include <cstring>
#include <algorithm>
#include <memory>

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

#define test(var) \
  std::cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << std::endl;

#define _STR(S) #S
#define STR(S) _STR(S)

using std::cout;
using std::cerr;
using std::endl;

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

enum isp_t { gg, gq, qq };
isp_t get_isp(Int_t id1, Int_t id2) noexcept {
  bool g1 = (id1 == 21);
  bool g2 = (id2 == 21);
  if (g1 && g1) return gg;
  else if ((!g1) && (!g2)) return qq;
  else return gq;
}

struct hist_bin {
  static std::vector<double> weights;
  static unsigned wi;
  struct bin {
    double w = 0, w2 = 0;
    size_t n = 0;
  };
  std::vector<bin> bins;
  hist_bin(): bins(weights.size()) { }

  inline hist_bin& operator++() noexcept {
    for (unsigned i=weights.size(); i!=0; ) {
      --i;
      const double weight = weights[i];
      bin& b = bins[i];
      b.w += weight;
      b.w2 += weight*weight;
      ++b.n;
    }
    return *this;
  }
  inline hist_bin& operator+=(const hist_bin& rhs) noexcept {
    for (unsigned i=weights.size(); i!=0; ) {
      --i;
      const bin& br = rhs.bins[i];
      bin& bl = bins[i];
      bl.w += br.w;
      bl.w2 += br.w2;
      bl.n += br.n;
    }
    return *this;
  }
};
std::vector<double> hist_bin::weights;
unsigned hist_bin::wi;

namespace ivanp { namespace root {
template <> class bin_converter<hist_bin> {
  inline const auto& get(const hist_bin& bin) const noexcept {
    return bin.bins[hist_bin::wi];
  }
public:
  inline const auto& weight(const hist_bin& b) const noexcept {
    return get(b).w;
  }
  inline const auto&  sumw2(const hist_bin& b) const noexcept {
    return get(b).w2;
  }
  inline const auto&    num(const hist_bin& b) const noexcept {
    return get(b).n;
  }
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

  cout << "\033[36mBinning\033[0m: " << bins_file << '\n' << endl;
  re_axes ra(bins_file); // read axes

  // Open input ntuples root file ===================================
  TChain chain(tree_name);
  for (const char* name : ntuples) {
    if (!chain.Add(name,0)) return 1;
    cout << "\033[36mInput ntuple\033[0m: " << name << endl;
  }
  cout << endl;

  boost::optional<TChain> weights_chain;
  if (weights.size()) {
    weights_chain.emplace("weights");
    for (const char* name : weights) {
      if (!weights_chain->Add(name,0)) return 1;
      cout << "\033[36mInput weights\033[0m: " << name << endl;
    }
    cout << endl;
  }

  // Set up branches for reading
  if (weights_chain) chain.AddFriend(&*weights_chain);
  TTreeReader reader(&chain);

  TTreeReaderValue<Int_t> _id(reader,"id");
  TTreeReaderValue<Int_t> _nparticle(reader,"nparticle");
  TTreeReaderArray<Int_t> _kf(reader,"kf");

  float_or_double_array_reader _px(reader,"px");
  float_or_double_array_reader _py(reader,"py");
  float_or_double_array_reader _pz(reader,"pz");
  float_or_double_array_reader _E (reader,"E" );

  boost::optional<TTreeReaderValue<Int_t>> _ncount;
  // boost::optional<TTreeReaderValue<Char_t>> _part;
  for ( auto bo : *reader.GetTree()->GetListOfBranches() ) {
    if (!strcmp(bo->GetName(),"ncount")) _ncount.emplace(reader,"ncount");
    // else if (!strcmp(bo->GetName(),"part")) _part.emplace(reader,"part");
  }
  TTreeReaderValue<Int_t> _id1(reader,"id1");
  TTreeReaderValue<Int_t> _id2(reader,"id2");

  // handle multiple weights
  std::vector<float_or_double_value_reader> _weights;
  if (!weights_chain) _weights.emplace_back(reader,"weight");
  else {
    const TObjArray *bb = weights_chain->GetListOfBranches();
    _weights.reserve(bb->GetEntriesFast());
    for (const auto* b : *bb) {
      _weights.emplace_back(reader,static_cast<const TBranch*>(b)->GetName());
    }
  }
  hist_bin::weights.resize(_weights.size());

  // Define histograms ==============================================
  size_t ncount = 0, num_events = 0, num_selected = 0;

  hist<int> h_Njets({need_njets+2u,0,int(need_njets+2)});

#define a_(name) auto a_##name = ra[#name];
#define h_(name) re_hist<1> h_##name(#name,ra[#name]);

  a_(y) a_(phi)

  h_(HT) h_(H_pT) h_(H_y) h_(H_eta) h_(H_phi) h_(H_mass)

  auto h_jet_pT   = reserve<re_hist<1>>(need_njets+1);
  auto h_jet_y    = reserve<re_hist<1>>(need_njets+1);
  auto h_jet_eta  = reserve<re_hist<1>>(need_njets+1);
  auto h_jet_phi  = reserve<re_hist<1>>(need_njets+1);
  auto h_jet_mass = reserve<re_hist<1>>(need_njets+1);

  for (unsigned i=0; i<need_njets+1; ++i) {
    auto name = cat("jet",i+1,"_pT");
    h_jet_pT.emplace_back(name,ra[name]);
  }
  for (unsigned i=0; i<need_njets+1; ++i) {
    auto name = cat("jet",i+1,"_y");
    h_jet_y.emplace_back(name,a_y);
  }
  for (unsigned i=0; i<need_njets+1; ++i) {
    auto name = cat("jet",i+1,"_eta");
    h_jet_eta.emplace_back(name,a_y);
  }
  for (unsigned i=0; i<need_njets+1; ++i) {
    auto name = cat("jet",i+1,"_phi");
    h_jet_phi.emplace_back(name,a_phi);
  }
  for (unsigned i=0; i<need_njets+1; ++i) {
    auto name = cat("jet",i+1,"_mass");
    h_jet_mass.emplace_back(name,ra[name]);
  }

  h_(jjpT_dpT )  h_(jjfb_dpT )
  h_(jjpT_dy  )  h_(jjfb_dy  )
  h_(jjpT_deta)  h_(jjfb_deta)
  h_(jjpT_dphi)  h_(jjfb_dphi)
  h_(jjpT_mass)  h_(jjfb_mass)

  // histograms from HGamma analysis
  h_(hgam_pT_yy)
  h_(hgam_yAbs_yy)
  // h_(hgam_cosTS_yy)
  // h_(hgam_pTt_yy)

  h_(hgam_HT)

  h_(hgam_pT_j1) h_(hgam_pT_j2)
  h_(hgam_yAbs_j1) h_(hgam_yAbs_j2)

  // h_(hgam_sumTau_yyj)
  // h_(hgam_maxTau_yyj)

  h_(hgam_Dphi_j_j) //h_(hgam_Dphi_j_j_signed)
  h_(hgam_Dy_j_j)
  h_(hgam_m_jj)

  h_(hgam_pT_yyjj)
  h_(hgam_Dphi_yy_jj)

  // histograms for mtop study
  a_(_x) a_(_HT) a_(_maxdy) a_(_pT) a_(_x2)

  re_hist<1,0>
    h_xH_HT(a__x, a__HT),
    h_x1_HT(a__x, a__HT),
    h_x2_HT(a__x, a__HT),

    h_gg_xH_HT(a__x, a__HT),
    h_gg_x1_HT(a__x, a__HT),
    h_gg_x2_HT(a__x, a__HT),

    h_gq_xH_HT(a__x, a__HT),
    h_gq_x1_HT(a__x, a__HT),
    h_gq_x2_HT(a__x, a__HT),

    h_qq_xH_HT(a__x, a__HT),
    h_qq_x1_HT(a__x, a__HT),
    h_qq_x2_HT(a__x, a__HT);

  re_hist<1,0,0>
    h_xH_HT_maxdy(a__x, a__HT, a__maxdy),
    h_x1_HT_maxdy(a__x, a__HT, a__maxdy),
    h_x2_HT_maxdy(a__x, a__HT, a__maxdy),

    h_gg_xH_HT_maxdy(a__x, a__HT, a__maxdy),
    h_gg_x1_HT_maxdy(a__x, a__HT, a__maxdy),
    h_gg_x2_HT_maxdy(a__x, a__HT, a__maxdy),

    h_gq_xH_HT_maxdy(a__x, a__HT, a__maxdy),
    h_gq_x1_HT_maxdy(a__x, a__HT, a__maxdy),
    h_gq_x2_HT_maxdy(a__x, a__HT, a__maxdy),

    h_qq_xH_HT_maxdy(a__x, a__HT, a__maxdy),
    h_qq_x1_HT_maxdy(a__x, a__HT, a__maxdy),
    h_qq_x2_HT_maxdy(a__x, a__HT, a__maxdy);

  // p1's pT in bins of p2's x
  std::array<std::array< re_hist<1,0>, 3>,3> h_p1pT_p2x;
  for (auto& a : h_p1pT_p2x) for (auto& h : a) h = {a__pT,a__x2};

  // xi vs xj in bins of HT
  a_(_x3)
  re_hist<1,1,0>
    h_xH_x1_HT(a__x3,a__x3,a__HT),
    h_xH_x2_HT(a__x3,a__x3,a__HT),
    h_x1_x2_HT(a__x3,a__x3,a__HT);

  // maxdy vs maxdphi in bins of HT
  a_(_maxdy2) a_(_maxdphi2)
  re_hist<1,1,0>
    h_maxdy_maxdphi_HT(a__maxdy2,a__maxdphi2,a__HT);

  h_(Hjets_mass)

  // ================================================================

  std::vector<fj::PseudoJet> particles;
  Int_t prev_id = -1, curr_id;

  fastjet::ClusterSequence::print_banner(); // get it out of the way
  cout << jet_def.description() << endl;
  cout << "Expecting \033[36m" << need_njets
       << "\033[0m or more jets per event\n" << endl;

  // LOOP ===========================================================
  using tc = ivanp::timed_counter<Long64_t>;
  for (tc ent(reader.GetEntries(true)); reader.Next(); ++ent) {
    for (unsigned i=_weights.size(); i!=0; ) { // get weights
      --i; hist_bin::weights[i] = *_weights[i];
    }

    // Keep track of multi-entry events -----------------------------
    curr_id = *_id;
    if (prev_id!=curr_id) {
      prev_id = curr_id;
      ncount += ( _ncount ? **_ncount : 1);
      ++num_events;
    }
    // --------------------------------------------------------------

    particles.clear();
    particles.reserve(*_nparticle);
    boost::optional<TLorentzVector> higgs;

    // Read particles -----------------------------------------------
    for (size_t i=0, n=_kf.GetSize(); i<n; ++i) {
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
    // --------------------------------------------------------------

    ++num_selected;

    // Fill histograms ----------------------------------------------
    h_HT(HT);

    h_H_pT(H_pT);
    h_H_y(H_y);
    h_H_eta(higgs->Eta());
    h_H_phi(H_phi);
    h_H_mass(higgs->M());

    h_hgam_pT_yy(H_pT);
    h_hgam_yAbs_yy(std::abs(H_y));
    h_hgam_HT(HT);

    // jet histograms
    for (unsigned i=0, n=std::min(njets,need_njets+1); i<n; ++i) {
      h_jet_pT  [i](jets[i].pT  );
      h_jet_y   [i](jets[i].y   );
      h_jet_eta [i](jets[i].eta );
      h_jet_phi [i](jets[i].phi );
      h_jet_mass[i](jets[i].mass);
    }

    // find maximum rapidity separation in the event
    double max_dy = std::abs(jets[0].y-H_y);
    for (unsigned i=1; i<njets; ++i) {
      for (unsigned j=0; j<i; ++j) {
        larger(max_dy,std::abs(jets[i].y-jets[j].y));
      }
      larger(max_dy,std::abs(jets[i].y-H_y));
    }

    // find maximum phi separation in the event
    double max_dphi = dphi(jets[0].phi,H_phi);
    for (unsigned i=1; i<njets; ++i) {
      for (unsigned j=0; j<i; ++j) {
        larger(max_dphi,dphi(jets[i].phi,jets[j].phi));
      }
      larger(max_dphi,dphi(jets[i].phi,H_phi));
    }

    h_Hjets_mass(Hjets.M());

    if (njets<1) continue; // 111111111111111111111111111111111111111

    h_hgam_pT_j1(jets[0].pT);
    h_hgam_yAbs_j1(std::abs(jets[0].y));

    if (njets<2) continue; // 222222222222222222222222222222222222222

    // HT fraction x in bins on HT and max rapidity separation
    const double xH = H_pT/HT, x1 = jets[0].pT/HT, x2 = jets[1].pT/HT;

    const auto xH_HT_bin = h_xH_HT(xH, HT);
    const auto x1_HT_bin = h_x1_HT(x1, HT);
    const auto x2_HT_bin = h_x2_HT(x2, HT);

    const auto xH_HT_maxdy_bin = h_xH_HT_maxdy(xH, HT, max_dy);
    const auto x1_HT_maxdy_bin = h_x1_HT_maxdy(x1, HT, max_dy);
    const auto x2_HT_maxdy_bin = h_x2_HT_maxdy(x2, HT, max_dy);

    // p1's pT in bins of p2's x
    h_p1pT_p2x[0][0](H_pT, xH);
    h_p1pT_p2x[1][0](jets[0].pT, xH);
    h_p1pT_p2x[2][0](jets[1].pT, xH);
    h_p1pT_p2x[0][1](H_pT, x1);
    h_p1pT_p2x[1][1](jets[0].pT, x1);
    h_p1pT_p2x[2][1](jets[1].pT, x1);
    h_p1pT_p2x[0][2](H_pT, x2);
    h_p1pT_p2x[1][2](jets[0].pT, x2);
    h_p1pT_p2x[2][2](jets[1].pT, x2);

    h_xH_x1_HT(xH,x1,HT);
    h_xH_x2_HT(xH,x2,HT);
    h_x1_x2_HT(x1,x2,HT);

    h_maxdy_maxdphi_HT(max_dy,max_dphi,HT);

    const auto isp = get_isp(*_id1, *_id2);
    if (isp == gg) {
      h_gg_xH_HT.fill_bin_check(xH_HT_bin);
      h_gg_x1_HT.fill_bin_check(x1_HT_bin);
      h_gg_x2_HT.fill_bin_check(x2_HT_bin);

      h_gg_xH_HT_maxdy.fill_bin_check(xH_HT_maxdy_bin);
      h_gg_x1_HT_maxdy.fill_bin_check(x1_HT_maxdy_bin);
      h_gg_x2_HT_maxdy.fill_bin_check(x2_HT_maxdy_bin);
    } else if (isp == gq) {
      h_gq_xH_HT.fill_bin_check(xH_HT_bin);
      h_gq_x1_HT.fill_bin_check(x1_HT_bin);
      h_gq_x2_HT.fill_bin_check(x2_HT_bin);

      h_gq_xH_HT_maxdy.fill_bin_check(xH_HT_maxdy_bin);
      h_gq_x1_HT_maxdy.fill_bin_check(x1_HT_maxdy_bin);
      h_gq_x2_HT_maxdy.fill_bin_check(x2_HT_maxdy_bin);
    } else {
      h_qq_xH_HT.fill_bin_check(xH_HT_bin);
      h_qq_x1_HT.fill_bin_check(x1_HT_bin);
      h_qq_x2_HT.fill_bin_check(x2_HT_bin);

      h_qq_xH_HT_maxdy.fill_bin_check(xH_HT_maxdy_bin);
      h_qq_x1_HT_maxdy.fill_bin_check(x1_HT_maxdy_bin);
      h_qq_x2_HT_maxdy.fill_bin_check(x2_HT_maxdy_bin);
    }

    // Jet pair with highest pT .....................................
    const dijet jjpT(jets[0],jets[1]);

    h_jjpT_dpT (jjpT.dpT );
    h_jjpT_dy  (jjpT.dy  );
    h_jjpT_deta(jjpT.deta);
    h_jjpT_dphi(jjpT.dphi);
    h_jjpT_mass(jjpT.mass);
    // ..............................................................

    // Two most forward-backward jets ...............................
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

    h_jjfb_dpT (jjfb.dpT );
    h_jjfb_dy  (jjfb.dy  );
    h_jjfb_deta(jjfb.deta);
    h_jjfb_dphi(jjfb.dphi);
    h_jjfb_mass(jjfb.mass);
    // ..............................................................

    h_hgam_pT_j2(jets[1].pT);
    h_hgam_yAbs_j2(std::abs(jets[1].y));

    h_hgam_Dphi_j_j(jjpT.dphi);
    h_hgam_Dy_j_j(jjpT.dy);
    h_hgam_m_jj(jjpT.mass);

    const auto Hjj = *higgs + jjpT.p;

    h_hgam_pT_yyjj(Hjj.Pt());
    h_hgam_Dphi_yy_jj(dphi(H_phi,jjpT.p.Phi()));

    // --------------------------------------------------------------
  } // END EVENT LOOP
  // ****************************************************************

  cout << "Selected entries: " << num_selected << endl;
  cout << "Processed events: " << num_events << endl;
  cout << "ncount: " << ncount << '\n' << endl;

  // Open output root file for histograms
  auto fout = std::make_unique<TFile>(output_file_name,"recreate");
  if (fout->IsZombie()) return 1;

  // write root historgrams
  hist_bin::wi = 0;
  for (const auto& _w : _weights) {
    auto* dir = fout->mkdir(cat(_w.GetBranchName(),"_Jet",
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

    slice_to_root<2>(h_xH_x1_HT,"xH_x1","HT");
    slice_to_root<2>(h_xH_x2_HT,"xH_x2","HT");
    slice_to_root<2>(h_x1_x2_HT,"x1_x2","HT");

    slice_to_root<2>(h_maxdy_maxdphi_HT,"maxdy_maxdphi","HT");

    slice_to_root<1>(h_p1pT_p2x[0][0],   "H_pT","xH");
    slice_to_root<1>(h_p1pT_p2x[0][1],   "H_pT","x1");
    slice_to_root<1>(h_p1pT_p2x[0][2],   "H_pT","x2");
    slice_to_root<1>(h_p1pT_p2x[1][0],"jet1_pT","xH");
    slice_to_root<1>(h_p1pT_p2x[1][1],"jet1_pT","x1");
    slice_to_root<1>(h_p1pT_p2x[1][2],"jet1_pT","x2");
    slice_to_root<1>(h_p1pT_p2x[2][0],"jet2_pT","xH");
    slice_to_root<1>(h_p1pT_p2x[2][1],"jet2_pT","x1");
    slice_to_root<1>(h_p1pT_p2x[2][2],"jet2_pT","x2");

    slice_to_root<1>(h_xH_HT,"xH","HT");
    slice_to_root<1>(h_x1_HT,"x1","HT");
    slice_to_root<1>(h_x2_HT,"x2","HT");

    slice_to_root<1>(h_xH_HT_maxdy,"xH","xH_x1","HT");
    slice_to_root<1>(h_x1_HT_maxdy,"x1","xH_x1","HT");
    slice_to_root<1>(h_x2_HT_maxdy,"x2","xH_x1","HT");

    slice_to_root<1>(h_gg_xH_HT,"gg_xH","HT");
    slice_to_root<1>(h_gg_x1_HT,"gg_x1","HT");
    slice_to_root<1>(h_gg_x2_HT,"gg_x2","HT");

    slice_to_root<1>(h_gq_xH_HT,"gq_xH","HT");
    slice_to_root<1>(h_gq_x1_HT,"gq_x1","HT");
    slice_to_root<1>(h_gq_x2_HT,"gq_x2","HT");

    slice_to_root<1>(h_qq_xH_HT,"qq_xH","HT");
    slice_to_root<1>(h_qq_x1_HT,"qq_x1","HT");
    slice_to_root<1>(h_qq_x2_HT,"qq_x2","HT");

    slice_to_root<1>(h_gg_xH_HT_maxdy,"gg_xH","xH_x1","HT");
    slice_to_root<1>(h_gg_x1_HT_maxdy,"gg_x1","xH_x1","HT");
    slice_to_root<1>(h_gg_x2_HT_maxdy,"gg_x2","xH_x1","HT");

    slice_to_root<1>(h_gq_xH_HT_maxdy,"gq_xH","xH_x1","HT");
    slice_to_root<1>(h_gq_x1_HT_maxdy,"gq_x1","xH_x1","HT");
    slice_to_root<1>(h_gq_x2_HT_maxdy,"gq_x2","xH_x1","HT");

    slice_to_root<1>(h_qq_xH_HT_maxdy,"qq_xH","xH_x1","HT");
    slice_to_root<1>(h_qq_x1_HT_maxdy,"qq_x1","xH_x1","HT");
    slice_to_root<1>(h_qq_x2_HT_maxdy,"qq_x2","xH_x1","HT");

    ++hist_bin::wi;
  }

  fout->cd();
  TH1D *h_N = new TH1D("N","N",1,0,1);
  h_N->SetBinContent(1,ncount);
  h_N->SetEntries(num_selected);
  fout->Write();
  cout << "\n\033[32mOutput\033[0m: " << fout->GetName() << endl;

  return 0;
}

