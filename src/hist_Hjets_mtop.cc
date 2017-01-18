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
#include <TLorentzVector.h>

#include <fastjet/ClusterSequence.hh>

#include "math.hh"
#include "binner.hh"
#include "re_axes.hh"
#include "timed_counter.hh"
#include "catstr.hh"
#include "exception.hh"

#define test(var) \
  std::cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << std::endl;

#include "burst.hh"

using std::cout;
using std::cerr;
using std::endl;

using ivanp::cat;

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

struct hist_bin {
  static double weight;
  double w, w2;
  size_t n;
  hist_bin(): w(0.), w2(0.), n(0) { }
  inline void operator++() noexcept {
    w += weight;
    w2 += weight*weight;
    ++n;
  }
  inline hist_bin& operator+=(const hist_bin& b) noexcept {
    w += b.w;
    w2 += b.w2;
    n += b.n;
    return *this;
  }
};
double hist_bin::weight;

template <typename... Axes>
using hist_t = ivanp::binner<hist_bin,
  std::tuple<ivanp::axis_spec<Axes>...>>;
template <typename T>
using hist = hist_t<ivanp::uniform_axis<T>>;

using re_axis = typename re_axes::axis_type;
template <size_t N>
using re_hist = ivanp::binner<hist_bin,
  ivanp::tuple_of_same_t<ivanp::axis_spec<re_axis>,N>>;

template <typename A, typename B>
TH1D* make_TH1D(const char* name, const A& axis, B begin, B end) {
  TH1D* h = new TH1D(name,"",axis.nbins(),axis.min(),axis.max());
  h->Sumw2();
  TArrayD& sumw2 = *h->GetSumw2();
  size_t n_total = 0, i = 0;
  for (auto bin=begin; bin!=end; ++bin, ++i) {
    h->SetBinContent(i,bin->w);
    sumw2[i] = bin->w2;
    n_total += bin->n;
  }
  h->SetEntries(n_total);
  return h;
}
template <size_t D, typename... T>
std::enable_if_t<D==1,TH1D*>
make_TH(const std::string& name, const ivanp::binner_slice<T...>& s) {
  TH1D* h = new TH1D(name.c_str(),"", s.end-s.begin-2,
    std::get<0>(s.ranges)[0], std::get<0>(s.ranges)[1]);
  h->Sumw2();
  TArrayD& sumw2 = *h->GetSumw2();
  size_t n_total = 0, i = 0;
  for (auto bin=s.begin; bin!=s.end; ++bin, ++i) {
    h->SetBinContent(i,bin->w);
    sumw2[i] = bin->w2;
    n_total += bin->n;
  }
  h->SetEntries(n_total);
  return h;
}

template <size_t D, typename... A>
auto burst(
  const ivanp::binner<hist_bin,std::tuple<A...>>& h
) {
  return ivanp::burst<D>( h.axes(), h.bins().begin(), h.bins().begin() );
}

template <size_t D, typename... A>
void make_root_hists(
  const ivanp::binner<hist_bin,std::tuple<A...>>& h,
  std::initializer_list<std::string> names
) {
  if (names.size() == sizeof...(A)-D) throw ivanp::exception(
    *names.begin(),": number of names doesn't match dimensions");
  const auto slices = ivanp::burst<D>(
    h.axes(), h.bins().begin(), h.bins().begin() );
  for (const auto& s : slices) {
    auto it = names.begin();
    std::string name = *(it++);
    for ( ; it!=names.end(); ++it) {
      name += cat('_',*it,
        "_[",std::get<0>(s.slices)[0],',', std::get<0>(s.slices)[1],')');
    }
    make_TH<1>(name,s);
  }
}

template <typename A>
TH1D* root_hist(const ivanp::binner<hist_bin,std::tuple<A>>& h,
  const std::string& name
) {
  return make_TH1D(name.c_str(),h.axis(),h.bins().begin(),h.bins().end());
}

// template <typename A1, typename A2>
// void root_hist(const ivanp::binner<hist_bin,std::tuple<A1,A2>>& h,
//   const std::string& name,
void root_hist(const re_hist<2>& h, const std::string& name,
  const std::string& var2
) {
  const auto& a1 = h.axis<0>();
  const auto& a2 = h.axis<1>();
  const auto nbins1 = h.nbins<0>();
  const auto nbins2 = h.nbins<1>();
  auto it = h.bins().begin() + nbins1;
  for (unsigned i=1; i<nbins2-1; ++i) {
    make_TH1D(cat(name,'_',var2,"_[",a2.lower(i),',',a2.upper(i),')').c_str(),
      a1,it,it+nbins1);
    it += nbins1;
  }
}
void root_hist(const re_hist<3>& h, const std::string& name,
  const std::string& var2, const std::string& var3
) {
  const auto& a1 = h.axis<0>();
  const auto& a2 = h.axis<1>();
  const auto& a3 = h.axis<2>();
  const auto nbins1 = h.nbins<0>();
  const auto nbins2 = h.nbins<1>();
  const auto nbins3 = h.nbins<2>();
  auto it = h.bins().begin() + nbins1*nbins2;
  for (unsigned j=1; j<nbins3-1; ++j) {
    it += nbins1;
    for (unsigned i=1; i<nbins2-1; ++i) {
      make_TH1D(cat( name,
          '_',var2,"_[",a2.lower(i),',',a2.upper(i),')',
          '_',var3,"_[",a3.lower(j),',',a3.upper(j),')'
        ).c_str(), a1,it,it+nbins1);
      it += nbins1;
    }
    it += nbins1;
  }
}

int main(int argc, char* argv[])
{
  fj::JetDefinition jet_def(fj::antikt_algorithm,0.4);
  const double jet_pt_cut = 30.;
  const double jet_eta_cut = 4.4;
  constexpr unsigned need_njets = 2;

  // Define histograms ==============================================
  size_t ncount = 0, num_events = 0, num_selected = 0;

  hist<int> h_Njets({need_njets+2,0,need_njets+2});

  re_axes ra("Hjets_mtop.bins");
#define a_(name) auto a_##name = ra[#name];
#define h_(name) re_hist<1> h_##name(#name,ra[#name]);

  a_(y) a_(phi)

  h_(HT) h_(H_pT) h_(H_y) h_(H_eta) h_(H_phi) h_(H_mass)

  std::vector<re_hist<1>> h_jet_pT, h_jet_y, h_jet_eta, h_jet_phi, h_jet_mass;
    h_jet_pT.reserve(need_njets+1);
    h_jet_y.reserve(need_njets+1);
    h_jet_eta.reserve(need_njets+1);
    h_jet_phi.reserve(need_njets+1);
    h_jet_mass.reserve(need_njets+1);
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

  // histograms for mtop study
  a_(_x) a_(_HT) a_(_maxdy) a_(_pT) a_(_x2)

  re_hist<3>
    h_xH_HT_maxdy(a__x, a__HT, a__maxdy),
    h_x1_HT_maxdy(a__x, a__HT, a__maxdy),
    h_x2_HT_maxdy(a__x, a__HT, a__maxdy);

  // p1's pT in bins of p2's x
  std::array<std::array< re_hist<2>, 3>,3> h_p1pT_p2x;
  for (auto& a : h_p1pT_p2x) for (auto& h : a) h = {a__pT,a__x2};

  // ================================================================

  // Open input ntuple root file
  auto fin = std::make_unique<TFile>(argv[1],"read");
  if (fin->IsZombie()) return 1;

  // Set up branches for reading
  TTreeReader reader("t3",fin.get());
  if (!reader.GetTree()) {
    cerr << "\033[31mNo tree \"t3\" in file "
         << argv[1] <<"\033[0m"<< endl;
    return 1;
  }

  using p4_t = Float_t;
  TTreeReaderValue<Int_t> _id(reader,"id");
  TTreeReaderValue<Int_t> _nparticle(reader,"nparticle");
  TTreeReaderArray<Int_t> _kf(reader,"kf");
  TTreeReaderArray<p4_t>  _px(reader,"px");
  TTreeReaderArray<p4_t>  _py(reader,"py");
  TTreeReaderArray<p4_t>  _pz(reader,"pz");
  TTreeReaderArray<p4_t>  _E (reader,"E");
  TTreeReaderValue<Double_t> _weight(reader,"weight");
  boost::optional<TTreeReaderValue<Int_t>> _ncount;
  // boost::optional<TTreeReaderValue<Char_t>> _part;
  for ( auto bo : *reader.GetTree()->GetListOfBranches() ) {
    if (!strcmp(bo->GetName(),"ncount")) _ncount.emplace(reader,"ncount");
    // else if (!strcmp(bo->GetName(),"part")) _part.emplace(reader,"part");
  }

  std::vector<fj::PseudoJet> particles;
  Int_t prev_id = -1, curr_id;

  // LOOP ***********************************************************
  using tc = ivanp::timed_counter<Long64_t>;
  for (tc ent(reader.GetEntries(true)); reader.Next() /* && ent < 1 */; ++ent) {
    hist_bin::weight = *_weight; // Read weight

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
      if (abs(it->eta()) > jet_eta_cut) fj_jets.erase(it);
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
    for (const auto& j : jets) HT += j.pT;

    double H_y = higgs->Rapidity();
    // --------------------------------------------------------------

    ++num_selected;

    // Fill histograms ----------------------------------------------
    h_HT(HT);

    h_H_pT(H_pT);
    h_H_y(H_y);
    h_H_eta(higgs->Eta());
    h_H_phi(higgs->Phi());
    h_H_mass(higgs->M());

    // jet histograms
    for (unsigned i=0, n=std::min(njets,need_njets+1); i<n; ++i) {
      h_jet_pT  [i](jets[i].pT  );
      h_jet_y   [i](jets[i].y   );
      h_jet_eta [i](jets[i].eta );
      h_jet_phi [i](jets[i].phi );
      h_jet_mass[i](jets[i].mass);
    }

    // find maximum rapidity separation in the event
    double max_dy = abs(jets[0].y-H_y);
    for (unsigned i=1; i<njets; ++i) {
      for (unsigned j=0; j<i; ++j) {
        larger(max_dy,abs(jets[i].y-jets[j].y));
      }
      larger(max_dy,abs(jets[i].y-H_y));
    }

    if (njets<2) continue; // 222222222222222222222222222222222222222

    // HT fraction x in bins on HT and max rapidity separation
    const double xH = H_pT/HT, x1 = jets[0].pT/HT, x2 = jets[1].pT/HT;
    h_xH_HT_maxdy(xH, HT, max_dy);
    h_x1_HT_maxdy(x1, HT, max_dy);
    h_x2_HT_maxdy(x2, HT, max_dy);

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

    // --------------------------------------------------------------
  } // END EVENT LOOP
  // ****************************************************************

  cout << "Selected entries: " << num_selected << endl;
  cout << "Processed events: " << num_events << endl;
  cout << "ncount: " << ncount << endl;

  // Open input ntuple root file
  auto fout = std::make_unique<TFile>(argv[2],"recreate");
  if (fout->IsZombie()) return 1;

  fout->mkdir("weight2_JetAntiKt4")->cd();

  // write root historgrams
  root_hist(h_Njets,"jets_N_excl");
  h_Njets.integrate_left();
  root_hist(h_Njets,"jets_N_incl");

  for (auto& h : re_hist<1>::all) root_hist(*h,h.name);

  root_hist(h_xH_HT_maxdy,"xH","HT","maxdy");
  root_hist(h_x1_HT_maxdy,"x1","HT","maxdy");
  root_hist(h_x2_HT_maxdy,"x2","HT","maxdy");

  // make_root_hists<1>(h_p1pT_p2x[0][0],{"test_H_pT","xH"});

  root_hist(h_p1pT_p2x[0][0],   "H_pT","xH");
  root_hist(h_p1pT_p2x[0][1],   "H_pT","x1");
  root_hist(h_p1pT_p2x[0][2],   "H_pT","x2");
  root_hist(h_p1pT_p2x[1][0],"jet1_pT","xH");
  root_hist(h_p1pT_p2x[1][1],"jet1_pT","x1");
  root_hist(h_p1pT_p2x[1][2],"jet1_pT","x2");
  root_hist(h_p1pT_p2x[2][0],"jet2_pT","xH");
  root_hist(h_p1pT_p2x[2][1],"jet2_pT","x1");
  root_hist(h_p1pT_p2x[2][2],"jet2_pT","x2");

  make_root_hists<1>(h_p1pT_p2x[0][0],{   "test_H_pT","xH"});
  make_root_hists<1>(h_p1pT_p2x[0][1],{   "test_H_pT","x1"});
  make_root_hists<1>(h_p1pT_p2x[0][2],{   "test_H_pT","x2"});
  make_root_hists<1>(h_p1pT_p2x[1][0],{"test_jet1_pT","xH"});
  make_root_hists<1>(h_p1pT_p2x[1][1],{"test_jet1_pT","x1"});
  make_root_hists<1>(h_p1pT_p2x[1][2],{"test_jet1_pT","x2"});
  make_root_hists<1>(h_p1pT_p2x[2][0],{"test_jet2_pT","xH"});
  make_root_hists<1>(h_p1pT_p2x[2][1],{"test_jet2_pT","x1"});
  make_root_hists<1>(h_p1pT_p2x[2][2],{"test_jet2_pT","x2"});

  // test( burst<1>(h_p1pT_p2x[0][0]) );
  // test( burst<1>(h_p1pT_p2x[0][1]) );
  // test( burst<1>(h_p1pT_p2x[0][2]) );
  // test( burst<1>(h_p1pT_p2x[1][0]) );
  // test( burst<1>(h_p1pT_p2x[1][1]) );
  // test( burst<1>(h_p1pT_p2x[1][2]) );
  // test( burst<1>(h_p1pT_p2x[2][0]) );
  // test( burst<1>(h_p1pT_p2x[2][1]) );
  // test( burst<1>(h_p1pT_p2x[2][2]) );

  fout->cd();
  (new TH1D("N","N",1,0,1))->SetBinContent(1,ncount);
  fout->Write();

  return 0;
}

