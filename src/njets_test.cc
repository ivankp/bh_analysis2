#include <iostream>
#include <iomanip>
#include <algorithm>
#include <memory>
#include <array>
#include <cmath>

#include <TFile.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TH1.h>

#include <fastjet/ClusterSequence.hh>

#include "math.hh"
#include "timed_counter.hh"
#include "float_or_double_reader.hh"
#include "binner_root.hh"
#include "parse_args/jetdef.hh"

#define test(var) \
  std::cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << std::endl;

using std::cout;
using std::cerr;
using std::endl;
using std::setw;
using ivanp::math::sq;

namespace fj = fastjet;

struct hist_bin {
  static double weight;
  static int current_id;

  int id = 0;
  double wtmp = 0, w = 0, w2 = 0;
  size_t n = 0;

  inline hist_bin& operator++() noexcept {
    if (id == current_id) wtmp += weight;
    else {
      id = current_id;
      w2 += wtmp*wtmp;
      wtmp = weight;
    }
    w += weight;
    ++n;

    return *this;
  }
};
double hist_bin::weight;
int hist_bin::current_id;

struct _w {
  inline auto weight(const hist_bin& b) const noexcept { return b.w;  }
  inline auto  sumw2(const hist_bin& b) const noexcept { return b.w2 + sq(b.wtmp); }
  inline auto    num(const hist_bin& b) const noexcept { return b.n;  }
};
struct _n {
  inline auto weight(const hist_bin& b) const noexcept { return b.n;  }
  inline auto    num(const hist_bin& b) const noexcept { return b.n;  }
};

constexpr size_t nbins = 8;
using hist = ivanp::binner<hist_bin, std::tuple<
  ivanp::axis_spec< ivanp::index_axis<>, false, false, true>>,
  std::array<hist_bin,nbins> >;

int main(int argc, char* argv[]) {
  if (argc != 4) {
    cout << "usage: " << argv[0] << " jetdef in.root out.root" << endl;
    return 0;
  }

  fj::JetDefinition jet_def;
  if (!parse::jetdef(argv[1],jet_def)) {
    cerr << "\033[31mCannot parse jet definition\033[0m: " << argv[1] << endl;
    return 1;
  }
  const auto alg_str = ivanp::cat(
      jet_def.jet_algorithm() == fj::antikt_algorithm ? "AntiKt"
    : jet_def.jet_algorithm() == fj::kt_algorithm ? "Kt"
    : jet_def.jet_algorithm() == fj::cambridge_algorithm ? "CA"
    : "", std::setprecision(2), jet_def.R()*10.
  );

  auto fin = std::make_unique<TFile>(argv[2],"read");
  if (fin->IsZombie()) return 1;

  TTreeReader reader("t3",fin.get());

  TTreeReaderValue<Int_t> _nparticle(reader,"nparticle");
  TTreeReaderValue<Int_t> _id(reader,"id");
  TTreeReaderArray<Int_t> _kf(reader,"kf");

  float_or_double_array_reader _px(reader,"px");
  float_or_double_array_reader _py(reader,"py");
  float_or_double_array_reader _pz(reader,"pz");
  float_or_double_array_reader _E (reader,"E" );
  
  float_or_double_value_reader _weight(reader, "weight");

#define h_(name) hist h_particles_##name("particles_" #name,{0,nbins});
  h_(nocut)
  h_(pt30)
  h_(pt30_y44)
  h_(y44)
#undef h_
#define h_(name) hist h_jets_##name("jets_"+alg_str+"_"+#name,{0,nbins});
  h_(nocut)
  h_(pt30)
  h_(pt30_y44)
  h_(y44)
#undef h_

  fastjet::ClusterSequence::print_banner(); // get it out of the way
  cout << jet_def.description() << endl;
  std::vector<fj::PseudoJet> particles;
  particles.reserve(nbins);

  // LOOP ***********************************************************
  using tc = ivanp::timed_counter<Long64_t>;
  for (tc ent(reader.GetEntries(true)); reader.Next(); ++ent) {
    hist_bin::weight = *_weight; // Read weight
    hist_bin::current_id = *_id;

    // Read particles -----------------------------------------------
    particles.clear();
    size_t n_pt30 = 0, n_pt30_y44 = 0, n_y44 = 0;
    for (size_t i=0, n=_kf.GetSize(); i<n; ++i) {
      if (_kf[i] != 25) {
        particles.emplace_back(_px[i],_py[i],_pz[i],_E[i]);
        if (particles.back().pt() > 30.) {
          ++n_pt30;
          if (particles.back().rap() < 4.4) ++n_pt30_y44;
        }
        if (particles.back().rap() < 4.4) ++n_y44;
      }
    }
    // --------------------------------------------------------------

    h_particles_nocut(particles.size());
    h_particles_pt30(n_pt30);
    h_particles_pt30_y44(n_pt30_y44);
    h_particles_y44(n_y44);

    // AntiKt4 ------------------------------------------------------
    const auto fj_jets = fj::ClusterSequence(particles,jet_def) // cluster
      .inclusive_jets();

    n_pt30 = 0; n_pt30_y44 = 0; n_y44 = 0;
    for (const auto& jet : fj_jets) {
      if (jet.pt() > 30.) {
        ++n_pt30;
        if (jet.rap() < 4.4) ++n_pt30_y44;
      }
      if (jet.rap() < 4.4) ++n_y44;
    }
    // --------------------------------------------------------------

    h_jets_nocut(fj_jets.size());
    h_jets_pt30(n_pt30);
    h_jets_pt30_y44(n_pt30_y44);
    h_jets_y44(n_y44);

    // --------------------------------------------------------------
  }
  // ****************************************************************

  auto fout = std::make_unique<TFile>(argv[3],"recreate");
  if (fout->IsZombie()) return 1;

  for (const auto& h : hist::all) {
    ivanp::root::to_root(*h,"w_"+h.name,_w{});
    ivanp::root::to_root(*h,"n_"+h.name,_n{});
  }

  fout->Write();
  cout << "\n\033[32mOutput\033[0m: " << fout->GetName() << endl;

  return 0;
}
