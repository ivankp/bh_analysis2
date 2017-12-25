#include <iostream>
// #include <vector>
// #include <array>
#include <complex>

#include <TFile.h>
#include <TDirectory.h>
#include <TH1.h>
#include <TF1.h>
#include <TFitResult.h>

#include "program_options.hh"
#include "tkey.hh"
#include "tc_msg.hh"
#include "string_alg.hh"
#include "math.hh"

#define TEST(VAR) \
  std::cout << tc::cyan << #VAR << tc::reset << " = " << VAR << std::endl;

using std::cout;
using std::cerr;
using std::endl;
namespace tc = termcolor;
using namespace ivanp;
using namespace ivanp::math;

// template <unsigned... N>
// std::array<double,sizeof...(N)> LegendrePIntegrals(double x) noexcept {
//   std::array<double,sizeof...(N)> pows { std::pow(x,N)... };
// }

/*
std::array<double,4> LegendrePIntegrals(double x) noexcept { // 0, 2, 4, 6
  const double x2 = x*x, x3 = x*x2, x5 = x3*x2, x7 = x5*x2;
  return {
    x,
    0.5*(x3-x),
    0.875*x5-1.25*x3+0.375*x,
    2.0625*x7-3.9375*x5+2.1875*x3-0.3125*x
  };
}
*/

double fitf(double* _x, double* c) {
  const double x2 = sq(*_x), x4 = x2*x2, x6 = x4*x2;

  const double p2 = 1.5*x2 - 0.5;
  const double p4 = 4.375*x4 - 3.75*x2 + 0.375;
  const double p6 = 14.4375*x6 - 19.6875*x4 + 6.5625*x2 - 0.3125;

  const auto phase = exp(std::polar<double>(1.,c[2]));

  return norm( c[0] + c[1]*phase*p2 + c[3]*p4 + c[4]*p6 );
}

TF1 *fcn;

void loop(TDirectory* dir) { // LOOP
  for (TKey& key : get_keys(dir)) {
    const TClass* key_class = get_class(key);

    if (key_class->InheritsFrom(TH1::Class())) { // HIST

      TH1 *h = read_key<TH1>(key);
      if (!starts_with(h->GetName(),"H_absCosTheta_Hj_mass")) continue;

      info("Fitting", h->GetName());
      fcn->SetParameter(2,0);
      auto result = h->Fit(fcn,"WLSV");
      TEST( result->MinFcnValue() );

    } else if (key_class->InheritsFrom(TDirectory::Class())) { // DIR
      loop(read_key<TDirectory>(key));
    }
  }
}

int main(int argc, char* argv[]) {
  const char *ifname;
  // unsigned nbins;
  bool fix_phase = false;

  try {
    using namespace ivanp::po;
    if (program_options()
        (ifname,'i',"input file",req(),pos())
        // (nbins,"--nbins","number of bins",req())
        (fix_phase,'p',"fix phase")
        .parse(argc,argv,true)) return 0;
  } catch (const std::exception& e) {
    cerr << e << endl;
    return 1;
  }
  // ================================================================

/*
  std::vector<std::array<double,4>> poly(nbins+1);
  poly.reserve(nbins+1);
  const double width = 1./nbins;
  for (unsigned i=nbins+1; i; ) { // integrals
    --i;
    poly[i] = LegendrePIntegrals(width*i);
  }
  for (unsigned i=nbins; i; ) { // partial sums
    --i;
    for (unsigned j=4; j; ) {
      --j;
      poly[i+1][j] -= poly[i][j];
    }
  }
*/

  fcn = new TF1("fit",fitf,0,1,5);
  fcn->SetParNames("c0","c2","phi2","c4","c6");

  if (fix_phase) fcn->FixParameter(2,0);
  else fcn->SetParLimits(2,-M_PI,M_PI);

  TFile fin(ifname);

  loop(&fin);
}
