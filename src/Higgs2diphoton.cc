#include <chrono>
#include <cmath>

#include "Higgs2diphoton.hh"

Higgs2diphoton::Higgs2diphoton()
: gen(std::chrono::system_clock::now().time_since_epoch().count()),
  phi_dist(0.,2*M_PI), cts_dist(-1.,1.)
{ }

std::pair<TLorentzVector,TLorentzVector>
Higgs2diphoton::operator()(const TLorentzVector& Higgs) {
  const double phi = phi_dist(gen);
  const double cts = cts_dist(gen);

  const double sts = std::sin(std::acos(cts));
  const double cos_phi = std::cos(phi);
  const double sin_phi = std::sin(phi);

  const double E  = Higgs.M()/2;
  const TVector3 boost = Higgs.BoostVector();

  TVector3 photon( cos_phi*sts, sin_phi*sts, cts );
  if (boost.Mag()!=0.) {
    TVector3 axis = boost;
    axis.SetMag(1.);
    photon.RotateUz(axis);
  }
  photon *= E;

  std::pair<TLorentzVector,TLorentzVector> diphoton( {photon,E}, {-photon,E} );
  diphoton.first.Boost(boost);
  diphoton.second.Boost(boost);

  return diphoton;
}
