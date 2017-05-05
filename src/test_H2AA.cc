#include <iostream>
#include <cmath>

#include "Higgs2diphoton.hh"

using std::cout;
using std::endl;

std::ostream& operator<<(std::ostream& os, const TLorentzVector& v) {
  os << '(' << v.Px() << ',' << v.Py() << ',' << v.Pz() << ',' << v.E() << ')';
  return os;
}

int main(int argc, char* argv[]) {
  
  Higgs2diphoton Hdecay;
  const auto diphoton = Hdecay({0.,1e2,10.,std::sqrt(1.+1e4+1e2)});

  cout << diphoton.first  << endl;
  cout << diphoton.second << endl;
  cout << diphoton.first+diphoton.second << endl;
  cout << diphoton.first.Angle(diphoton.second.Vect()) << endl;
}
