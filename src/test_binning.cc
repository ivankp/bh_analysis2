#include <iostream>
#include <vector>

#include "binner.hh"
#include "re_axes.hh"

using std::cout;
using std::endl;

using re_axis = typename re_axes::axis_type;
template <size_t N>
using re_hist = ivanp::binner<int,
  ivanp::tuple_of_same_t<ivanp::axis_spec<re_axis,true,true>,N>>;
template <size_t N>
using re_hist_nof = ivanp::binner<int,
  ivanp::tuple_of_same_t<ivanp::axis_spec<re_axis,false,false,true>,N>>;

int main(int argc, char* argv[])
{
  re_axes ra("test.bins");

  re_hist<1> h1(ra["x1"]);
  re_hist<2> h2(ra["x1"],ra["x2"]);
  re_hist_nof<2> h3(ra["x1"],ra["x2"]);

  for (auto& b : h2.bins()) b = 0;

  h1(2);
  h2(2,5);
  h2(1,4,7);
  h2(3,6,2);
  h3(2,5);
  h3(1,5.5,7);

  cout << "h1:";
  for (const auto& b : h1.bins()) cout << ' ' << b;
  cout << endl;

  cout << "h2:";
  for (unsigned b1=0, n1=h2.nbins<1>(); b1<n1; ++b1) {
    static unsigned b = 0;
    if (b1) cout << "   ";
    for (unsigned b0=0, n0=h2.nbins<0>(); b0<n0; ++b0) {
      cout << ' ' << h2.bins()[b++];
    }
    cout << endl;
  }

  cout << "h3:";
  for (unsigned b1=0, n1=h3.nbins<1>(); b1<n1; ++b1) {
    static unsigned b = 0;
    if (b1) cout << "   ";
    for (unsigned b0=0, n0=h3.nbins<0>(); b0<n0; ++b0) {
      cout << ' ' << h3.bins()[b++];
    }
    cout << endl;
  }

  h3(3,5); // meant to throw on overflow

  return 0;
}
