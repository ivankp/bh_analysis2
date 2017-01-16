#include <iostream>

#include "re_axes.hh"

int main(int argc, char* argv[])
{
  re_axes ra("binning.txt");

  std::cout << ra["H_eta"].nbins() << std::endl;
  std::cout << ra["H_eta_test"].nbins() << std::endl;
  std::cout << ra["H_eta2_test"].nbins() << std::endl;
  std::cout << ra["H_phi"].nbins() << std::endl;

  return 0;
}
