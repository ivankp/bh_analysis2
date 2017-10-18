#ifndef IVANP_REGEX_AXES_HH
#define IVANP_REGEX_AXES_HH

#include <string>
#include <memory>

#include "axis.hh"

class re_axes {
public:
  using axis_type = ivanp::ref_axis<double,
    std::shared_ptr<ivanp::abstract_axis<double>>>;

private:
  struct store;
  store *_store;

public:
  re_axes(const std::string& filename);
  ~re_axes();
  axis_type operator[](const std::string& name) const;
};

#endif
