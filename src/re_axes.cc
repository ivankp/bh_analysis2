#include "re_axes.hh"

#include <fstream>
#include <tuple>
#include <vector>
#include <regex>
#include <utility>
#include <cctype>

#include "exception.hh"

struct re_axes::store: public
  std::vector<std::pair<std::regex,axis_type>> { };

re_axes::re_axes(const std::string& filename): _store(new store) {

  std::ifstream f(filename);
  if (!f) throw std::runtime_error("failed to open file: "+filename);

  char c; // char buffer
  bool e = false, // expression complete
       l = false, // hit left brace
       u = false, // uniform axis
       comment = false;
  std::string re, num_str;
  std::vector<double> nums;

  const auto push_num = [&](){
    if (num_str.size()!=0) {
      nums.push_back(std::stod(num_str));
      num_str.clear();
    }
  };

  while (f.get(c)) {
    if (comment) {
      if (c=='\n') comment = false;
      continue;
    } else if (c=='#') {
      comment = true;
      continue;
    }
    if (!e) {
      if (!isspace(c)) { e = true; re += c; }
    } else if (!l) {
      if (c!='{') {
        re += c;
      } else {
        l = true;
        while (isspace(re.back())) re.pop_back();
      }
    } else {
      if (u && c==':') throw ivanp::exception("more than 1 \':\'");
      if (c!='}') {
        if (isspace(c)) {
          push_num();
        } else if (c==':') {
          if ( ( num_str.size() && nums.size()==0) ||
               (!num_str.size() && nums.size()==1) ) {
            u = true;
            push_num();
          } else throw ivanp::exception("out of place \':\'");
        } else num_str += c;
      } else {
        if (u && nums.size()!=3) throw ivanp::exception(
          "more than 3 arguments for uniform axis");

        using axis_ptr = ivanp::abstract_axis<double>*;

        _store->emplace_back( std::piecewise_construct,
          std::make_tuple( std::move(re),
            std::regex::nosubs | std::regex::optimize | std::regex::extended ),
          std::make_tuple( u
            ? static_cast<axis_ptr>(
                new ivanp::uniform_axis<double,true>(
                nums[0], nums[1], nums[2]))
            : static_cast<axis_ptr>(
                new ivanp::container_axis<std::vector<double>,true>(
                std::move(nums)))
          )
        );

        if (u) nums.clear();

        e = false;
        l = false;        
        u = false;        
      }
    }
  } // end while c

}

re_axes::~re_axes() { delete _store; }

re_axes::axis_type re_axes::operator[](const std::string& name) const {
  for (const auto& ra : *_store) {
    if (std::regex_match( name, ra.first, std::regex_constants::match_any ))
      return ra.second;
  }
  throw ivanp::exception("No binning found for ", name);
}
