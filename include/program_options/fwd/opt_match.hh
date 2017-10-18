#ifndef IVANP_OPT_MATCH_FWD_HH
#define IVANP_OPT_MATCH_FWD_HH

namespace ivanp { namespace po {
namespace detail {

struct opt_match_base {
  virtual bool operator()(const char* arg) const noexcept = 0;
  virtual ~opt_match_base() { }
  virtual std::string str() const noexcept = 0;
};

enum opt_type { long_opt, short_opt, context_opt };

// defined in .cc
bool opt_match_impl_chars(const char* arg, const char* m) noexcept;

}
}}

#endif
