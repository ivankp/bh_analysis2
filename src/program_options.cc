#include <iostream>
#include <cstring>
#include <cctype>
#include <stdexcept>

#define IVANP_PROGRAM_OPTIONS_CC
#include "program_options.hh"

#define TEST(var) \
  std::cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << std::endl;

using std::cout;
using std::cerr;
using std::endl;
using namespace std::string_literals;

bool is_number(const char* str) noexcept {
// https://stackoverflow.com/q/4654636/2640636
#ifdef PROGRAM_OPTIONS_BOOST_LEXICAL_CAST
  static double d;
  return boost::conversion::try_lexical_convert(str,d);
#else
  char* p;
  std::strtod(str,&p);
  return !*p;
#endif
}

namespace ivanp { namespace po {

namespace detail {

bool opt_match_impl_chars(const char* arg, const char* m) noexcept {
  int i = 0;
  for (; m[i]!='\0' && arg[i]!='\0'; ++i)
    if ( arg[i]!=m[i] ) return false;
  return m[i]=='\0' && arg[i]=='\0';
}

opt_type get_opt_type(const char* arg) noexcept {
  unsigned char n = 0;
  char c;
  for (c=arg[n]; c=='-'; c=arg[++n]) ;
  if (c=='\0') n = 0;
  switch (n) {
    case  1: return is_number(arg) ? context_opt : short_opt;
    case  2: return long_opt;
    default: return context_opt;
  }
}

}

void check_count(detail::opt_def* opt) {
  if (!opt->is_multi() && opt->count)
    throw error("too many options " + opt->name);
}

bool program_options::parse(int argc, char const * const * argv,
                            bool help_if_no_args) {
  using namespace ::ivanp::po::detail;
  // --help takes precedence over everything else
  if (help_if_no_args && argc==1) {
    help();
    return true;
  }
  for (int i=1; i<argc; ++i) {
    for (const char* h : help_flags) {
      if (!strcmp(h,argv[i])) {
        help();
        return true;
      }
    }
  }

  opt_def *opt = nullptr;
  const char* val = nullptr;
  std::string tmp; // use view here
  bool last_was_val = false;

  for (int i=1; i<argc; ++i) {
    const char* arg = argv[i];
    last_was_val = false;

    const auto opt_type = get_opt_type(arg);
#ifdef PROGRAM_OPTIONS_DEBUG
    cout << arg << ' ' << opt_type << endl;
#endif

    // ==============================================================

    if (opt_type!=context_opt) {
      if (opt) {
        if (!opt->count) opt->as_switch();
        opt = nullptr;
      }
      if (opt_type==long_opt) { // long: split by '='
        if ((val = strchr(arg,'='))) arg = tmp.assign(arg,val).c_str(), ++val;
      } else { // short: allow spaceless
        if (arg[2]!='\0') val = arg+2;
      }
    }

    // ==============================================================

    if (!opt || (opt->is_multi() && opt->count)) {
      for (auto& m : matchers[opt_type]) {
        if ((*m.first)(arg)) { // match
          opt = m.second;
#ifdef PROGRAM_OPTIONS_DEBUG
          cout << arg << " matched: " << opt->name << endl;
#endif
          check_count(opt);
          if (opt_type==context_opt) val = arg;
          if (opt->is_switch()) {
            if (val) {
              if (opt_type!=context_opt) throw po::error(
                "switch " + opt->name + " does not take arguments");
              else val = nullptr;
            }
            opt->as_switch(), opt = nullptr;
          } else if (val) {
            opt->parse(val), val = nullptr;
            last_was_val = true;
            if (!opt->is_multi()) opt = nullptr;
          }
          goto next_arg;
        }
      }
    }

    if (opt) {
#ifdef PROGRAM_OPTIONS_DEBUG
      cout << arg << " arg of: " << opt->name << endl;
#endif
      opt->parse(arg);
      last_was_val = true;
      if (!opt->is_multi()) opt = nullptr;
      continue;
    }

    // handle positional options
    if (opt_type==context_opt && !opt && pos.size()) {
      auto *pos_opt = pos.front();
      check_count(pos_opt);
#ifdef PROGRAM_OPTIONS_DEBUG
      cout << arg << " pos: " << pos_opt->name << endl;
#endif
      pos_opt->parse(arg);
      last_was_val = true;
      if (!pos_opt->is_pos_end()) pos.pop();
      continue;
    }

    throw po::error("unexpected option ",arg);
    next_arg: ;
  } // end arg loop
  if (opt) {
    if (!opt->count) opt->as_switch();
    else if (!last_was_val) throw po::error("dangling option " + opt->name);
    opt = nullptr;
  }

  for (opt_def *opt : req) // check required passed
    if (!opt->count) throw error("missing required option "+opt->name);

  for (opt_def *opt : default_init) // init with default values
    if (!opt->count) opt->default_init();

  return false;
}

inline bool streq_ignorecase(const char* str, const char* s1) {
  const char *a=str, *b=s1;
  for (; *a!='\0' && *b!='\0'; ++a, ++b)
    if (toupper(*a) != toupper(*b)) return false;
  return *a=='\0' && *b=='\0';
}
template <typename S1>
inline bool streq_any_ignorecase(const char* str, S1 s1) {
  return streq_ignorecase(str,s1);
}
template <typename S1, typename... Ss>
inline bool streq_any_ignorecase(const char* str, S1 s1, Ss... ss) {
  return streq_ignorecase(str,s1) || streq_any_ignorecase(str,ss...);
}

namespace detail {

void arg_parser_impl_bool(const char* arg, bool& var) {
  if (streq_any_ignorecase(arg,"1","TRUE","YES","ON","Y")) var = true;
  else if (streq_any_ignorecase(arg,"0","FALSE","NO","OFF")) var = false;
  else throw po::error('\"',arg,"\" cannot be interpreted as bool");
}

}

unsigned utf_len(const char* s) {
  // https://stackoverflow.com/a/4063229/2640636
  unsigned len = 0;
  while (*s) len += (*s++ & 0xc0) != 0x80;
  return len;
}

std::vector<unsigned> wrap(std::string& str, unsigned w) {
  std::vector<unsigned> br;
  unsigned d = 0, l = 0;
  for (unsigned i=0, n=str.size(); i<n; ++i, ++l) {
    char c = str[i];
    if (c==' '||c=='\t') {
      if (l>w) {
        if (l==i-d) c = '\n', l = 0, br.push_back(i);
        else str[d] = '\n', l = i-d, br.push_back(d);
      }
      d = i;
    } else if (c=='\n') {
      char& cd = str[d];
      if (l>w && (cd==' '||cd=='\t')) cd = '\n', br.push_back(d);
      d = i, l = 0;
      br.push_back(d);
    }
  }
  if (l>w) str[d] = '\n', br.push_back(d);
  return br;
}

const std::string& fmt_descr(std::string& str, unsigned t, unsigned w) {
  const auto br = wrap(str,w-t); // replace spaces with \n for line wrapping
  if (br.size()==0) return str;

  std::string buff;
  buff.reserve(str.size()+br.size()*t);
  buff.append(str,0,br.front()+1);
  size_t pos;
  for (size_t i=1, n=br.size(); i<n; ++i) {
    pos = br[i-1];
    buff.append(t,' ');
    buff.append(str,pos+1,br[i]-pos);
  }
  pos = br.back();
  buff.append(t,' ');
  buff.append(str,pos+1,str.size()-pos);

  str = std::move(buff);
  return str;
}

void program_options::help() {
  static constexpr unsigned line_width = 80;
  if (help_prefix_str.size()) {
    wrap(help_prefix_str,line_width);
    cout << help_prefix_str << "\n\n";
  }

  cout << "Options:\n";

  const unsigned ndefs = opt_defs.size();
  std::array<unsigned,2> w{0,0}; // name and marks widths
  std::vector<unsigned> lens(ndefs);
  std::vector<std::string> marks(ndefs);

  for (unsigned d=0; d<ndefs; ++d) {
    const auto& opt = opt_defs[d];
    unsigned len = lens[d] = utf_len(opt->name.c_str());
    if (len > w[0]) w[0] = len;
    len = 0;
    if (opt->is_req()        ) marks[d] += '*', ++len;
    if (opt->is_switch()     ) marks[d] += '-', ++len;
    if (opt->is_switch_init()) marks[d] += '?', ++len;
    if (opt->is_pos()        ) marks[d] += '^', ++len;
    if (len > w[1]) w[1] = len;
  }

  w[0] += 1;
  const unsigned tab = w[0]+w[1]+3;
  for (unsigned d=0; d<ndefs; ++d) {
    const auto& opt = opt_defs[d];
    // name
    cout << "  " << opt->name;
    for (unsigned i=0, n=w[0]-lens[d]; i<n; ++i) cout << ' ';
    // marks
    if (w[1]) {
      cout << marks[d];
      for (unsigned i=0, n=w[1]-marks[d].size()+1; i<n; ++i) cout << ' ';
    }
    // description
    cout << fmt_descr(opt->descr,tab,line_width) << '\n';
  }

  if (w[1]) {
    cout << "\nannotation:\n"
      "  * required\n"
      "  - switch\n"
      "  ? argument is optional\n"
      "  ^ positional\n";
  }

  if (help_suffix_str.size()) {
    wrap(help_suffix_str,line_width);
    cout <<'\n'<< help_suffix_str << '\n';
  }
  cout.flush();
}

}} // end namespace ivanp
