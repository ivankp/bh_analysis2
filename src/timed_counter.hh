// Written by Ivan Pogrebnyak

#ifndef IVANP_TIMED_COUNTER_CHRONO_HH
#define IVANP_TIMED_COUNTER_CHRONO_HH

#include <iostream>
#include <iomanip>
#include <sstream>
#include <chrono>

namespace ivanp {

class comma_numpunct: public std::numpunct<char> {
protected:
  virtual char do_thousands_sep() const { return '\''; }
  virtual std::string do_grouping() const { return "\03"; }
};

template <typename CntType, typename Compare = std::less<CntType>>
class timed_counter {
  static_assert( std::is_integral<CntType>::value,
    "Cannot instantiate timed_counter of non-integral type");
public:
  using value_type    = CntType;
  using compare_type  = Compare;
  using clock_type    = std::chrono::system_clock;
  using time_type     = std::chrono::time_point<clock_type>;
  using duration_type = std::chrono::duration<double>;

private:
  value_type cnt, cnt_start, cnt_end;
  time_type start, last;
  bool past1s;
  int nb;
  compare_type cmp;

  static const std::locale cnt_locale;

public:
  void print() {
    using std::cout;
    using std::setw;
    using std::setfill;
    std::ios::fmtflags f(cout.flags());
    auto prec = cout.precision();

    const auto dt = duration_type((last = clock_type::now()) - start).count();
    const int hours   = dt/3600;
    const int minutes = (dt-hours*3600)/60;
    const int seconds = (dt-hours*3600-minutes*60);

    for (int i=0; i<nb; ++i) cout << '\b';
    std::stringstream cnt_ss;
    cnt_ss.imbue(cnt_locale);
    cnt_ss << setw(14) << cnt;
    cout << cnt_ss.rdbuf() << " | ";
    cout.precision(2);
    cout << std::fixed << setw(6) << (
      cnt==cnt_start ? 0. : 100.*float(cnt-cnt_start)/float(cnt_end-cnt_start)
    ) <<'%'<< " | ";
    if (hours) {
      if (nb<38) nb = 38;
      cout << setw(5) << hours << ':'
      << setfill('0') << setw(2) << minutes << ':'
      << setw(2) << seconds << setfill(' ');
    } else if (minutes) {
      if (nb<32) nb = 32;
      cout << setw(2) << minutes << ':'
      << setfill('0') << setw(2) << seconds << setfill(' ');
    } else {
      cout << setw(2) << seconds <<'s';
    }

    cout.flush();
    cout.flags(f);
    cout.precision(prec);
  }
  void print_check() {
    if ( duration_type(clock_type::now()-last).count() > 1 ) {
      print();
      if (!past1s) past1s = true;
    }
  }
  void print_last() {
    print();
    if (!past1s) {
      std::chrono::duration<double,std::milli> ms(clock_type::now() - start);
      std::cout << "\b\b\b" << int(ms.count()) << "ms";
    }
    std::cout << std::endl;
  }

  timed_counter() { }
  timed_counter(value_type i, value_type n)
  : cnt(i), cnt_start(i), cnt_end(n), start(clock_type::now()), last(start),
    past1s(false), nb(30)
  { print(); }
  timed_counter(value_type n)
  : cnt(0), cnt_start(0), cnt_end(n), start(clock_type::now()), last(start),
    past1s(false), nb(30)
  { print(); }
  ~timed_counter() { print_last(); }

  void set(value_type i, value_type n) {
    cnt = i;
    cnt_start = i;
    cnt_end = n;
    last = start = clock_type::now();
    past1s = false;
    nb = 30;
  }
  void set(value_type n) {
    cnt = 0;
    cnt_start = 0;
    cnt_end = n;
    last = start = clock_type::now();
    past1s = false;
    nb = 30;
  }

  inline bool ok() const noexcept { return cmp(cnt,cnt_end); }

  // prefix
  inline value_type operator++() { print_check(); return ++cnt; }
  inline value_type operator--() { print_check(); return --cnt; }

  // postfix
  inline value_type operator++(int) { print_check(); return cnt++; }
  inline value_type operator--(int) { print_check(); return cnt--; }

  template <typename T>
  inline value_type operator+= (T i) { print_check(); return cnt += i; }
  template <typename T>
  inline value_type operator-= (T i) { print_check(); return cnt -= i; }

  template <typename T>
  inline bool operator== (T i) const noexcept { return cnt == i; }
  template <typename T>
  inline bool operator!= (T i) const noexcept { return cnt != i; }
  template <typename T>
  inline bool operator<  (T i) const noexcept { return cnt <  i; }
  template <typename T>
  inline bool operator<= (T i) const noexcept { return cnt <= i; }
  template <typename T>
  inline bool operator>  (T i) const noexcept { return cnt >  i; }
  template <typename T>
  inline bool operator>= (T i) const noexcept { return cnt >= i; }

  // cast to integral type
  template <typename T> inline operator T () {
    static_assert( std::is_integral<T>::value,
      "Cannot cast timed_counter to a non-integral type" );
    return cnt;
  }

  friend inline std::ostream&
  operator<<(std::ostream& os, const timed_counter& tc)
  noexcept(noexcept(os << tc.cnt)) {
    return (os << tc.cnt);
  }
};

template <typename T, typename L>
const std::locale timed_counter<T,L>::cnt_locale(
  std::locale(), new comma_numpunct() );

} // end namespace

#endif
