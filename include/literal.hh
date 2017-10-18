#ifndef IVANP_LITERAL_HH
#define IVANP_LITERAL_HH

#include <stdexcept>

class literal {
  const char* const str;
  const size_t n;

public:
  typedef const char* const_iterator;

  template <size_t N>
  constexpr literal(const char(&a)[N]) noexcept : str(a), n(N-1) { }
  constexpr literal(const char* str, size_t n) noexcept : str(str), n(n) { }

  constexpr const char* data() const noexcept { return str; }
  constexpr size_t size() const noexcept { return n; }

  constexpr const_iterator begin() const noexcept { return str; }
  constexpr const_iterator end()   const noexcept { return str + n; }

  constexpr char operator[](size_t n) const {
    return n < this->n ? str[n] : throw std::out_of_range("literal");
  }
};

inline std::ostream& operator<<(std::ostream& os, const literal& s) {
  return os.write(s.data(), s.size());
}

#endif
