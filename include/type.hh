#ifndef IVANP_TYPE_HH
#define IVANP_TYPE_HH

#include "literal.hh"

// https://stackoverflow.com/a/20170989/2640636

template <typename T>
constexpr literal type_str() {
#ifdef __clang__
  literal p = __PRETTY_FUNCTION__;
  return literal(p.data() + 24, p.size() - 25);
#elif defined(__GNUC__)
  literal p = __PRETTY_FUNCTION__;
# if __cplusplus < 201402
  return literal(p.data() + 29, p.size() - 30);
# else
  return literal(p.data() + 39, p.size() - 40);
# endif
#else
# error type function does not work with this compiler
#endif
}

template <typename T> void prt_type() {
  std::cout << "\033[34;1m" << type_str<T>() << "\033[0m ("
            << sizeof(T) << ')' << std::endl;
}

#endif
