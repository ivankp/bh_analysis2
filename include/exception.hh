#ifndef IVANP_EXCEPTION
#define IVANP_EXCEPTION

#include <string>
#include <exception>

#include "catstr.hh"

namespace ivanp {

class exception : public std::exception {
  std::string str;
public:
  exception() noexcept { }

  template <typename... Args>
  exception(Args&& ...args) noexcept
  : str(cat(std::forward<Args>(args)...)) { };

  exception(const char* cstr) noexcept : str(cstr) { }
  exception& operator= (const char* cstr) noexcept {
    str = cstr;
    return *this;
  }

  exception(const std::string& s) noexcept : str(s) { }
  exception(std::string&& s) noexcept : str(std::move(s)) { }
  exception& operator= (const std::string& s) noexcept {
    str = s;
    return *this;
  }
  exception& operator= (std::string&& s) noexcept {
    str = std::move(s);
    return *this;
  }

  exception(const ivanp::exception& e) noexcept : str(e.str) { }
  exception(ivanp::exception&& e) noexcept : str(std::move(e.str)) { }
  exception& operator= (const ivanp::exception& e) noexcept {
    str = e.str;
    return *this;
  }
  exception& operator= (ivanp::exception&& e) noexcept {
    str = std::move(e.str);
    return *this;
  }

  exception(const std::exception& e) noexcept : str(e.what()) { }
  exception& operator= (const std::exception& e) noexcept {
    str = e.what();
    return *this;
  }

  virtual const char* what() const noexcept {
    return str.c_str();
  }
};

}

#endif
