#ifndef ENUM_TRAITS_HH
#define ENUM_TRAITS_HH

#include <cstring>
#include <string>
#include <array>
#include <stdexcept>

#include <boost/preprocessor/seq/enum.hpp>
#include <boost/preprocessor/seq/for_each.hpp>
#include <boost/preprocessor/seq/transform.hpp>

// MAKE_ENUM(num,(one)(two)(tree))
//
// num : enum class num { one, two, three };
// enum_traits<num>::type : num
// enum_traits<num>::size : 3
// enum_traits<num>::all : std::integer_sequence<num,one,two,tree>
// enum_traits<num>::str(one) : "one"
// enum_traits<num>::val("one") : num::one
// enum_traits<num>::all_str() : std::array<const char*,size>
// enum_traits<num>::all_val() : std::array<num,size>

#ifdef _GLIBCXX_OSTREAM
#define MAKE_ENUM_OSTREAM_OP(NAME) \
  std::ostream& operator<<(std::ostream& os, const NAME& val) { \
    return os << enum_traits<NAME>::str(val); \
  }
#else
#define MAKE_ENUM_OSTREAM_OP(NAME) { }
#endif

#ifdef _GLIBCXX_ISTREAM
#define MAKE_ENUM_ISTREAM_OP(NAME) \
  std::istream& operator>>(std::istream& is, NAME& val) { \
    std::string str; \
    is >> str; \
    val = enum_traits<NAME>::val(str.c_str()); \
    return is; \
  }
#else
#define MAKE_ENUM_ISTREAM_OP(NAME) { }
#endif

template <typename Enum> struct enum_traits;

#define PP_SEQ_TRANSFORM_ENUM(macro, data, seq) \
  BOOST_PP_SEQ_ENUM( BOOST_PP_SEQ_TRANSFORM(macro, data, seq) )

#define PP_STRINGIZE_OP(r, data, elem) #elem
#define PP_PRECEDE_OP(r, data, elem) data elem

#define ENUM_TO_STR_CASE(r, data, elem) \
  case data elem: return #elem;

#define ENUM_FROM_STR_CASE(r, data, elem) \
  if (!strcmp( str, #elem )) return data elem; else

#define MAKE_ENUM(NAME, VALUES) \
  enum class NAME { BOOST_PP_SEQ_ENUM(VALUES) }; \
  template <> struct enum_traits<NAME> { \
    using type = NAME; \
    using all = std::integer_sequence<type, \
      PP_SEQ_TRANSFORM_ENUM(PP_PRECEDE_OP, NAME::, VALUES) >; \
    static constexpr size_t size = all::size(); \
    \
    static constexpr const char* str(type val) noexcept { \
      switch (val) { \
        BOOST_PP_SEQ_FOR_EACH(ENUM_TO_STR_CASE, NAME::, VALUES) \
        default: return nullptr; \
      } \
    }; \
    \
    static type val(const char* str) { \
      BOOST_PP_SEQ_FOR_EACH( ENUM_FROM_STR_CASE, NAME::, VALUES ) \
      throw std::runtime_error("string cannot be converted to enum " #NAME); \
    } \
    \
    static constexpr std::array<type,size> all_val() noexcept { \
      return { PP_SEQ_TRANSFORM_ENUM(PP_PRECEDE_OP, NAME::, VALUES) }; \
    }; \
    \
    static constexpr std::array<const char*,size> all_str() noexcept { \
      return { PP_SEQ_TRANSFORM_ENUM(PP_STRINGIZE_OP, ,VALUES) }; \
    }; \
  }; \
  MAKE_ENUM_OSTREAM_OP(NAME) \
  MAKE_ENUM_ISTREAM_OP(NAME)

#endif
