#ifndef ENUM_TRAITS_HH
#define ENUM_TRAITS_HH

#include <cstring>
#include <string>
#include <array>
#include <stdexcept>

#include <boost/preprocessor/stringize.hpp>
#include <boost/preprocessor/seq/enum.hpp>
#include <boost/preprocessor/seq/for_each.hpp>
#include <boost/preprocessor/seq/transform.hpp>

#ifdef _GLIBCXX_OSTREAM
#define MAKE_ENUM_OSTREAM_OP(NAME) \
  std::ostream& operator<<(std::ostream& os, const NAME::type& val) { \
    return os << enum_traits<NAME>::str(val); \
  }
#else
#define MAKE_ENUM_OSTREAM_OP(NAME) { }
#endif

#ifdef _GLIBCXX_ISTREAM
#define MAKE_ENUM_ISTREAM_OP(NAME) \
  std::istream& operator>>(std::istream& is, NAME::type& val) { \
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

#define PP_STRINGIZE_OP(r, data, elem) BOOST_PP_STRINGIZE(elem)
#define PP_PRECEDE_OP(r, data, elem) data elem

#define ENUM_TO_STR_CASE(r, data, elem) \
  case data elem: return BOOST_PP_STRINGIZE(elem);

#define ENUM_FROM_STR_CASE(r, data, elem) \
  if (!strcmp( str, BOOST_PP_STRINGIZE(elem) )) return data elem; else

#define MAKE_ENUM(NAME, VALUES) \
  struct NAME { enum type { BOOST_PP_SEQ_ENUM(VALUES) }; }; \
  template <> struct enum_traits<NAME> { \
    using type = NAME::type; \
    using all = std::integer_sequence<type, \
      PP_SEQ_TRANSFORM_ENUM(PP_PRECEDE_OP,NAME::,VALUES) >; \
    static constexpr size_t size = BOOST_PP_SEQ_SIZE(VALUES); \
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
      throw std::runtime_error( \
        "string cannot be converted to enum " BOOST_PP_STRINGIZE(NAME) ); \
    } \
    \
    static constexpr std::array<type,size> all_val() noexcept { \
      return { PP_SEQ_TRANSFORM_ENUM(PP_PRECEDE_OP,NAME::,VALUES) }; \
    }; \
    \
    static constexpr std::array<const char*,size> all_str() noexcept { \
      return { PP_SEQ_TRANSFORM_ENUM(PP_STRINGIZE_OP,_,VALUES) }; \
    }; \
  }; \
  MAKE_ENUM_OSTREAM_OP(NAME) \
  MAKE_ENUM_ISTREAM_OP(NAME)

#endif
