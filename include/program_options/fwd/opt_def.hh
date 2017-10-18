#ifndef IVANP_OPT_DEF_FWD_HH
#define IVANP_OPT_DEF_FWD_HH

namespace ivanp { namespace po {
namespace detail {

struct opt_def {
  std::string name, descr;
  unsigned count = 0;

  opt_def(std::string&& descr): descr(std::move(descr)) { }
  virtual ~opt_def() { }
  virtual void parse(const char* arg) = 0;
  virtual void as_switch() = 0;
  virtual void default_init() = 0;

  virtual bool is_switch() const noexcept = 0;
  virtual bool is_multi() const noexcept = 0;
  virtual bool is_pos() const noexcept = 0;
  virtual bool is_pos_end() const noexcept = 0;
  virtual bool is_req() const noexcept = 0;
  virtual bool is_signed() const noexcept = 0;
  virtual bool is_named() const noexcept = 0;
  virtual bool is_switch_init() const noexcept = 0;

  template <typename Props>
  inline void set_name(Props&&, nothing) noexcept { }
  template <typename Props, size_t I>
  inline void set_name(Props&& props,
    just<std::integral_constant<size_t,I>>
  ) noexcept { name = std::move(std::get<I>(std::move(props)).name); }
};

}
}}

#endif
