#ifndef IVANP_FLOAT_OR_DOUBLE_HH
#define IVANP_FLOAT_OR_DOUBLE_HH

class float_or_double {
  bool is_double;
  union { float f; double d; } _v;

public:
  float_or_double(float  f): is_double(false), _v.f(f) { }
  float_or_double(double d): is_double(true ), _v.d(d) { }

  inline double operator*() const noexcept {return (is_double ? _v.d : _v.f);}
  inline operator  double() const noexcept {return (is_double ? _v.d : _v.f);}
  inline operator   float() const noexcept {return (is_double ? _v.d : _v.f);}
};

#endif
