#ifndef IVANP_FLOAT_OR_DOUBLE_READER_HH
#define IVANP_FLOAT_OR_DOUBLE_READER_HH

#include <TLeaf.h>

#include "exception.hh"

bool branch_is_double(TTree* t, const char* branchname) {
  const char *branchtype = t->GetLeaf(branchname)->GetTypeName();
  if (!strcmp(branchtype,"Double_t")) return true;
  else if (!strcmp(branchtype,"Float_t")) return false;
  else throw ivanp::exception(
    "The type of branch ",branchname," is ",branchtype);
}

class float_or_double_value_reader {
public:
  using double_reader_type = TTreeReaderValue<Double_t>;
  using float_reader_type = TTreeReaderValue<Float_t>;

private:
  bool is_double;
  char data[ sizeof(double_reader_type) > sizeof(float_reader_type)
           ? sizeof(double_reader_type) : sizeof(float_reader_type) ];
  inline double_reader_type* d_ptr() noexcept {
    return reinterpret_cast<double_reader_type*>(data);
  }
  inline float_reader_type* f_ptr() noexcept {
    return reinterpret_cast<float_reader_type*>(data);
  }
  inline const double_reader_type* d_ptr() const noexcept {
    return reinterpret_cast<const double_reader_type*>(data);
  }
  inline const float_reader_type* f_ptr() const noexcept {
    return reinterpret_cast<const float_reader_type*>(data);
  }

public:
  float_or_double_value_reader(
    TTreeReader& tr, const char* branchname
  ) : is_double(branch_is_double(tr.GetTree(),branchname))
  {
    if (is_double) new(data) double_reader_type(tr,branchname);
    else new(data) float_reader_type(tr,branchname);
  }

  ~float_or_double_value_reader() {
    if (is_double) d_ptr()->~double_reader_type();
    else           f_ptr()->~float_reader_type();
  }

  inline double operator*() noexcept {
    return (is_double ? **d_ptr() : **f_ptr());
  }

  const char* GetBranchName() const noexcept {
    return (is_double ? d_ptr()->GetBranchName() : f_ptr()->GetBranchName());
  }
};

class float_or_double_array_reader {
public:
  using double_reader_type = TTreeReaderArray<Double_t>;
  using float_reader_type = TTreeReaderArray<Float_t>;

private:
  bool is_double;
  char data[ sizeof(double_reader_type) > sizeof(float_reader_type)
           ? sizeof(double_reader_type) : sizeof(float_reader_type) ];
  inline double_reader_type* d_ptr() noexcept {
    return reinterpret_cast<double_reader_type*>(data);
  }
  inline float_reader_type* f_ptr() noexcept {
    return reinterpret_cast<float_reader_type*>(data);
  }
  inline const double_reader_type* d_ptr() const noexcept {
    return reinterpret_cast<const double_reader_type*>(data);
  }
  inline const float_reader_type* f_ptr() const noexcept {
    return reinterpret_cast<const float_reader_type*>(data);
  }

public:
  float_or_double_array_reader(
    TTreeReader& tr, const char* branchname
  ) : is_double(branch_is_double(tr.GetTree(),branchname))
  {
    if (is_double) new(data) double_reader_type(tr,branchname);
    else new(data) float_reader_type(tr,branchname);
  }

  ~float_or_double_array_reader() {
    if (is_double) d_ptr()->~double_reader_type();
    else           f_ptr()->~float_reader_type();
  }

  inline double operator[](size_t i) noexcept {
    return (is_double ? (*d_ptr())[i] : (*f_ptr())[i]);
  }

  const char* GetBranchName() const noexcept {
    return (is_double ? d_ptr()->GetBranchName() : f_ptr()->GetBranchName());
  }
};

#endif
