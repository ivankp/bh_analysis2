#include <cmath>

#include <LHAPDF/LHAPDF.h>

#include "reweighter.hh"

constexpr std::array<int,10> quarks {
  1,-1, 2,-2, 3,-3, 4,-4, 5,-5
};

entry::entry(TTree& tree, Long64_t cacheSize) {
  tree.SetCacheSize(cacheSize);
  tree.SetBranchStatus("*",0);
  branch(tree, "nparticle",   &nparticle);
  branch(tree, "px",           px);
  branch(tree, "py",           py);
  branch(tree, "pz",           pz);
  branch(tree, "E",            E);
  branch(tree, "kf",           kf);
  branch(tree, "alphas",      &alphas);
  branch(tree, "weight2",     &weight2);
  branch(tree, "me_wgt",      &me_wgt);
  branch(tree, "me_wgt2",     &me_wgt2);
  branch(tree, "x1",          &x[0]);
  branch(tree, "x2",          &x[1]);
  branch(tree, "x1p",         &xp[0]);
  branch(tree, "x2p",         &xp[1]);
  branch(tree, "id1",         &id[0]);
  branch(tree, "id2",         &id[1]);
  branch(tree, "fac_scale",   &fac_scale);
  branch(tree, "ren_scale",   &ren_scale);
  branch(tree, "usr_wgts",     usr_wgts);
  branch(tree, "alphasPower", &alphas_power);
  branch(tree, "part",         part);
  tree.StopCacheLearningPhase();
}

inline std::vector<std::unique_ptr<LHAPDF::PDF>>
mkpdfs(const std::string& name, bool all) {
  if (all) {
    const auto pdfs = LHAPDF::mkPDFs(name);
    return { pdfs.begin(), pdfs.end() };
  } else {
    std::vector<std::unique_ptr<LHAPDF::PDF>> pdfs;
    pdfs.emplace_back( LHAPDF::mkPDF(name) );
    return std::move(pdfs);
  }
}

reweighter::reweighter(
  TTree& tree, const scale_defs& sd,
  const std::string& pdf_name, bool all,
  Long64_t cacheSize
) : entry(tree,cacheSize), _pdfs(mkpdfs(pdf_name,all)),
    sd(sd), scale_values(sd.scale_fcns.size()),
    new_weights(sd.scales.size()+_pdfs.size()-1),
    _fac_vars(sd.scales_fac.size()), _ren_vars(sd.scales_ren.size())
{ }
reweighter::~reweighter() { }

double reweighter::fr1(unsigned r, double muF) const {
  const double x = this->x[r];
  if (id[r] != 21) return pdf->xfxQ(id[r], x, muF)/x;
  else {
    double f = 0.;
    for (int q : quarks) f += pdf->xfxQ(q, x, muF)/x;
    return f;
  }
};

double reweighter::fr2(unsigned r, double muF) const {
  const double x  = this->x[r];
  const double xp = x/this->xp[r];
  if (id[r] != 21) return pdf->xfxQ(id[r], xp, muF)/x;
  else {
    double f = 0.;
    for (int q : quarks) f += pdf->xfxQ(q, xp, muF)/x;
    return f;
  }
};

double reweighter::fr3(unsigned r, double muF) const {
  return pdf->xfxQ(21, x[r], muF)/x[r];
}
double reweighter::fr4(unsigned r, double muF) const {
  return pdf->xfxQ(21, x[r]/xp[r], muF)/x[r];
}

void reweighter::fac_calc(unsigned i) {
  fac_vars& v = _fac_vars[i];
  const double muF = scale_values[sd.scales_fac[i]];

  double f[2];
  for (int i=0; i<2; ++i) f[i] = pdf->xfxQ(id[i], x[i], muF)/x[i];
  v.ff = f[0]*f[1];

  if (part[0]=='I') {
    const double lf = 2.*std::log(muF/fac_scale);
    double w[8];
    for (int i=0; i<8; ++i) w[i] = usr_wgts[i+2] + usr_wgts[i+10]*lf;

    v.m = ( fr1(0,muF)*w[0] + fr2(0,muF)*w[1]
          + fr3(0,muF)*w[2] + fr4(0,muF)*w[3] )*f[1]
        + ( fr1(1,muF)*w[4] + fr2(1,muF)*w[5]
          + fr3(1,muF)*w[6] + fr4(1,muF)*w[7] )*f[0];
  } else v.m = 0.;
}

void reweighter::ren_calc(unsigned i) {
  ren_vars& v = _ren_vars[i];
  const double muR = scale_values[sd.scales_ren[i]];

  v.ar = std::pow(pdf->alphasQ(muR)/alphas, alphas_power);

  if (part[0]=='V' || part[0]=='I') {
    const double lr = 2.*std::log(muR/ren_scale);
    v.w0 = me_wgt + lr*usr_wgts[0] + 0.5*lr*lr*usr_wgts[1];
  } else {
    v.w0 = me_wgt2;
  }
}

double reweighter::combine(unsigned i) { // i is scale index
  const auto& scale = sd.scales[i];
  double w;

  if (scale.fac) {
    const fac_vars& fac = _fac_vars[*scale.fac];
    w = fac.m;
    if (scale.ren) {
      const ren_vars& ren = _ren_vars[*scale.ren];
      w += ren.w0 * fac.ff;
    } else {
      w += me_wgt2 * fac.ff;
    }
  } else {
    w = weight2;
  }
  if (scale.ren) {
    const ren_vars& ren = _ren_vars[*scale.ren];
    w *= ren.ar;
  }

  return w;
}

void reweighter::operator()() {
  for (unsigned i=0; i<sd.scale_fcns.size(); ++i)
    scale_values[i] = sd.scale_fcns[i](*this);

  pdf = _pdfs[0].get();
  for (unsigned i=0; i<sd.scales_fac.size(); ++i) fac_calc(i);
  for (unsigned i=0; i<sd.scales_ren.size(); ++i) ren_calc(i);

  // scale variations
  for (unsigned i=0; i<sd.scales.size(); ++i)
    new_weights[i] = combine(i);

  // pdf variations
  for (unsigned i=1; i<_pdfs.size(); ++i) {
    pdf = _pdfs[i].get();
    fac_calc(0);
    ren_calc(0);

    new_weights[sd.scales.size()+i] = combine(0);
  }
}

const int reweighter::pdf_id(unsigned i) const {
  const auto *pdf = _pdfs[i].get();
  return pdf->lhapdfID() + pdf->memberID();
}

