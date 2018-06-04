#ifdef HIST_GLOBAL // ===============================================

#define BINS_FILE "config/Hjets_ang.bins"

MAKE_ENUM(H_y_cat,(all)(central))

#define CATEGORIES (H_y_cat)(photon_cuts)(isp)

#endif
#ifdef HIST_DEFS // ============================================

h_(H_y)
h_(H_cosTheta)
h_(Hj_mass)

h_(H_absCosTheta,Hj_mass)

#endif
#ifdef HIST_VARS // ==================================================

const auto H_y = higgs->Rapidity();

cat_bin::id<H_y_cat>() = ( std::abs(H_y) < 0.1 );

#endif
#ifdef HIST_FILL // ============================================

h_H_y(H_y);

const TLorentzVector jet1(
  fj_jets.front()[0],
  fj_jets.front()[1],
  fj_jets.front()[2],
  fj_jets.front()[3]
);
const auto Q = *higgs + jet1;
const double Q2 = Q*Q;
const double M = sqrt(Q2);

const TLorentzVector Z(0,0,Q.E(),Q.Pz());
const auto ell = ((Q*jet1)/Q2)**higgs - ((Q**higgs)/Q2)*jet1;

const double cos_theta = (ell*Z) / sqrt(sq(ell)*sq(Z));

h_H_cosTheta(cos_theta);
h_Hj_mass(M);

h_H_absCosTheta_Hj_mass(std::abs(cos_theta),M);

#endif
#ifdef HIST_INFO // =================================================

note("central","|H_y| < 0.1");

#endif
