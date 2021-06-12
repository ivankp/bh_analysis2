#ifdef HIST_GLOBAL // ===============================================

#define BINS_FILE "config/Hjets_ang.bins"

struct ang_t {
  double M, cos_theta;
};
ang_t ang(const TLorentzVector& higgs, const TLorentzVector& jet) {
  const auto Q = higgs + jet;
  const double Q2 = Q*Q;

  const TLorentzVector Z(0,0,Q.E(),Q.Pz());
  const auto ell = ((Q*jet)/Q2)*higgs - ((Q*higgs)/Q2)*jet;

  return {
    sqrt(Q2),
    (ell*Z) / sqrt(sq(ell)*sq(Z))
  };
};

MAKE_ENUM(H_rapidity_cut,(all)(central_higgs))
MAKE_ENUM(fat_jet,(all)(nsubjets_1)(nsubjets_2))

#define CATEGORIES (isp)(photon_cuts)(H_rapidity_cut)(fat_jet)

#endif
#ifdef HIST_DEFS // =================================================

h_(j1_nsub)

h_(H_y)
h_(H_cosTheta)
h_(Hj_mass)

h_(j1_mass)
h_(j2_mass)

h_(j1_x)

h_(H_cosTheta,Hj_mass)

h_(Hj_mass,j1_x)
h_(j1_mass,j1_x)

h_(H_pT,Hj_mass)
h_(j1_pT,Hj_mass)
h_(j2_pT,Hj_mass)

h_(pp_pT_rat,Hj_mass)
h_(pp_dphi,Hj_mass)

#endif
#ifdef HIST_VARS // =================================================

std::sort( particles.begin(), particles.end(),
  [](const fj::PseudoJet& a, const fj::PseudoJet& b){
    return ( a.pt() > b.pt() );
  });

const auto H_y = higgs->Rapidity();
cat_bin::id<H_rapidity_cut>() = ( std::abs(H_y) < 0.1 );

#endif
#ifdef HIST_FILL // =================================================

const auto j1_nsub = fj_jets[0].constituents().size();
cat_bin::id<fat_jet>() = j1_nsub ;

h_H_y(H_y);

const auto ang_Hj1 = ang(*higgs,
  { fj_jets[0][0], fj_jets[0][1], fj_jets[0][2], fj_jets[0][3] }
);

h_H_cosTheta(ang_Hj1.cos_theta);
h_Hj_mass(ang_Hj1.M);

h_H_cosTheta_Hj_mass(ang_Hj1.cos_theta,ang_Hj1.M);

h_j1_nsub(j1_nsub);

const double j1_mass = fj_jets[0].m();
h_j1_mass(j1_mass);

const double j1_x = j1_mass/(jet_def.R()*fj_jets[0].pt());
h_j1_x(j1_x);

h_Hj_mass_j1_x(ang_Hj1.M,j1_x);
h_j1_mass_j1_x(j1_mass,j1_x);



h_H_pT_Hj_mass(higgs->Pt(),ang_Hj1.M);
h_j1_pT_Hj_mass(fj_jets[0].pt(),ang_Hj1.M);

h_pp_pT_rat_Hj_mass(particles[0].pt()/particles[1].pt(),ang_Hj1.M);
h_pp_dphi_Hj_mass(std::abs(particles[0].delta_phi_to(particles[1])),ang_Hj1.M);


if (fj_jets.size() < 2) continue; // --------------------------------

h_j2_mass(fj_jets[1].m());

#endif
#ifdef HIST_INFO // =================================================

note("central_higgs","|H_y| < 0.1");

#endif
