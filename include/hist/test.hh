#ifdef HIST_GLOBAL // ===============================================

#define BINS_FILE "config/Hjets.bins"

// MAKE_ENUM(pt_cut_cat,(all)(AA_pt_50_100))
// MAKE_ENUM(mass_cut_cat,(all)(M_121_129))

// #define CATEGORIES (pt_cut_cat)(mass_cut_cat)(photon_cuts)(isp)

#endif
#ifdef HIST_DEFS // =================================================

/*
h_(A1_rap)
h_(A2_rap)

h_(A1j_dR_unsorted)
h_(A2j_dR_unsorted)
h_(A1j_dR)
h_(A2j_dR)
h_(Hj_dR)

h_(A_dE)

h_(A1j_dphi)
h_(A2j_dphi)
h_(A1j_dy)
h_(A2j_dy)
*/

h_(H_eta)
h_(H_y)

// h_(Hj_mass)
// h_(j1_x)
// h_(Hj_mass,j1_x)

// h_(AA_pT)
// h_(j1_pT)
// h_(j1_y)

// h_(AA_dR)

#endif
#ifdef HIST_VARS // =================================================

/*
const auto pt = higgs->Pt();
const auto mass = higgs->M();

cat_bin::id<pt_cut_cat>() = ( 50 < pt && pt < 100 );
cat_bin::id<mass_cut_cat>() = ( 121 < mass && mass < 129 );
*/

#endif
#ifdef HIST_FILL // =================================================

// if (njets != 1) continue;

// const auto& jet = fj_jets.front();

h_H_eta(higgs->Eta());
h_H_y(higgs->Rapidity());

/*
const double dR = particles.at(0).delta_R(particles.at(1));

if (dR < 0.2) continue;
if (dR > 0.7) continue;

const auto jet = particles[0] + particles[1];

const double M = (*higgs+jet).M();
const double jet_x = jet.m()/(jet_def.R()*jet.pt());

h_Hj_mass(M);
h_j1_x(jet_x);
h_Hj_mass_j1_x(M,jet_x);
*/

/*
h_A1_rap(diphoton.first.Rapidity());
h_A2_rap(diphoton.second.Rapidity());

const TLorentzVector j1(fj_jets[0][0], fj_jets[0][1], fj_jets[0][2], fj_jets[0][3]);
h_A1j_dR_unsorted(diphoton.first.DeltaR(j1));
h_A2j_dR_unsorted(diphoton.second.DeltaR(j1));
h_A1j_dR(A1->DeltaR(j1));
h_A2j_dR(A2->DeltaR(j1));
h_Hj_dR(higgs->DeltaR(j1));

h_A_dE(A1->E() - A2->E());

h_A1j_dphi(dphi(A1->Phi(),j1.Phi()));
h_A2j_dphi(dphi(A2->Phi(),j1.Phi()));

h_A1j_dy(std::abs(A1->Rapidity()-j1.Rapidity()));
h_A2j_dy(std::abs(A2->Rapidity()-j1.Rapidity()));
*/

// h_Hj_mass((diphoton.first+diphoton.second+fj_jets[0]).M());

/*
h_AA_pT(higgs->Pt());
h_j1_pT(fj_jets[0].pt());
h_j1_y(fj_jets[0].rap());
*/

/*
const double HM = higgs->M();
if (HM<121) continue;
if (HM>129) continue;

h_AA_dR(A1->DeltaR(*A2));
*/

#endif
#ifdef HIST_POST // =================================================

// auto h_Hj_mass_j1_x_incl = hist<1,0>("h_Hj_mass-j1_x_incl",h_Hj_mass_j1_x);
// h_Hj_mass_j1_x_incl.integrate_left<1>();

#endif
