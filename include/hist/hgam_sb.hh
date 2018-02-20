#ifdef IMPL_GLOBAL // ===============================================

#define BINS_FILE "config/sb.bins"

#endif
#ifdef IMPL_HIST_DEFS // ============================================

h_(m_yy) h_(pT_yy)
h_(pT_yy_105_160) h_(pT_yy_121_129)
hj_(pT)
h_(cosTS_yy) h_(cos_y1) h_(cos_y2)

h_(m_yy,pT_yy)
h_(m_yy,pT_j1)

h_(dR_yy,pT_yy)
h_(cosTS_yy,pT_yy)
h_(cos_y1,pT_yy)

h_(m_yyj,eta_yy) h_(m_yyj,eta_j1)

#endif
#ifdef IMPL_VARS // ==================================================

#endif
#ifdef IMPL_HIST_FILL // ============================================

const double H_mass = higgs->M();
const double H_pT = higgs->Pt();
const auto jets_pT = fj_jets | [](const auto& jet){ return jet.pt(); };

const double cosTS_yy =
  std::sinh(std::abs(A1_eta-A2_eta)) * A1_pT * A2_pT * 2
  / ( std::sqrt(1.+sq(H_pT/H_mass)) * sq(H_mass) );

TLorentzVector AAf = *higgs, A1f = *A1, A2f = *A2;
const auto boost = -AAf.BoostVector();
AAf.Boost(boost);
A1f.Boost(boost);
A2f.Boost(boost);
const double A1_cos = std::abs(A1f.CosTheta());
const double A2_cos = std::abs(A2f.CosTheta());

for (unsigned i=0, n=std::min(njets,njets_expected+1); i<n; ++i) {
  h_jet_pT[i](jets_pT[i]);
}

h_m_yy(H_mass);
h_pT_yy(H_pT);

if (105.<H_mass && H_mass<160.) {
  h_pT_yy_105_160(H_pT);
  if (121.<H_mass && H_mass<129.)
    h_pT_yy_121_129(H_pT);
}

h_cosTS_yy(cosTS_yy);
h_cos_y1(A1_cos);
h_cos_y2(A2_cos);

h_m_yy_pT_yy(H_mass,H_pT);
h_m_yy_pT_j1(H_mass,jets_pT[0]);

h_dR_yy_pT_yy(diphoton.first.DeltaR(diphoton.second),H_pT);

h_cosTS_yy_pT_yy(cosTS_yy,H_pT);
h_cos_y1_pT_yy(A1_cos,H_pT);

const auto Hj = *higgs + fj_jets[0];
const double Hj_mass = Hj.M();
const double H_eta = higgs->Eta();
const double j1_eta = fj_jets[0].eta();

h_m_yyj_eta_yy(Hj_mass,H_eta);
h_m_yyj_eta_j1(Hj_mass,j1_eta);

#endif
#ifdef IMPL_INFO // =================================================

#endif
