#ifdef HIST_GLOBAL // ===============================================

#define BINS_FILE "config/sb.bins"

#include <boost/preprocessor/seq.hpp>

#define CAT_OP(s, state, x) BOOST_PP_CAT(state##_,x)
#define CAT_(pref, seq) BOOST_PP_SEQ_FOLD_LEFT(CAT_OP, pref, seq)

#define SEQ_m_yy ((105)(160))((117)(133))((121)(129))((124)(126))

#endif
#ifdef HIST_DEFS // ============================================

#define SEQ_h_(r, pref, seq) h_(CAT_(pref,seq))

h_(m_yy) h_(pT_yy)
BOOST_PP_SEQ_FOR_EACH(SEQ_h_, pT_yy, SEQ_m_yy)
hj_(pT)
h_(cos_cs) h_(cos_hr)

h_(m_yy,pT_yy)
h_(m_yy,pT_j1)
h_(m_yy,pT_j2)

h_(dR_yy,pT_yy)
h_(cos_cs,pT_yy)
h_(cos_hr,pT_yy)

h_(m_yyj,eta_yy) h_(m_yyj,eta_j1)

#endif
#ifdef HIST_VARS // ==================================================

#endif
#ifdef HIST_FILL // ============================================

const double H_mass = higgs->M();
const double H_pT = higgs->Pt();
const auto jets_pT = fj_jets | [](const auto& jet){ return jet.pt(); };

const double cos_cs =
  std::sinh(std::abs(A1_eta-A2_eta)) * A1_pT * A2_pT * 2
  / ( std::sqrt(1.+sq(H_pT/H_mass)) * sq(H_mass) );

TLorentzVector A1_hf = *A1;
A1_hf.Boost(-higgs->BoostVector());
const double cos_hr = std::abs(A1_hf.CosTheta());

const unsigned nj = std::min(njets,njets_expected+1);
for (unsigned i=0; i<nj; ++i) {
  h_jet_pT[i](jets_pT[i]);
}

h_m_yy(H_mass);
h_pT_yy(H_pT);

#define h_pT_yy_mass_case(s, state, seq) \
  if (BOOST_PP_SEQ_ELEM(0,seq)<H_mass && H_mass<BOOST_PP_SEQ_ELEM(1,seq)) \
  { CAT_(h_pT_yy,seq)(H_pT); state }

BOOST_PP_SEQ_FOLD_RIGHT(h_pT_yy_mass_case,,SEQ_m_yy)

h_cos_cs(cos_cs);
h_cos_hr(cos_hr);

h_m_yy_pT_yy(H_mass,H_pT);
h_m_yy_pT_j1(H_mass,jets_pT[0]);
if (nj>1) h_m_yy_pT_j2(H_mass,jets_pT[1]);

h_dR_yy_pT_yy(diphoton.first.DeltaR(diphoton.second),H_pT);

h_cos_cs_pT_yy(cos_cs,H_pT);
h_cos_hr_pT_yy(cos_hr,H_pT);

const auto Hj = *higgs + fj_jets[0];
const double Hj_mass = Hj.M();
const double H_eta = higgs->Eta();
const double j1_eta = fj_jets[0].eta();

h_m_yyj_eta_yy(Hj_mass,H_eta);
h_m_yyj_eta_j1(Hj_mass,j1_eta);

#endif
#ifdef HIST_INFO // =================================================

#endif
