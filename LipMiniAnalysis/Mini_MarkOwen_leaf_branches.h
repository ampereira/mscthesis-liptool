
   // Declaration of leaf types
   Int_t           runNumber;
   Int_t           eventNumber;
   Int_t           channelNumber;
   Int_t           rndRunNumber;
   Int_t           dataPeriod;
   Float_t         mu;
   Float_t         mcWeight;
   Float_t         eventWeight_PRETAG;
   Float_t         eventWeight_BTAG;
   Int_t           pvxp_n;
   Int_t           hfor;
   Float_t         vxp_z;
   Float_t         mc_vxp_z;
   Double_t        mcevt_pdf1;
   Double_t        mcevt_pdf2;
   Int_t           mcevt_pdf_id1;
   Int_t           mcevt_pdf_id2;
   Double_t        mcevt_pdf_scale;
   Double_t        mcevt_pdf_x1;
   Double_t        mcevt_pdf_x2;
   Float_t         scaleFactor_PILEUP;
   Float_t         scaleFactor_ELE;
   Float_t         scaleFactor_MUON;
   Float_t         scaleFactor_DILEP;
   Float_t         scaleFactor_BTAG;
   Float_t         scaleFactor_WJETSNORM;
   Float_t         scaleFactor_WJETSSHAPE;
   Float_t         scaleFactor_JVFSF;
   Float_t         scaleFactor_ZVERTEX;
   Float_t         scaleFactor_ALLPRETAG;
   Float_t         scaleFactor_ALLBTAG;
   Int_t           truE;
   Int_t           truM;
   Bool_t          trigE;
   Bool_t          trigM;
   Bool_t          passGRL;
   Bool_t          cosmicEvent;
   Bool_t          isOS;
   Bool_t          hasGoodVertex;
   Float_t         ht;
   Float_t         mass;
   Float_t         massInv_LL;
   Float_t         massInv_BQQ;
   Float_t         massTrans_BLMet;
   Float_t         massTrans_LMet;
   UInt_t          flag_DL;
   UInt_t          flag_TTZ;
   UInt_t          channel_DL;
   UInt_t          channel_TTZ;
   Float_t         scaledWeight;
   UInt_t          lep_n;
   Bool_t          lep_truthMatched[3];   //[lep_n]
   Bool_t          lep_trigMatched[3];   //[lep_n]
   Float_t         lep_pt[3];   //[lep_n]
   Float_t         lep_eta[3];   //[lep_n]
   Float_t         lep_phi[3];   //[lep_n]
   Float_t         lep_E[3];   //[lep_n]
   Float_t         lep_z0[3];   //[lep_n]
   Float_t         lep_charge[3];   //[lep_n]
   Bool_t          lep_isTight[3];   //[lep_n]
   UInt_t          lep_type[3];   //[lep_n]
   UInt_t          lep_flag[3];   //[lep_n]
   Float_t         el_cl_eta[3];   //[lep_n]
   Float_t         lep_ptcone30[3];   //[lep_n]
   Float_t         lep_etcone20[3];   //[lep_n]
   Float_t         massTrans_LMet_Vec[3];   //[lep_n]
   UInt_t          lepPair_n;
   Bool_t          isSameFlavor_LL_Vec[3];   //[lepPair_n]
   Bool_t          isOppSign_LL_Vec[3];   //[lepPair_n]
   Float_t         massInv_LL_Vec[3];   //[lepPair_n]
   Float_t         met_sumet;
   Float_t         met_et;
   Float_t         met_phi;
   UInt_t          jet_n;
   UInt_t          alljet_n;
   Float_t         jet_pt[14];   //[alljet_n]
   Float_t         jet_eta[14];   //[alljet_n]
   Float_t         jet_phi[14];   //[alljet_n]
   Float_t         jet_E[14];   //[alljet_n]
   Float_t         jet_jvf[14];   //[alljet_n]
   Int_t           jet_trueflav[14];   //[alljet_n]
   Int_t           jet_truthMatched[14];   //[alljet_n]
   Float_t         jet_SV0[14];   //[alljet_n]
   Float_t         jet_MV1[14];   //[alljet_n]

   // List of branches
   TBranch        *b_runNumber;   //!
   TBranch        *b_eventNumber;   //!
   TBranch        *b_channelNumber;   //!
   TBranch        *b_rndRunNumber;   //!
   TBranch        *b_dataPeriod;   //!
   TBranch        *b_mu;   //!
   TBranch        *b_mcWeight;   //!
   TBranch        *b_eventWeight_PRETAG;   //!
   TBranch        *b_eventWeight_BTAG;   //!
   TBranch        *b_pvxp_n;   //!
   TBranch        *b_hfor;   //!
   TBranch        *b_vxp_z;   //!
   TBranch        *b_mc_vxp_z;   //!
   TBranch        *b_mcevt_pdf1;   //!
   TBranch        *b_mcevt_pdf2;   //!
   TBranch        *b_mcevt_pdf_id1;   //!
   TBranch        *b_mcevt_pdf_id2;   //!
   TBranch        *b_mcevt_pdf_scale;   //!
   TBranch        *b_mcevt_pdf_x1;   //!
   TBranch        *b_mcevt_pdf_x2;   //!
   TBranch        *b_m_scaleFactor_PILEUP;   //!
   TBranch        *b_m_scaleFactor_ELE;   //!
   TBranch        *b_m_scaleFactor_MUON;   //!
   TBranch        *b_m_scaleFactor_DILEP;   //!
   TBranch        *b_m_scaleFactor_BTAG;   //!
   TBranch        *b_m_scaleFactor_WJETSNORM;   //!
   TBranch        *b_m_scaleFactor_WJETSSHAPE;   //!
   TBranch        *b_m_scaleFactor_JVFSF;   //!
   TBranch        *b_m_scaleFactor_ZVERTEX;   //!
   TBranch        *b_m_scaleFactor_ALLPRETAG;   //!
   TBranch        *b_m_scaleFactor_ALLBTAG;   //!
   TBranch        *b_truE;   //!
   TBranch        *b_truM;   //!
   TBranch        *b_trigE;   //!
   TBranch        *b_trigM;   //!
   TBranch        *b_passGRL;   //!
   TBranch        *b_cosmicEvent;   //!
   TBranch        *b_isOS;   //!
   TBranch        *b_hasGoodVertex;   //!
   TBranch        *b_ht;   //!
   TBranch        *b_mass;   //!
   TBranch        *b_massInv_LL;   //!
   TBranch        *b_massInv_BQQ;   //!
   TBranch        *b_massTrans_BLMet;   //!
   TBranch        *b_massTrans_LMet;   //!
   TBranch        *b_flag_DL;   //!
   TBranch        *b_flag_TTZ;   //!
   TBranch        *b_channel_DL;   //!
   TBranch        *b_channel_TTZ;   //!
   TBranch        *b_scaledWeight;   //!
   TBranch        *b_lep_n;   //!
   TBranch        *b_lep_truthMatched;   //!
   TBranch        *b_lep_trigMatched;   //!
   TBranch        *b_lep_pt;   //!
   TBranch        *b_lep_eta;   //!
   TBranch        *b_lep_phi;   //!
   TBranch        *b_lep_E;   //!
   TBranch        *b_lep_z0;   //!
   TBranch        *b_lep_charge;   //!
   TBranch        *b_lep_isTight;   //!
   TBranch        *b_lep_type;   //!
   TBranch        *b_lep_flag;   //!
   TBranch        *b_el_cl_eta;   //!
   TBranch        *b_lep_ptcone30;   //!
   TBranch        *b_lep_etcone20;   //!
   TBranch        *b_massTrans_LMet_Vec;   //!
   TBranch        *b_lepPair_n;   //!
   TBranch        *b_isSameFlavor_LL_Vec;   //!
   TBranch        *b_isOppSign_LL_Vec;   //!
   TBranch        *b_massInv_LL_Vec;   //!
   TBranch        *b_met_sumet;   //!
   TBranch        *b_met_et;   //!
   TBranch        *b_met_phi;   //!
   TBranch        *b_jet_n;   //!
   TBranch        *b_alljet_n;   //!
   TBranch        *b_jet_pt;   //!
   TBranch        *b_jet_eta;   //!
   TBranch        *b_jet_phi;   //!
   TBranch        *b_jet_E;   //!
   TBranch        *b_jet_jvf;   //!
   TBranch        *b_jet_trueflav;   //!
   TBranch        *b_jet_truthMatched;   //!
   TBranch        *b_jet_SV0;   //!
   TBranch        *b_jet_MV1;   //!
