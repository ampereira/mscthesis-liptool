// #################################################################################
//
//    _      _         __  __ _       _                        _           _     
//   | |    (_)       |  \/  (_)     (_)     /\               | |         (_)    
//   | |     _ _ __   | \  / |_ _ __  _     /  \   _ __   __ _| |_   _ ___ _ ___ 
//   | |    | | '_ \  | |\/| | | '_ \| |   / /\ \ | '_ \ / _` | | | | / __| / __|
//   | |____| | |_) | | |  | | | | | | |  / ____ \| | | | (_| | | |_| \__ \ \__ \
//   |______|_| .__/  |_|  |_|_|_| |_|_| /_/    \_\_| |_|\__,_|_|\__, |___/_|___/
//            | |                                                 __/ |          
//            |_|                                                |___/           
//                                                                              
// #################################################################################
//                                                                              
//   TopD3PDMaker                               \      | | / |  /    |   |  ||  
//                                               '-___-+-+'  |_/     |   |  |   
//   by Antonio Onofre                          `--___-++--^'|       |   | /|   
//      (antonio.onofre@cern.ch)                       ||    |       |   |' |   
//   date: 23.Nov.2012                          --___  ||    |       | _/|  |   
//                                              ==---:-++___ |  ___--+'  |  |   
//                                                '--_'|'---'+--___  |   |  |   
//                                                    '+-_   |     '-+__ |  |   
//                                              ._.          ._.       ._____   
//                                              | |          | |       | ___ \  
//                                              | |_.        | |       | .___/  
//                                              |___|        |_|       |_|
//      
//   Note:      this code was built from classes developed by 
//              Nuno Castro   (nfcastro@lipc.fis.uc.pt)  and 
//              Filipe Veloso (fveloso@lipc.fis.uc.pt)             
//
//
// #################################################################################

#ifndef  MiniTTHReader_h
#define  MiniTTHReader_h

#include "Ntu.h"
#include <TRandom3.h>

class Top_MET_D3PD;

class MiniTTHReader : public Ntu {
public :

   Int_t itmp;

   // Trigger match
   Int_t PassTriggerMatching;

   // counting truth
   Int_t Truth0_El_N;
   Int_t Truth0_Mu_N;
   Int_t Truth0_Tau_N;
   Int_t Truth0_nuEl_N;
   Int_t Truth0_nuMu_N;
   Int_t Truth0_nuTau_N;
   Int_t Truth0_lepTau_N;
   Int_t Truth0_elTau_N;
   Int_t Truth0_muTau_N;

   Int_t Truth0_Nu_El_N;
   Int_t Truth0_Nu_Mu_N;
   Int_t Truth0_Nu_Tau_N;

   Int_t Truth0_El_k;
   Int_t Truth0_Mu_k;
   Int_t Truth0_Tau_k;

   Int_t Truth0_Whad;
   Int_t Truth0_Wlep;
   Int_t Truth0_Wtau;


  // ===== AO 8 Oct 2010 ===================== below =========
  //.................t
  Int_t ITQ;       	//...line for top quark
  Int_t IQ1;       	//...counter for top quark
  TLorentzVector t;
  //.................tbar
  Int_t ITB;       	//...line for anti top quark
  Int_t IQ2;       	//...counter for anti top quark
  TLorentzVector tb;
  //
  //.................W+
  Int_t IWP;       	//...line for W+
  Int_t IW1;       	//...counter for W+
  TLorentzVector Wp;
  //.................W-
  Int_t IWN;       	//...line for W-
  Int_t IW2;       	//...counter for W-
  TLorentzVector Wn;
  //
  //.................b
  Int_t IBQ;       	//...line for b
  Int_t IB1;       	//...counter for b
  TLorentzVector b;
  //.................bb
  Int_t IBB;       	//...line for bb
  Int_t IB2;       	//...counter for bb
  TLorentzVector bb;
  //
  //.................s
  Int_t ISQ;       	//...line for s
  Int_t IS1;       	//...counter for s
  TLorentzVector s;
  //.................sb
  Int_t ISB;       	//...line for sb
  Int_t IS2;       	//...counter for sb
  TLorentzVector sb;
  //
  //.................down
  Int_t IDWQ;       	//...line for down
  Int_t IDW1;       	//...counter for down
  TLorentzVector dw;
  //.................downbar
  Int_t IDWB;       	//...line for down bar
  Int_t IDW2;       	//...counter down bar
  TLorentzVector dwb;
  //
  //.................W+->f1f2
  Int_t IWPf1;       	//...line for f1
  Int_t IWPf2;       	//...line for f2
  Int_t IWPf1_Coun;  	//...counter for f1
  Int_t IWPf2_Coun;  	//...counter for f2
  Int_t IWPtau_Neu;    	//...counter for Tau Neutrinos
  Int_t IWPtau_elNeu;  	//...counter for ele Neutrinos (from tau decay)
  Int_t IWPtau_muNeu;  	//...counter for muo Neutrino  (from tau decay)
  TLorentzVector Wpf1;
  TLorentzVector Wpf2;
  Int_t pdgID_Wp_dw;    //...Code of 1st W+ Daughter
  Int_t pdgID_Wp_up;  	//...Code of 2nd W+ Daughter
  //.................W-->f1f2
  Int_t IWNf1;       	//...line for f1
  Int_t IWNf2;       	//...line for f2
  Int_t IWNf1_Coun;  	//...counter for f1
  Int_t IWNf2_Coun;  	//...counter for f2
  Int_t IWNtau_Neu;    	//...counter for Tau Neutrinos
  Int_t IWNtau_elNeu;  	//...counter for ele Neutrinos (from tau decay)
  Int_t IWNtau_muNeu;  	//...counter for muo Neutrino  (from tau decay)
  TLorentzVector Wnf1;
  TLorentzVector Wnf2;
  Int_t pdgID_Wn_dw;  	//...Code of 1st W- Daughter
  Int_t pdgID_Wn_up;    //...Code of 2nd W- Daughter
  //................. index for truth particle corresponding to tau->lepton
  Int_t tru_id_leptauWp;          //...line of tau(W+)->mu+(e+)
  Int_t tru_id_leptauWn;          //...line of tau(W-)->mu-(e-)
  //
  // ===== AO 8 Oct 2010 ===================== above =========


  // ========================================
  // -- samor 7.Nov.2012 (following lines) -- 
  // --   mini-ntuples leafs and branches  --
  // ========================================
     #include "Mini_MarkOwen_leaf_branches.h"
  // ========================================  
  // -- samor 7.Nov.2012 (above lines)     -- 
  // ========================================


  // Declare methods in TOPD3PD class (all of them)
  MiniTTHReader(int);
  virtual ~MiniTTHReader();
  virtual void Init();
  virtual void Input(char*);
  virtual Int_t Isub();
  virtual Int_t LumiBlock();
  virtual Int_t RunNumber();
  virtual Int_t EveNumber();
  virtual Int_t TruthEleNumber();
  virtual Int_t TruthMuonNumber();
  virtual Int_t ElectronTrigger();
  virtual Int_t MuonTrigger();
  virtual Int_t Cosmic();
  virtual Int_t HforFlag();
  virtual Double_t Ht_Mini();
  virtual Double_t massInv_LL_Mini();
  virtual Int_t jet_n_Mini();
  virtual Int_t EleMuoOverlap();
  virtual Int_t GoodRL();
  virtual Double_t MissPx();
  virtual Double_t MissPy();
  virtual Double_t MissPt();
  virtual Double_t Weight();
  virtual Double_t Sphericity();
  virtual Double_t Aplanarity();
  virtual Double_t Planarity();
  virtual Double_t GeV();

  virtual void FillTruthVec();
  virtual void FillObjects();

  // parameters
  Bool_t FillPhotons;
  Bool_t FillElectrons;
  Bool_t FillJets;
  Bool_t FillMuons;

  // general
  Int_t PhRecoType;
  Int_t PhIsoCut;
  Int_t SystID;
  TRandom3 rnd;
  Int_t isAF2;

  // trigger matching
  std::vector<Int_t> LeptonMatchTrigger;
  std::vector<Int_t> LeptonMatchL2Trigger;


  // variables for QCd estimation
  Int_t jets25inevent;

  std::vector<TLorentzVectorWFlags> use_Photon_Vec;
  std::vector<TLorentzVectorWFlags> use_Electron_Vec;
  std::vector<TLorentzVectorWFlags> use_Muon_Vec;
  std::vector<TLorentzVectorWFlags> use_Jet_Vec;


};

#endif

