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

#define MiniTTHReader_cxx
#include "MiniTTHReader.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>

#include <cmath>

using std::cout ;
using std::endl ;

// #############################################################################
MiniTTHReader::MiniTTHReader(int i_isData) {
// #############################################################################
//
//  Purpose: parameters and tools setting
//
//  authors: fveloso
//  first version: 21.fev.2007
//
//  last change: 10.Nov.2012
//  by: A.Onofre
//
// #############################################################################

  // :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  // parameters
  // :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  FillPhotons   = false;
  FillElectrons = true;
  FillJets      = true;
  FillMuons     = true;

  //i_isdata = 1, mc = 0, af2 = -1
  if (i_isData == -1) {
     isData = 0;
     isAF2 = 1;
  } else {
     isData = i_isData;
     isAF2 = 0;
  }

  // AO 19 Nov 2012 ==== Make sure you are reading the right TChain (in this case a mini ntuple)
  tree  = new TChain("mini");

  // :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  // tool initialisation:
  // :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // !!!!!!!! No tools are initialized in this class !!!!!!!!!!!!!
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Init();
}

// #############################################################################
MiniTTHReader::~MiniTTHReader() {
// #############################################################################
//
//  Purpose: delete variables and make sure things are properly handled
//
//  authors: fveloso
//  first version: 21.fev.2007
//
//  last change: 10.Nov.2012
//  by: A.Onofre
//
// #############################################################################

  delete tree;
}


// #############################################################################
void MiniTTHReader::Init() {
// #############################################################################
//
//  Purpose: 	The Init() function is called when the selector needs to initialize
//  		a new tree or chain. Typically here the branch addresses and branch
//  		pointers of the tree will be set.
//  		It is normaly not necessary to make changes to the generated
//  		code, but the routine can be extended by the user if needed.
//  		Init() will be called many times when running on PROOF
//  		(once per file to be processed).
//
//  authors: fveloso
//  first version: 21.fev.2007
//
//  last change: 10.Nov.2012
//  by: A.Onofre
//
// #############################################################################


   // Set object pointer
   // ========================================
   // -- samor 9.Nov.2012 (following lines) -- 
   // --   mini-ntuples leafs and branches  --
   // ========================================
     #include "Mini_MarkOwen_set_pointers.h"
   // ========================================  
   // -- samor 9.Nov.2012 (above lines)     -- 
   // ========================================

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   // ========================================
   // -- samor 9.Nov.2012 (following lines) -- 
   // --   mini-ntuples leafs and branches  --
   // ========================================
     #include "Mini_MarkOwen_fChains.h"
   // ========================================  
   // -- samor 9.Nov.2012 (above lines)     -- 
   // ========================================

   Notify();
}

// #############################################################################
void MiniTTHReader::Input(char *file){
// #############################################################################
//
//  purpose: open nTuple files
//
//  authors: fveloso
//  first version: 21.fev.2007
//
//  last change: 10.Nov.2012
//  by: A.Onofre
//
// #############################################################################

  tree->Add(file);

}

// #############################################################################
Int_t MiniTTHReader::Isub() {
// #############################################################################
//
//  Purpose: return data Run number or MC process number
//
//  authors: fveloso
//  first version: 21.fev.2007
//
//  last change: 10.Nov.2012
//  by: A.Onofre
//
// #############################################################################

  if (isData) return runNumber; 
  else return channelNumber;

}

// #############################################################################
Int_t MiniTTHReader::LumiBlock() {
// #############################################################################
//
//  Purpose: return lumi block
//
//  authors: fveloso
//  first version: 21.fev.2007
//
//  last change: 10.Nov.2012
//  by: A.Onofre
//
// #############################################################################

  int myLumi = -999;
  return myLumi;

}

// #############################################################################
Int_t MiniTTHReader::RunNumber() {
// #############################################################################
//
//  Purpose: return run number
//
//  authors: A.Onofre
//  first version: 14.nov.2012
//
//  last change: 14.Nov.2012
//  by: A.Onofre
//
// #############################################################################

  return runNumber;

}

// #############################################################################
Int_t MiniTTHReader::EveNumber() {
// #############################################################################
//
//  Purpose: return event number
//
//  authors: A.Onofre
//  first version: 14.nov.2012
//
//  last change: 14.Nov.2012
//  by: A.Onofre
//
// #############################################################################

  return eventNumber;

}

// #############################################################################
Int_t MiniTTHReader::TruthEleNumber() {
// #############################################################################
//
//  Purpose: return number of TRUTH electrons
//
//  authors: A.Onofre
//  first version: 14.nov.2012
//
//  last change: 14.Nov.2012
//  by: A.Onofre
//
// #############################################################################

  return truE;

}

// #############################################################################
Int_t MiniTTHReader::TruthMuonNumber() {
// #############################################################################
//
//  Purpose: return number of TRUTH muons
//
//  authors: A.Onofre
//  first version: 14.nov.2012
//
//  last change: 14.Nov.2012
//  by: A.Onofre
//
// #############################################################################

  return truM;

}

// #############################################################################
Int_t MiniTTHReader::ElectronTrigger() {
// #############################################################################
//
//  Purpose: return electron trigger
//
//  authors: fveloso
//  first version: 21.fev.2007
//
//  last change: 10.Nov.2012
//  by: A.Onofre
//
// #############################################################################

  int myTrigE = 0;
  if ( trigE ) myTrigE = 1;
  return myTrigE;

}

// #############################################################################
Int_t MiniTTHReader::MuonTrigger() {
// #############################################################################
//
//  Purpose: return muon trigger
//
//  authors: fveloso
//  first version: 21.fev.2007
//
//  last change: 10.Nov.2012
//  by: A.Onofre
//
// #############################################################################

  int myTrigM = 0;
  if ( trigM ) myTrigM = 1;
  return myTrigM;

}


// #############################################################################
Int_t MiniTTHReader::Cosmic() {
// #############################################################################
//
//  Purpose: verify if event is a cosmic ray  
//           ( Cosmic = 0, it IS NOT a cosmic event)
//
//  authors: fveloso
//  first version: 21.fev.2007
//
//  last change: 10.Nov.2012
//  by: A.Onofre
//
// #############################################################################

  int isCosmic = 0;
  if ( cosmicEvent ) isCosmic = 1;
  return isCosmic;

}

// #############################################################################
Int_t MiniTTHReader::HforFlag() {
// #############################################################################
//
//  Purpose: veto (special) events in case HF==4  
//
//  authors: A.Onofre
//  first version: 21.Dec.2012
//
//  last change: 27.Dec.2012
//  by: A.Onofre
//
// #############################################################################

  int isHforFlag = 0;
  // check hfor for MC events
  if ( ( isData == 0 ) && ( hfor == 4 ) ){ 
	// ttbar+LF Alpgen
	if ( runNumber == 195847 ) isHforFlag = 1;
  }

  return isHforFlag;

}

// #############################################################################
Double_t MiniTTHReader::Ht_Mini() {
// #############################################################################
//
//  Purpose: get Ht from Minintuple
//
//  authors: A.Onofre
//  first version: 21.Dec.2012
//
//  last change: 30.Dec.2012
//  by: A.Onofre
//
// #############################################################################

  return ht;

}

// #############################################################################
Double_t MiniTTHReader::massInv_LL_Mini() {
// #############################################################################
//
//  Purpose: get mLL from Minintuple
//
//  authors: A.Onofre
//  first version: 21.Dec.2012
//
//  last change: 30.Dec.2012
//  by: A.Onofre
//
// #############################################################################

  return massInv_LL;

}

// #############################################################################
Int_t MiniTTHReader::jet_n_Mini() {
// #############################################################################
//
//  Purpose: get number of jets from Minintuple
//
//  authors: A.Onofre
//  first version: 21.Dec.2012
//
//  last change: 30.Dec.2012
//  by: A.Onofre
//
// #############################################################################

  return jet_n;

}

// #############################################################################
Int_t MiniTTHReader::EleMuoOverlap() {
// #############################################################################
//
//  Purpose: check electron-muon overlap
//           ( EleMuoOverlap = 0, it is OK, no overlap)
//
//  authors: fveloso
//  first version: 21.fev.2007
//
//  last change: 10.Nov.2012
//  by: A.Onofre
//
// #############################################################################

  int isEleMuoOverlap = 0;
  return isEleMuoOverlap;

}

// #############################################################################
Int_t MiniTTHReader::GoodRL() {
// #############################################################################
//
//  Purpose: check event passed a goodrunlist
//           ( GoodRL = 1, it is OK)
//
//  authors: A.Onofre
//  first version: 11.Nov.2012
//
//  last change: 11.Nov.2012
//  by: A.Onofre
//
// #############################################################################

  int myPassGRL = 0;
  if ( passGRL ) myPassGRL = 1;
  return myPassGRL;

}


// #############################################################################
Double_t MiniTTHReader::GeV() {
// #############################################################################

  double myGeV = 1000;
  return myGeV;

}

// #############################################################################
Double_t MiniTTHReader::MissPx() {
// #############################################################################
//
//  Purpose: calculate missing Px
//
//  authors: fveloso
//  first version: 21.fev.2007
//
//  last change: 10.Nov.2012
//  by: A.Onofre
//
// #############################################################################

  return  met_et * cos(met_phi) / GeV();  

}

// #############################################################################
Double_t MiniTTHReader::MissPy() {
// #############################################################################
//
//  Purpose: calculate missing Py
//
//  authors: fveloso
//  first version: 21.fev.2007
//
//  last change: 10.Nov.2012
//  by: A.Onofre
//
// #############################################################################

  return  met_et * sin(met_phi) / GeV();  

}

// #############################################################################
Double_t MiniTTHReader::MissPt() {
// #############################################################################
//
//  Purpose: calculate missing Pt
//
//  authors: fveloso
//  first version: 21.fev.2007
//
//  last change: 10.Nov.2012
//  by: A.Onofre
//
// #############################################################################

  return met_et;

}

// #############################################################################
Double_t MiniTTHReader::Sphericity() {
// #############################################################################

  return 0.;

}

// #############################################################################
Double_t MiniTTHReader::Aplanarity() {
// #############################################################################

  return 0.;

}

// #############################################################################
Double_t MiniTTHReader::Planarity() {
// #############################################################################

  return 0.;

}

// #############################################################################
Double_t MiniTTHReader::Weight() {
// #############################################################################
//
//  Purpose: evaluate event weight. Implement here all scale factors
//
//  authors: fveloso
//  first version: 21.fev.2007
//
//  last change: 10.Nov.2012
//  by: A.Onofre
//
// #############################################################################

  double myWeight = 0;

  if (isData) {
    myWeight = 1.;
  } else {
    myWeight = mcWeight;
  }
  return myWeight;

}

// #############################################################################
void MiniTTHReader::FillTruthVec(){
// #############################################################################
//
//  purpose: fill vectors of TLorentzVectors-with-index for truth particles
//
//  authors: nfcastro, fveloso
//  first version: 30.???.??
//
//  last change: 10.Nov.2012
//  by: A.Onofre
//
// #############################################################################

  Truth0_El_N  = 0;
  Truth0_Mu_N  = 0;
  Truth0_Tau_N = 0;
  Truth0_nuEl_N  = 0;
  Truth0_nuMu_N  = 0;
  Truth0_nuTau_N = 0;
  Truth0_lepTau_N = 0;
  Truth0_elTau_N = 0;
  Truth0_muTau_N = 0;
  Truth0_Nu_El_N  = 0;
  Truth0_Nu_Mu_N  = 0;
  Truth0_Nu_Tau_N = 0;
  Truth0_El_k=-1;
  Truth0_Mu_k=-1;
  Truth0_Tau_k=-1;

  // ===== AO 8 Oct 2010 ===================== below =========  

  Int_t IST=0;       		//...particle status code

  //.................t
  ITQ=0;       			//...line for top quark
  IQ1=0;       			//...counter for top quark
  t.SetPtEtaPhiM(0., 0., 0., 0.);
  //.................tbar
  ITB=0;       			//...line for anti top quark
  IQ2=0;       			//...counter for anti top quark
  tb.SetPtEtaPhiM(0., 0., 0., 0.);

  //.................W+
  IWP=0;       			//...line for W+
  IW1=0;       			//...counter for W+
  Wp.SetPtEtaPhiM(0., 0., 0., 0.);
  //.................W-
  IWN=0;       //...line for W-
  IW2=0;       //...counter for W-
  Wn.SetPtEtaPhiM(0., 0., 0., 0.);

  //.................b
  IBQ=0;       			//...line for b
  IB1=0;       			//...counter for bb
  b.SetPtEtaPhiM(0., 0., 0., 0.);
  //.................bb
  IBB=0;       			//...line for bb
  IB2=0;       			//...counter for bb
  bb.SetPtEtaPhiM(0., 0., 0., 0.);

  //.................s from t
  ISQ=0;       			//...line for s
  IS1=0;       			//...counter for s
  s.SetPtEtaPhiM(0., 0., 0., 0.);
  //.................sb from tbar
  ISB=0;       			//...line for sb
  IS2=0;       			//...counter for sb
  sb.SetPtEtaPhiM(0., 0., 0., 0.);

  //.................dw from t
  IDWQ=0;       		//...line for down
  IDW1=0;       		//...counter for down
  dw.SetPtEtaPhiM(0., 0., 0., 0.);
  //.................dwb from tbar
  IDWB=0;       		//...line for down bar
  IDW2=0;       		//...counter for down bar
  dwb.SetPtEtaPhiM(0., 0., 0., 0.);

  //.................W+->f1f2
  IWPf1=0;       		//...line for f1
  IWPf2=0;       		//...line for f2
  IWPf1_Coun=0;  		//...counter for f1 
  IWPf2_Coun=0;  		//...counter for f2 
  IWPtau_Neu=0;    		//...counter for Tau Neutrinos 
  IWPtau_elNeu=0;  		//...counter for ele Neutrinos (from tau decay)
  IWPtau_muNeu=0;  		//...counter for muo Neutrino  (from tau decay)
  Wpf1.SetPtEtaPhiM(0., 0., 0., 0.);
  Wpf2.SetPtEtaPhiM(0., 0., 0., 0.);
  pdgID_Wp_dw=0;       		//...Code of 1st W+ Daughter
  pdgID_Wp_up=0;       		//...Code of 2nd W+ Daughter
  //.................W-->f1f2
  IWNf1=0;       		//...line for f1
  IWNf2=0;       		//...line for f2
  IWNf1_Coun=0;  		//...counter for f1 
  IWNf2_Coun=0;  		//...counter for f2 
  IWNtau_Neu=0;  	  	//...counter for Tau Neutrinos 
  IWNtau_elNeu=0;  		//...counter for ele Neutrinos (from tau decay)
  IWNtau_muNeu=0;  		//...counter for muo Neutrino  (from tau decay)
  Wnf1.SetPtEtaPhiM(0., 0., 0., 0.);
  Wnf2.SetPtEtaPhiM(0., 0., 0., 0.);
  pdgID_Wn_dw=0;       		//...Code of 1st W- Daughter
  pdgID_Wn_up=0;       		//...Code of 2nd W- Daughter
  //.................Tau->lep
  //
  tru_id_leptauWp = 0;		//...line of tau(W+)->mu+(e+)  
  tru_id_leptauWn = 0;          //...line of tau(W-)->mu-(e-)
  //
  // ===== AO 8 Oct 2010 ===================== above =========

  TruthVec->clear();


}

// #############################################################################
void MiniTTHReader::FillObjects() {
// #############################################################################
//
//  WARNING: READ THE FOLLOWING LINES!!!!!!
//  purpose: although formelly this function was used to fill the
//           electrons and muons vectors, due to pratical reasons it is the
//           only function called: it fills (in this order) the electrons,
//           jets and muons vectors
//
//  authors: nfcastro, fveloso
//  first version: ??.???.??
//
//  last change: 10.Nov.2012
//  by: A.Onofre
//
// #############################################################################

  // :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  // RecoType
  // :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  //RecoType = i_RecoType;
  PhRecoType 	= RecoType 	 - RecoType/100*100;
  PhIsoCut 	= RecoType/100 	 - RecoType/10000*100;
  SystID 	= RecoType/10000 - RecoType/1000000*100;

  // reset output
  use_Photon_Vec.clear();
  use_Electron_Vec.clear();
  use_Muon_Vec.clear();
  use_Jet_Vec.clear();
  LeptonVec->clear();
  JetVec->clear();
  Vtx->clear();

  // usefull variables
  jets25inevent=0;

  // auxiliar stuff
  Int_t kf = 0;
  TLorentzVector Good;
  Float_t etaCl;

  // ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  // photons
  // ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  // Loop over all photons, if any!!!!
  if (FillPhotons) {
	// for the moment do nothing!!! minintuples do not have information
  }

  // ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  // electrons
  // ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  //For 2011:  https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/TopCommonObjects2011#Electrons
  Double_t elpT, elEta, elPhi, elEne;
  Int_t elTruthMatch, elTrigMatch;

  // Loop over all electrons.
  if (FillElectrons) for (Int_t k=0; k<lep_n; ++k){

    // ------------------------------check it is an electron
    if ( fabs(lep_type[k])  != 11 )  continue;  

    // ------------------------------accept electron
    elpT   = lep_pt[k];
    elEta  = lep_eta[k];   
    elPhi  = lep_phi[k];   
    elEne  = lep_E[k];
    Good.SetPtEtaPhiE(elpT, elEta, elPhi, elEne);
    if ( lep_charge[k]>0 ) kf = -11; else kf = 11;
    // truth and trigger match
    elTruthMatch = 0;
    if ( lep_truthMatched[k] ) elTruthMatch = 1 ;    
    elTrigMatch  = 0;
    if ( lep_trigMatched[k] )  elTrigMatch  = 1 ;    
    // creat TLorentzWFlags vector
    TLorentzVectorWFlags GoodWFlags(Good, k, kf, 999., elTruthMatch, elTrigMatch);
    use_Electron_Vec.push_back(GoodWFlags);

  }  // Fill Electrons done


  // ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  // take vertex z information (assumed to be ok)
  // ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  if( hasGoodVertex == 1 ){
    // ------------------------------accept vertex
    TVector3 v( -999., -999., vxp_z);
    Vtx->push_back(v);
  }


  // ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  // jets
  // ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  // for jet energy scale and resolution uncertainties:
  Double_t jetET, jetpT, jetEta, jetPhi, jetEne;
  Int_t jetTruthMatch;

  // Loop over all jets
  if (FillJets) for (Int_t k=0; k<jet_n; ++k)  {

    // ------------------------------if jets have MV1>0.795
    if ( BTagCut != 999  &&  jet_MV1[k] > BTagCut  ) kf = 5; else kf = 1;

    // ------------------------------accept jet 
    jetpT   = jet_pt[k];
    jetEta  = jet_eta[k];
    jetPhi  = jet_phi[k];
    jetEne  = jet_E[k];
    Good.SetPtEtaPhiE(jetpT, jetEta, jetPhi, jetEne );
    // truth and trigger match
    jetTruthMatch = 0;
    if ( jet_truthMatched[k] == 1 ) jetTruthMatch = 1 ;    
    // creat TLorentzWFlags vector
    TLorentzVectorWFlags GoodWFlags(Good, k, kf, 999., jetTruthMatch, -1);
    use_Jet_Vec.push_back(GoodWFlags);

  }  // Fill jets done


  // ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  // muons
  // ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  Double_t mupT, muEta, muPhi, muEne;
  Int_t muTruthMatch, muTrigMatch;

  // loop over all muons
  if (FillMuons) for (Int_t k=0; k<lep_n; ++k){

    // ------------------------------check it is an electron
    if ( fabs(lep_type[k])  != 13 )  continue; 
 
    // ------------------------------accept muon
    mupT   = lep_pt[k];
    muEta  = lep_eta[k];
    muPhi  = lep_phi[k]; 
    muEne  = lep_E[k];
    if ( lep_charge[k]>0 ) kf = -13; else kf = 13;
    Good.SetPtEtaPhiE(mupT, muEta, muPhi, muEne);   
    // truth and trigger match
    muTruthMatch = 0;
    if ( lep_truthMatched[k] ) muTruthMatch = 1 ;    
    muTrigMatch  = 0;
    if ( lep_trigMatched[k] )  muTrigMatch  = 1 ;    
    // creat TLorentzWFlags vector
    TLorentzVectorWFlags GoodWFlags(Good, k, kf, 999., muTruthMatch, muTrigMatch);
    use_Muon_Vec.push_back(GoodWFlags);

  }  // Fill muons done 


  // ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  // Fill objects
  // ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

   // Fill jets
   for (Int_t k=0; k< use_Jet_Vec.size(); ++k) {
        JetVec->push_back(use_Jet_Vec[k]);
        if ( use_Jet_Vec[k].Pt() > 25*GeV() ) jets25inevent++;
   }

   // Fill electrons
   for (Int_t k=0; k< use_Electron_Vec.size(); ++k) {
    	LeptonVec->push_back(use_Electron_Vec[k]);
   }
   // Fill muons
   for (Int_t k=0; k< use_Muon_Vec.size();     ++k) {
	LeptonVec->push_back(use_Muon_Vec[k]);
   }
}


