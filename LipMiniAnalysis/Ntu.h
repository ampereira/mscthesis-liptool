// #############################################################################
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
// #############################################################################
//                                                                              
//   LipMiniAnalysis                            \      | | / |  /    |   |  ||  
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
//   Note:      this code was built from LipCbrAnalysis developed by                                                       
//              Nuno Castro   (nfcastro@lipc.fis.uc.pt)  and 
//              Filipe Veloso (fveloso@lipc.fis.uc.pt)             
//
//
// #############################################################################
//
//  file: Ntu.h
//
//  purpose: this is the declaration of the Ntu class
//
//  author: fveloso
//  first version: 18.nov.2006
//
//  last change: 04.jul.2007
//  by: fveloso
//
// #############################################################################


#ifndef Ntu_h
#define Ntu_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <vector>
#include "TLorentzVectorWFlags.h"
#include "TVector3.h"

class Ntu {
public:
  TChain *tree;
  TChain *tree2;
  TChain *tree3;
  TChain *tree4;
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain

  Ntu();
  virtual ~Ntu();
  virtual Int_t GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void Init();
  virtual Bool_t Notify();

  std::vector<TLorentzVectorWFlags> * TruthVec;
  std::vector<TLorentzVectorWFlags> * TauVec;
  std::vector<TLorentzVectorWFlags> * LeptonVec;
  std::vector<TLorentzVectorWFlags> * PhotonVec;
  std::vector<TLorentzVectorWFlags> * JetVec;
  std::vector<TVector3> * Vtx;
  double LightJetCalib;
  double BJetCalib;
  double BTagCut;
  double PtCutJet, PtCutEle, PtCutMuo;
  int JetOverlap;
  Int_t RecoType;
  int isData;
  void FillVectors(
    	std::vector<TLorentzVectorWFlags>&, 
	std::vector<TLorentzVectorWFlags>&,
    	std::vector<TLorentzVectorWFlags>&,
 	std::vector<TLorentzVectorWFlags>&,
    	std::vector<TLorentzVectorWFlags>&, 
	std::vector<TVector3>&, 
	double, 
	double, 
	double, 
	int, 
	double, 
	double, 
	double, 
	Int_t);

  virtual void Input(char*);

  virtual int Isub();
  virtual int LumiBlock();
  virtual int RunNumber();
  virtual int EveNumber();
  virtual int TruthEleNumber();
  virtual int TruthMuonNumber();
  virtual int ElectronTrigger();
  virtual int MuonTrigger();
  virtual int Cosmic();
  virtual int HforFlag();
  virtual double Ht_Mini();
  virtual double massInv_LL_Mini();
  virtual int jet_n_Mini();
  virtual int EleMuoOverlap();
  virtual int GoodRL();
  virtual double MissPx();
  virtual double MissPy();
  virtual double MissPt();
  virtual double Weight();
  virtual double Sphericity();
  virtual double Aplanarity();
  virtual double Planarity();
  virtual double GeV();

  virtual void FillTruthVec();
  virtual void FillTauVec();
  virtual void FillObjects();
  virtual void FillPhotonVec();
  virtual void FillJetVec();


};

#endif
