
#ifndef ttH_dilep_h
#define ttH_dilep_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <vector>
#include "TLorentzVector.h"
#include <iostream>
#include <TH1.h>
#include <TH2.h>
#include "Riostream.h"
#include "../../LipMiniAnalysis/LipMiniAnalysis.h"

using std::cout;
using std::endl;


// #############################################################################
class ttH_dilep : public LipMiniAnalysis {
// #############################################################################
public :

  ttH_dilep();

  // signal pdf name
//  char* signalpdf;
  char signalpdf[100];

  // amostras para ler
  virtual void DefineSamples();
  virtual void UserValues();
  virtual void UserCommandLineOptions(char*, Int_t);
  virtual void DefineSignalPdfSamples();

  // stuff for selection levels
  virtual void DoCuts();

  // stuff for histograms
  virtual void BookHistograms(THistos&);
  virtual void FillHistograms(THistos&);

  // stuff for pdf
  virtual void BookPdfHistos(Int_t i_syst);
  virtual void PdfValues();

  // stuf for likelihoods
  virtual void BookLikeHistos();

  // stuff for calculations

  Int_t event_counter;
  Int_t divideSample;
  
  Int_t leptonSep;
  Int_t lepSample;
  Int_t pileup;
  Int_t CutTriggerEle;
  Int_t CutTriggerMuo;
  Int_t CutTriggerEleMuo;

  Int_t alldata_ele;
  Int_t alldata_muo;

  Int_t ss_lep;
  Int_t ss_jets;
  Int_t ss_bjets;
  Int_t ss_goodbtaggedjet;
  Int_t ss_hardntjet;
 
  Double_t WtMass ;
  Double_t ctln ;

  Int_t ngs;

  // #############################################
  // Define Reference Mass Values
  // #############################################
  Double_t mW_UserValue ;   // W mass
  Double_t mt_UserValue ;   // top mass
  Double_t mH_UserValue ;   // Higgs mass

  // #############################################
  // stuff for ttDilepKinFit 20 December 2011
  // Objects are identified with  _ttDKF
  // #############################################
  // Solution Choice
  int ttDKF_SolutionChoice;
  // Choice of Jet Combination Scheme
  int ttDKF_JetCombChoice; 
  // Choice of Jet Combination Scheme
  int ttDKF_HiggsChoice; 
  // Choice of Number of Jets for Permutations
  int ttDKF_njet_UserValue; 
  int ttDKF_njets; 
  // top quark 1
  std::vector<TLorentzVectorWFlags> b1_ttDKF;
  std::vector<TLorentzVectorWFlags> l1_ttDKF;
  std::vector<TLorentzVectorWFlags> n1_ttDKF;
  std::vector<TLorentzVectorWFlags> W1_ttDKF;
  std::vector<TLorentzVectorWFlags> t1_ttDKF;
  // top quark 2
  std::vector<TLorentzVectorWFlags> b2_ttDKF;
  std::vector<TLorentzVectorWFlags> l2_ttDKF;
  std::vector<TLorentzVectorWFlags> n2_ttDKF;
  std::vector<TLorentzVectorWFlags> W2_ttDKF;
  std::vector<TLorentzVectorWFlags> t2_ttDKF;
  // ttbar system
  std::vector<TLorentzVectorWFlags> ttbar_ttDKF;
  // Higgs
  std::vector<TLorentzVectorWFlags> b1_Higgs_ttDKF;
  std::vector<TLorentzVectorWFlags> b2_Higgs_ttDKF;
  std::vector<TLorentzVectorWFlags> Higgs_ttDKF;
  std::vector<double> mHiggsJet1_ttDKF;
  std::vector<double> mHiggsJet2_ttDKF;

  // ======================================================================================
  // Probabilities for reconstruction ....... samor 16.Dec.2012 ...........................
  // ======================================================================================
  std::vector<double> ProbHiggs_ttDKF;
  std::vector<double> ProbTTbar_ttDKF;
  std::vector<double> ProbTotal_ttDKF;


  // -----------------------------
  // Define pdf for Kinematcal Fit
  // -----------------------------
  int myNpdfKinFit;
  TH1D *pdfKinFit[100];
        double *pdfKinFitVec[100];
        int NBins[100];
  double LowerEdge[100], UpperEdge[100], Scale[100];

        TFile * fTruth;
  // -----------------------------


 // #############################################

        // #####################################################################
        // Miguel Fiolhais 10 Nov 2010  
        // #####################################################################
        int NbtagJet;
        std::vector<TLorentzVectorWFlags> MyGoodJetVec;
        std::vector<TLorentzVectorWFlags> MyGoodBtaggedJetVec;
        std::vector<TLorentzVectorWFlags> MyGoodNonBtaggedJetVec;
        std::vector<TLorentzVectorWFlags> MyChoiceJetVec;

	// ######################################################################################
	// introduced by Miki + AO, Ocotber 17
	// ######################################################################################

	// measured variables
	double Ht;
	double Hz; //samor

	double Ml1b1;
	double Ml2b2;
	double Ml1b2;
	double Ml2b1;

	double Ml1b1_truth;
	double Ml2b2_truth;
	double Ml1b2_truth;
	double Ml2b1_truth;

        double pnupanu[4];

        // ======================================================================================
	// truth variables.......................................................................
        // ======================================================================================
        //_____objects_____
	TLorentzVector 	TruthT,    TruthB,    TruthWp, TruthLepP,    TruthNeu, 
			TruthTbar, TruthBbar, TruthWn, TruthLepN,    TruthNeubar; 
        //_____Boost to top_____
	TLorentzVector 	TruthB_BoostedtoT,       TruthWp_BoostedtoT,    TruthLepP_BoostedtoT,       TruthNeu_BoostedtoT, 
			TruthBbar_BoostedtoTbar, TruthWn_BoostedtoTbar, TruthLepN_BoostedtoTbar,    TruthNeubar_BoostedtoTbar; 
        //_____Boost to ttbar_____
	TLorentzVector  TruthT_Boostedtottbar, TruthTbar_Boostedtottbar;
        //_____angles______
	double XX, YY, TruthHasSolution, 
			TruthCos_LepP_T_BoostedtoT,       TruthCos_Neu_T_BoostedtoT,          TruthCos_B_T_BoostedtoT,
			TruthCos_LepN_Tbar_BoostedtoTbar, TruthCos_Neubar_Tbar_BoostedtoTbar, TruthCos_Bbar_Tbar_BoostedtoTbar,
			TruthCos_LepP_B_BoostedtoWp, TruthCos_LepN_Bbar_BoostedtoWn;
        //_____Boost to W+-_____
	TLorentzVector 	TruthB_BoostedtoWp,    TruthLepP_BoostedtoWp,       TruthNeu_BoostedtoWp, 
			TruthBbar_BoostedtoWn, TruthLepN_BoostedtoWn,    TruthNeubar_BoostedtoWn; 

	double Truth_WCosTheta1, Truth_WCosTheta2;
        
        // ======================================================================================
	// reconstructed variables...............................................................
        // ======================================================================================
        double HasSolution, NumberGoodSolutions, NumberGoodEvents, CosTheta1, CosTheta2;
	TLorentzVector Lepton1, Lepton2, B1, B2, Neutrino, Antineutrino, Top1, Top2;

        //_____objects_____
	TLorentzVector 	RecT,    RecB,    RecWp, RecLepP,    RecNeu, 
			RecTbar, RecBbar, RecWn, RecLepN,    RecNeubar,
			RecTTbar, RecHiggs, RecHiggsB1, RecHiggsB2;
        //_____Boost to top_____
	TLorentzVector 	RecB_BoostedtoT,       RecWp_BoostedtoT,    RecLepP_BoostedtoT,       RecNeu_BoostedtoT, 
			RecBbar_BoostedtoTbar, RecWn_BoostedtoTbar, RecLepN_BoostedtoTbar,    RecNeubar_BoostedtoTbar; 
        //_____Boost to ttbar_____
	TLorentzVector  RecT_Boostedtottbar, RecTbar_Boostedtottbar;
        //_____angles______
	double 		RecCos_LepP_T_BoostedtoT,       RecCos_Neu_T_BoostedtoT,          RecCos_B_T_BoostedtoT,
			RecCos_LepN_Tbar_BoostedtoTbar, RecCos_Neubar_Tbar_BoostedtoTbar, RecCos_Bbar_Tbar_BoostedtoTbar,
			RecCos_LepP_B_BoostedtoWp, RecCos_LepN_Bbar_BoostedtoWn;
        //_____Boost to W+-_____
	TLorentzVector 	RecB_BoostedtoWp,    RecLepP_BoostedtoWp,       RecNeu_BoostedtoWp, 
			RecBbar_BoostedtoWn, RecLepN_BoostedtoWn,    RecNeubar_BoostedtoWn; 

	double WCosTheta1, WCosTheta2;

        //_____Higgs Masses______
	double 		RecMassHiggsJet1,	RecMassHiggsJet2, 	RecProbTotal_ttH;


        // ======================================================================================
        // Truth versus Rec variables............................................................
        // ======================================================================================
        double Truth_myCosMlb_lp, Truth_myCosMlb_ln,
                 Rec_myCosMlb_lp,   Rec_myCosMlb_ln;
        double Truth_mylbMass_lp, Truth_mylbMass_ln,
                 Rec_mylbMass_lp,   Rec_mylbMass_ln;

        // ======================================================================================
        // Truth Leptons Reconstructed Correctly.................................................
        // ======================================================================================
        Int_t  Dilep_Truth_Reco_Match;
        Int_t  myTruth_WtauDecay;
        std::vector<TLorentzVectorWFlags> myTruthLepVec;
        std::vector<TLorentzVectorWFlags> myTruthNeuVec;
        std::vector<TLorentzVectorWFlags> myTruth_t_Vec;
        std::vector<TLorentzVectorWFlags> myTruth_W_Vec;
        std::vector<TLorentzVectorWFlags> myTruth_b_Vec;

        // ======================================================================================
        // Matched Truth with Reconstructed Leptons, jets and missing for Resolution Studies.....
        // ======================================================================================
	// Truth information
        std::vector<TLorentzVectorWFlags> myResolution_Truth_Lep;
        std::vector<TLorentzVectorWFlags> myResolution_Truth_Jet;
        std::vector<TLorentzVectorWFlags> myResolution_Truth_Mis;
	// Reconstructed information
        std::vector<TLorentzVectorWFlags> myResolution_Reco_Lep;
        std::vector<TLorentzVectorWFlags> myResolution_Reco_Jet;
        std::vector<TLorentzVectorWFlags> myResolution_Reco_Mis;

	// Random variables
	TRandom3 rnd;

	// ######################################################################################

  double root1;

  Double_t corrlum;
  Double_t corrnorm_wjets;
  Double_t corrnorm_zjets;
  Double_t corrnorm_tt;
  Double_t corrnorm_wz;
  Double_t corrnorm_db;
  Double_t corrnorm_st;

  Double_t ptcut;
  Double_t cutwmass;
  Double_t cuttopmass;
  Double_t cutptjet;
  Double_t deltarmincut;
  Double_t cutsqrtchi2; 
  Double_t cutchi2;

  Int_t    nlepcut;
  Int_t    njetmin;

  Int_t ntagged;
  double ref_weight;


	// ######################################################################################
	// introduced by AO
	// ######################################################################################

        // =========================================================
        // ===== AO 8 Oct 2010 ===================== below =========
        // =========================================================
        //           Get MY OWN variables out of the box
        // =========================================================
        // ===== AO 8 Oct 2010 ===================== below =========
        // =========================================================

        Int_t myRunNumber;       //...run number

        //.................t
        Int_t my_ITQ;            //...line for top quark
        Int_t my_IQ1;            //...counter for top quark
        TLorentzVector my_t;
        //.................tbar
        Int_t my_ITB;            //...line for anti top quark
        Int_t my_IQ2;            //...counter for anti top quark
        TLorentzVector my_tb;


        //.................W+
        Int_t my_IWP;            //...line for W+
        Int_t my_IW1;            //...counter for W+
        TLorentzVector my_Wp;
        //.................W-
        Int_t my_IWN;            //...line for W-
        Int_t my_IW2;            //...counter for W-
        TLorentzVector my_Wn;


        //.................b
        Int_t my_IBQ;            //...line for b
        Int_t my_IB1;            //...counter for bb
        TLorentzVector my_b;
        //.................bb
        Int_t my_IBB;            //...line for bb
        Int_t my_IB2;            //...counter for bb
        TLorentzVector my_bb;


        //.................s
        Int_t my_ISQ;            //...line for s
        Int_t my_IS1;            //...counter for s
        TLorentzVector my_s;
        //.................sb
        Int_t my_ISB;            //...line for sb
        Int_t my_IS2;            //...counter for sb
        TLorentzVector my_sb;


        //.................d
        Int_t my_IDWQ;            //...line for d
        Int_t my_IDW1;            //...counter for d
        TLorentzVector my_dw;
        //.................db
        Int_t my_IDWB;            //...line for db
        Int_t my_IDW2;            //...counter for db
        TLorentzVector my_dwb;


        //.................W+->f1f2
        Int_t my_IWPf1;          //...line for f1
       	Int_t my_IWPf2;          //...line for f2
        Int_t my_IWPf1_Coun;     //...counter for f1
        Int_t my_IWPf2_Coun;     //...counter for f2
        Int_t my_IWPtau_Neu;     //...counter for Tau Neutrinos
        Int_t my_IWPtau_elNeu;   //...counter for ele Neutrinos (from tau decay)
        Int_t my_IWPtau_muNeu;   //...counter for muo Neutrino  (from tau decay)
        Int_t my_pdgID_Wp_dw;    //...Code of 1st W+ Daughter
        Int_t my_pdgID_Wp_up;    //...Code of 2nd W+ Daughter
        TLorentzVector my_Wpf1;
        TLorentzVector my_Wpf2;
        //.................W-->f1f2
        Int_t my_IWNf1;          //...line for f1
        Int_t my_IWNf2;          //...line for f2
        Int_t my_IWNf1_Coun;     //...counter for f1
        Int_t my_IWNf2_Coun;     //...counter for f2
        Int_t my_IWNtau_Neu;     //...counter for Tau Neutrinos
        Int_t my_IWNtau_elNeu;   //...counter for ele Neutrinos (from tau decay)
        Int_t my_IWNtau_muNeu;   //...counter for muo Neutrino  (from tau decay)
        Int_t my_pdgID_Wn_dw;    //...Code of 1st W- Daughter
        Int_t my_pdgID_Wn_up;    //...Code of 2nd W- Daughter
        TLorentzVector my_Wnf1;
        TLorentzVector my_Wnf2;
        // =========================================================
        // ===== AO 8 Oct 2010 ===================== above =========
        // =========================================================
        Int_t myEVE_Total;
        //....semileptonic....
        Int_t myEVE_semi_ele;
        Int_t myEVE_semi_mu;
        Int_t myEVE_semi_nontau_ele;
        Int_t myEVE_semi_nontau_mu;
        Int_t myEVE_semi_tau;
        Int_t myEVE_semi_leptau;
        Int_t myEVE_semi_eletau;
        Int_t myEVE_semi_mutau;
        //....dileptonic............with taus....
        Int_t myEVE_di_eleele;
        Int_t myEVE_di_mumu;
        Int_t myEVE_di_elemu;
        //....dileptonic............without taus....
        Int_t myEVE_di_nontau_eleele;
        Int_t myEVE_di_nontau_elemu;
        Int_t myEVE_di_nontau_eletau;
        Int_t myEVE_di_nontau_mumu;
        Int_t myEVE_di_nontau_mutau;
        Int_t myEVE_di_nontau_tautau;
        //....multileptons............
        Int_t myEVE_0leptons;
        Int_t myEVE_3leptons;
        Int_t myEVE_4leptons;
        Int_t myEVE_more4leptons;
        //...top quarks
	Int_t myEVE_ttbar;
        Int_t myEVE_4tops;
        Int_t myEVE_nonttbar;
        //.....................discriminate what ttbar events do we have.....
        Int_t myttbar_00;
        Int_t myttbar_01;
        Int_t myttbar_02;
        Int_t myttbar_03;
        Int_t myttbar_10;
        Int_t myttbar_11;
        Int_t myttbar_12;
        Int_t myttbar_13;
        Int_t myttbar_20;
        Int_t myttbar_21;
        Int_t myttbar_22;
        Int_t myttbar_23;
        Int_t myttbar_30;
        Int_t myttbar_31;
        Int_t myttbar_32;
        Int_t myttbar_33;
        //...W bosons
	Int_t myEVE_WW;
        Int_t myEVE_4W;
        Int_t myEVE_nonWW;
        //.....................discriminate what bbbar events do we have.....
        Int_t mybbbar_00;
        Int_t mybbbar_01;
        Int_t mybbbar_02;
        Int_t mybbbar_03;
        Int_t mybbbar_10;
        Int_t mybbbar_11;
        Int_t mybbbar_12;
        Int_t mybbbar_13;
        Int_t mybbbar_20;
        Int_t mybbbar_21;
        Int_t mybbbar_22;
        Int_t mybbbar_23;
        Int_t mybbbar_30;
        Int_t mybbbar_31;
        Int_t mybbbar_32;
        Int_t mybbbar_33;
        //...b quarks
        Int_t myEVE_bb;
        Int_t myEVE_4b;
        Int_t myEVE_nonbb;
        //.....................discriminate what ssbar events do we have.....
        Int_t myssbar_00;
        Int_t myssbar_01;
        Int_t myssbar_02;
        Int_t myssbar_03;
        Int_t myssbar_10;
        Int_t myssbar_11;
        Int_t myssbar_12;
        Int_t myssbar_13;
        Int_t myssbar_20;
        Int_t myssbar_21;
        Int_t myssbar_22;
        Int_t myssbar_23;
        Int_t myssbar_30;
        Int_t myssbar_31;
        Int_t myssbar_32;
        Int_t myssbar_33;
        //...s quarks
        Int_t myEVE_ss;
        Int_t myEVE_4s;
        Int_t myEVE_nonss;
        //.....................discriminate what ddbar events do we have.....
        Int_t myddbar_00;
        Int_t myddbar_01;
        Int_t myddbar_02;
        Int_t myddbar_03;
        Int_t myddbar_10;
        Int_t myddbar_11;
        Int_t myddbar_12;
        Int_t myddbar_13;
        Int_t myddbar_20;
        Int_t myddbar_21;
        Int_t myddbar_22;
        Int_t myddbar_23;
        Int_t myddbar_30;
        Int_t myddbar_31;
        Int_t myddbar_32;
        Int_t myddbar_33;
        //...d quarks
        Int_t myEVE_dd;
        Int_t myEVE_4d;
        Int_t myEVE_nondd;
        // =========================================================
        // ===== AO 8 Oct 2010 ===================== above =========
        // =========================================================
	// ######################################################################################

  Int_t ntruthlep;
  Int_t ntruthele;
  Int_t ntruthmu;
  Int_t ntrutheletau;
  Int_t ntruthmutau;
  Int_t ntruthtau;
  Int_t ntruthleptau;

  Double_t NonbJet1_weight;
  Double_t NonbJet2_weight;
  Double_t Leptonicb_weight;
  Double_t Hadronicb_weight;

  Double_t jet_weight[30];
  Int_t jet_ntracks[30];

  TLorentzVector ll;
  TLorentzVector llmiss;

  Int_t           PassedEF_mu20;
  Int_t           PassedEF_e20_loose;
  Int_t           PassedEF_mu10; //4 fev 2010 samor
  Int_t           PassedEF_e10_medium;

  Int_t           Passed_L1_EM2;//24.05.10 samor
  Int_t           Passed_L1_MU0;

  // AO 05.nov.2012 ================================================
  // New systematic error stuff
  double LumSys;
  double LumSysZj;
  double LumSysWj;
  double LumSysWZp;
  double LumSystt;
  double LumSysst;
  double SysPtLep;
  double SysPtJet;
  double SysPtPho;
  double SysPtMiss;
  double SysEtot;
  double SysPtGluon;
  double NBinFactor;
  // AO 05.nov.2012 ================================================

  ofstream of;
  ofstream f2;

  void Calculations();
  void Calculations2();
  void Calculations3();
  void ttDilepKinFit();

  //stuff for after events loop
  virtual void PostLoopCalculations();

  TH2D *histo2D[50];

  // New stuff for duplicated events
  bool CheckDoubleEvents(int run, int event);

  // New stuff for Truth Matching
  bool isTruthMatched_Electron( int el_idx );
  bool isTruthMatched_Muon( int mu_idx );

  // Returns Delta R*R .
  virtual double DeltaR2( double eta1, double eta2, double phi1, double phi2);

};

#endif
