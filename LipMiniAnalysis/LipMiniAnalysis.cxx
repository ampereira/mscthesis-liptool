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
//   Note: 	this code was built from LipCbrAnalysis developed by 
//     	   	Nuno Castro   (nfcastro@lipc.fis.uc.pt)  and 
//         	Filipe Veloso (fveloso@lipc.fis.uc.pt)             
//
//
// #################################################################################
//
//  file: LipMiniAnalysis.cxx
//
//  purpose: this is the code of the LipMiniAnalysis class, which is a skeleton
//  for particle physics analyses using ROOT
//
//  authors: nfcastro, fveloso
//  first version: nov.2005
//
//  last change: 23.may.2011
//  by: fveloso
//
// #################################################################################

#define LipMiniAnalysis_cxx

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <vector>
#include <TRandom.h>
#include <TROOT.h>
#include <TLorentzVector.h>
#include <iomanip>
#include <fstream>

#include <cmath>

#include <TChain.h>
#include <TFile.h>

using std::cout ;
using std::endl ;
using std::cerr ;
using std::string ;
using std::stringstream ;

#include "LipMiniAnalysis.h"

// #############################################################################
bool LorentzVecComp(TLorentzVectorWFlags a, TLorentzVectorWFlags b) {
// #############################################################################
//
//  purpose: this function is used to sort the particles by decreasing pT
//
//  authors: fveloso
//  first version: 27.fev.2007
//
// #############################################################################

   return a.Pt() > b.Pt();
}

// #############################################################################
LipMiniAnalysis::LipMiniAnalysis(){
// #############################################################################
//
//  purpose: this is the constructor of the LipMiniAnalysis class
//
//  authors: ROOT, fveloso
//  first version: ??.???.??
//
//  last change: 30.nov.2009
//  by: nfcastro
//
// #############################################################################

//  cout = ostream("LipMiniAnalysis.out");

  LipMiniAnalysisVersion = 1.00;

  // stuff for cuts
  MaxCuts = 1;

}

// #############################################################################
void LipMiniAnalysis::DefaultValues() {
// #############################################################################
//
//  purpose: 
//
//  authors: fveloso
//  first version: 26.fev.2007
//
//  last change: 19.Nov.2012
//  by: A.Onofre
//
// #############################################################################

  // units
  if (!nTuple) GeV = 1;
  else GeV = nTuple->GeV();
  keV = GeV*.000001;
  MeV = GeV*.001;
  TeV = GeV*1000.;

  // general stuff
  OutputFileName = "../Output/summary_LipMiniAnalysis.root";
  BackgroundPdfFile = "background_pdf.root";
  GRLXMLFile = "grl.xml";
  // SignalPdfFiles are defined in DefineSignalPdfSamples()
  Luminosity = 10.;
  Sample = 0;
  subperiod = -1;

  // masses
  mZ = 91.19*GeV;
  mW = 80.40*GeV;
  mt = 172.5*GeV;

  // ---------------------------------------------------------
  // NOTE: Here only parameter initializations are performed
  //       To set real values set them in the analysis code 
  //           void ttH_dilep::UserValues(){}  or equivalent
  // ---------------------------------------------------------
  // cuts for electrons
  EtaCutEle 	=  999.;
  PtCutEle  	= -999.*GeV;
  IsoCutEle 	=  999.;
  RemoveElCrack	= -1   ;

  // cuts for muons
  EtaCutMuo 	=  999.;
  PtCutMuo  	= -999.*GeV;
  IsoCutMuo 	=  999.;

  // cuts for photons
  EtaCutPho 	=  999.;
  PtCutPho 	= -999.*GeV;
  IsoCutPho 	=  999.;

  // cuts for jets
  EtaCutJet 	=  999.;
  PtCutJet 	= -999.*GeV;
  IsoCutJet 	=  999.;

  // jet energy calibration
  LightJetCalib = 1.;
  BJetCalib 	= 1.;

  // btag cut
  BTagCut = 999.; 

  // overlap of jets with photons and taus
  JetOverlap = 1;

  // Options for different reconstruction of objets within the same ntuples
  RecoType = 0;

  // discriminant analysis
  DoLike 		=    0;
  LogELikeMinCut 	= -999;
  LogELikeMaxCut 	=  999;
  Log10LikeMinCut 	= -999;
  Log10LikeMaxCut 	=  999;

  // flag for data events
  isData = 0;

  // flag to do good objects selection
  doGoodObjSelection = 0;

  // flags for data events
  isQCDantie = 0;
  isQCDmmm = 0;
  isPretag = 0;

  data_e_periodC = 0;
  data_e_periodD = 0;
  data_e_periodE = 0;
  data_e_periodF = 0;
  data_e_periodI = 0;

  data_m_periodB = 0;
  data_m_periodC = 0;
  data_m_periodD = 0;
  data_m_periodE = 0;
  data_m_periodF = 0;
  data_m_periodI = 0;

  data_e_periodA_F = 0;
  data_m_periodA_F = 0;

  data_e_periodG_H = 0;
  data_m_periodG_H = 0;

}

// #############################################################################
void LipMiniAnalysis::UserValues() {
// #############################################################################
//
//  purpose: 
//
//  notes: should be replaced by the user!
//
//  authors: fveloso
//  first version: 26.fev.2007
//
// #############################################################################

}

// #############################################################################
void LipMiniAnalysis::CommandLineOptions(int print) {
// #############################################################################
//
//  purpose: process command line options
// 	     Note: in this version cuts fro electrons separated from muons
//
//  notes: for help see for example
//    http://www.gnu.org/software/libc/manual/html_node/Getopt-Long-Options.html
//
//  authors: fveloso
//  first version: 21.jun.2007
//
//  last change: 23.Nov.2012
//  by: A.Onofre
//
// #############################################################################

  // define command line options
  static struct option long_options[] = {
    // user options
    {"User",                	required_argument, 0, 'a'}, // char[]

    // general stuff
    {"OutputFileName",      	required_argument, 0, 'b'}, // char[]
    {"BackgroundPdfFile",   	required_argument, 0, 'c'}, // char[]
    {"Luminosity",          	required_argument, 0, 'd'},
    {"Sample",              	required_argument, 0, 'e'}, // int

    // masses
    {"mZ",                  	required_argument, 0, 'f'}, // * GeV
    {"mW",                  	required_argument, 0, 'g'}, // * GeV
    {"mt",                  	required_argument, 0, 'h'}, // * GeV

    // cuts for electrons
    {"EtaCutEle",           	required_argument, 0, 'i'},
    {"PtCutEle",            	required_argument, 0, 'j'}, // * GeV
    {"IsoCutEle",           	required_argument, 0, 'k'},
    {"RemoveElCrack",       	required_argument, 0, 'l'}, // int

    // cuts for muons
    {"EtaCutMuo",           	required_argument, 0, '1'},
    {"PtCutMuo",            	required_argument, 0, '2'}, // * GeV
    {"IsoCutMuo",           	required_argument, 0, '3'},

    // cuts for photons
    {"EtaCutPho",           	required_argument, 0, 'm'},
    {"PtCutPho",            	required_argument, 0, 'n'}, // * GeV
    {"IsoCutPho",           	required_argument, 0, 'o'},

    // cuts for jets
    {"EtaCutJet",           	required_argument, 0, 'p'},
    {"PtCutJet",            	required_argument, 0, 'q'}, // * GeV
    {"IsoCutJet",           	required_argument, 0, 'r'},

    // jet energy calibration? and btag cut; both for fullsim topview >= 12-13
    {"LightJetCalib",       	required_argument, 0, 's'},
    {"BJetCalib",           	required_argument, 0, 't'},
    {"BTagCut",             	required_argument, 0, 'u'},

    // discriminant analysis
    {"DoLike",              	required_argument, 0, 'v'}, // int
    {"LogELikeMinCut",      	required_argument, 0, 'w'},
    {"LogELikeMaxCut",      	required_argument, 0, 'x'},
    {"Log10LikeMinCut",     	required_argument, 0, 'y'},
    {"Log10LikeMaxCut",     	required_argument, 0, 'z'},

    // overlap of jets with photons and taus
    {"JetOverlap",          	required_argument, 0, 'A'}, // int

    {"GRLXMLFile",          	required_argument, 0, 'B'}, // char[]

    {"isData",              	required_argument, 0, 'C'}, // int
    {"doGoodObjSelection",  	required_argument, 0, 'D'}, // int

    {"isQCDantie",          	required_argument, 0, 'E'}, // int
    {"isQCDmmm",      	    	required_argument, 0, 'F'}, // int
    {"data_e_periodC",      	required_argument, 0, 'G'}, // int
    {"data_e_periodD",      	required_argument, 0, 'H'}, // int
    {"data_e_periodE",      	required_argument, 0, 'I'}, // int
    {"data_e_periodF",      	required_argument, 0, 'J'}, // int

    {"isPretag",            	required_argument, 0, 'K'}, // int
    {"data_m_periodB",      	required_argument, 0, 'L'}, // int
    {"data_m_periodC",      	required_argument, 0, 'M'}, // int
    {"data_m_periodD",      	required_argument, 0, 'N'}, // int
    {"data_m_periodE",      	required_argument, 0, 'O'}, // int
    {"data_m_periodF",      	required_argument, 0, 'P'}, // int
    {"subperiod",           	required_argument, 0, 'Q'}, // int
   

    {"data_e_periodA_F",    	required_argument, 0, 'R'}, // int
    {"data_e_periodG_H",    	required_argument, 0, 'S'}, // int
    {"data_m_periodA_F",    	required_argument, 0, 'T'}, // int
    {"data_m_periodG_H",    	required_argument, 0, 'U'}, // int

    {"data_e_periodI",      	required_argument, 0, 'V'}, // int
    {"SetSystematicsFileName",  required_argument, 0, 'W'}, // char[]
    {"data_m_periodI",      	required_argument, 0, 'X'}, // int

    // Options for different reconstruction of objets within the same ntuples
    {"RecoType",            	required_argument, 0, 'Z'}, // int

    // Terminate the array with an element containing all zeros.
    {0, 0, 0, 0}
  };

  // process command line options
  optind = 1; // trick to process several times with getopt or getopt_long
  int c;
  while (1) {
    int o = 0;
    int option_index = 0;
    c = getopt_long (argc, argv, "a:b:c:d:e:f:g:h:i:j:k:l:1:2:3:m:n:o:p:q:r:s:t:u:v:w:x:y:z:A:B:C:D:E:F:G:H:I:J:K:L:M:N:O:P", long_options, &option_index);


    // Detect the end of the options.
    if (c == -1) break;

    // process each option
    switch (c) {
      case 'a':
        UserCommandLineOptions(optarg,print);
        break;
      case 'b':
        OutputFileName=optarg;
        o = 1;
        break;
      case 'c':
        BackgroundPdfFile = optarg;
        o = 1;
        break;
      case 'd':
        Luminosity = atof(optarg);
        o = 1;
        break;
      case 'e':
        Sample = atoi(optarg);
        o = 1;
        break;
      case 'f':
        mZ = atof(optarg)*GeV;
        o = 1;
        break;
      case 'g':
        mW = atof(optarg)*GeV;
        o = 1;
        break;
      case 'h':
        mt = atof(optarg)*GeV;
        o = 1;
        break;
      case 'i':
        EtaCutEle = atof(optarg);
        o = 1;
        break;
      case 'j':
        PtCutEle = atof(optarg)*GeV;
        o = 1;
        break;
      case 'k':
        IsoCutEle = atof(optarg);
        o = 1;
        break;
      case 'l':
        RemoveElCrack = atoi(optarg);
        o = 1;
        break;
      case '1':
        EtaCutMuo = atof(optarg);
        o = 1;
        break;
      case '2':
        PtCutMuo = atof(optarg)*GeV;
        o = 1;
        break;
      case '3':
        IsoCutMuo = atof(optarg);
        o = 1;
        break;
      case 'm':
        EtaCutPho = atof(optarg);
        o = 1;
        break;
      case 'n':
        PtCutPho = atof(optarg)*GeV;
        o = 1;
        break;
      case 'o':
        IsoCutPho = atof(optarg);
        o = 1;
        break;
      case 'p':
        EtaCutJet = atof(optarg);
        o = 1;
        break;
      case 'q':
        PtCutJet = atof(optarg)*GeV;
        o = 1;
        break;
      case 'r':
        IsoCutJet = atof(optarg);
        o = 1;
        break;
      case 's':
        LightJetCalib = atof(optarg);
        o = 1;
        break;
      case 't':
        BJetCalib = atof(optarg);
        o = 1;
        break;
      case 'u':
        BTagCut = atof(optarg);
        o = 1;
        break;
      case 'v':
        DoLike = atoi(optarg);
        o = 1;
        break;
      case 'w':
        LogELikeMinCut = atof(optarg);
        o = 1;
        break;
      case 'x':
        LogELikeMaxCut = atof(optarg);
        o = 1;
        break;
      case 'y':
        Log10LikeMinCut = atof(optarg);
        o = 1;
        break;
      case 'z':
        Log10LikeMaxCut = atof(optarg);
        o = 1;
        break;
      case 'A':
        JetOverlap = atoi(optarg);
        o = 1;
        break;
      case 'B':
        GRLXMLFile = optarg;
        o = 1;
        break;
      case 'C':
        isData = atoi(optarg);
        o = 1;
        break;
      case 'D':
        doGoodObjSelection = atoi(optarg);
        o = 1;
        break;
      case 'E':
        isQCDantie = atoi(optarg);
        o = 1;
        break;
      case 'F':
        isQCDmmm = atoi(optarg);
        o = 1;
        break;
      case 'G':
        data_e_periodC = atoi(optarg);
        o = 1;
        break;
      case 'H':
        data_e_periodD = atoi(optarg);
        o = 1;
        break;
      case 'I':
        data_e_periodE = atoi(optarg);
        o = 1;
        break;
      case 'J':
        data_e_periodF = atoi(optarg);
        o = 1;
        break;
      case 'K':
        isPretag = atoi(optarg);
        o = 1;
        break;
      case 'L':
        data_m_periodB = atoi(optarg);
        o = 1;
        break;
      case 'M':
        data_m_periodC = atoi(optarg);
        o = 1;
        break;
      case 'N':
        data_m_periodD = atoi(optarg);
        o = 1;
        break;
      case 'O':
        data_m_periodE = atoi(optarg);
        o = 1;
        break;
      case 'P':
        data_m_periodF = atoi(optarg);
        o = 1;
        break;
      case 'Q':
        subperiod = atoi(optarg);
        o = 1;
        break;
      case 'R':
         data_e_periodA_F = atoi(optarg);
        o = 1;
        break;
      case 'S':
         data_e_periodG_H = atoi(optarg);
        o = 1;
        break;
      case 'T':
         data_m_periodA_F = atoi(optarg);
        o = 1;
        break;
      case 'U':
         data_m_periodG_H = atoi(optarg);
        o = 1;
        break;
      case 'V':
         data_e_periodI = atoi(optarg);
        o = 1;
        break;
      case 'W':
        SetSystematicsFileName = optarg;
        o = 1;
        break;
      case 'X':
         data_m_periodI = atoi(optarg);
        o = 1;
        break;
      case 'Z':
        RecoType = atoi(optarg);
        o = 1;
        break;
      case '?':
        /* getopt_long already printed an error message. */
        break;
      default:
        abort ();
    }
    
    if (o==1 && print==1)
     cout << "  " << long_options[option_index].name <<" = " <<  optarg << endl;
  }
}

// #############################################################################
void LipMiniAnalysis::UserCommandLineOptions(char *options, int print) {
// #############################################################################
//
//  purpose: 
//
//  notes: should be replaced by the user!
//
//  authors: fveloso
//  first version: 24.jun.2007
//
// #############################################################################

}

// #############################################################################
void LipMiniAnalysis::SetSystematics() {
// #############################################################################

  Syst.clear();
  SystName.clear();

  if (SetSystematicsFileName!="") {
    ifstream f;
    f.open(SetSystematicsFileName.c_str());
    if (!f) {
      cerr << "  Error: it wasnt possible to open " << SetSystematicsFileName << " which sets Output Directories and Normal or Systematic studies" << endl;   
      exit(1);
    }
    // reads file line by line
    string Linha;
    Int_t recotype;
    string nomesys;
    while (getline(f, Linha)) {
      // process lines
      stringstream sLinha;
      sLinha << Linha;
      sLinha >> recotype >> nomesys;
      Syst.push_back(recotype);
      SystName.push_back(nomesys);
    }
    if (Syst.size()==0) {
      cerr << "  Error: it wasnt possible to open " << SetSystematicsFileName << " which sets Output Directories and Normal or Systematic studies" << endl;   
      exit(1);      
    }
  } else {  
    Syst.push_back(0); SystName.push_back("../Output/");
  }

}

// #############################################################################
void LipMiniAnalysis::Start(int i_argc, char *const *i_argv) {
// #############################################################################
//
//  purpose: 
//
//  authors: fveloso
//  first version: 16.nov.2006
//
//  last change: 19.Nov.2012
//  by: A.Onofre
//
// #############################################################################

  argc = i_argc;
  argv = i_argv; 

  // Set variables' values
  DefaultValues();
  UserValues();
  CommandLineOptions(0);
  SetSystematics();

  // Redirect standard output
  char LogName[1000];
  char ErrName[1000];
  sprintf(LogName, "%s%s.out", SystName[0].c_str(), OutputFileName.c_str());
  sprintf(ErrName, "%s%s.err", SystName[0].c_str(), OutputFileName.c_str());

  // .....comment the two lines below if you would like to have the output writen to screen
  freopen(LogName, "w", stdout);
  freopen(ErrName, "w", stderr);

  // Prints logo
  cout << endl;
  cout << " #################################################################################" << endl;
  cout << "										 " << endl;
  cout << "    _      _         __  __ _       _                        _           _    " << endl; 
  cout << "   | |    (_)       |  \\/  (_)     (_)     /\\               | |         (_)   " << endl; 
  cout << "   | |     _ _ __   | \\  / |_ _ __  _     /  \\   _ __   __ _| |_   _ ___ _ ___ " << endl;
  cout << "   | |    | | '_ \\  | |\\/| | | '_ \\| |   / /\\ \\ | '_ \\ / _` | | | | / __| / __|" << endl;
  cout << "   | |____| | |_) | | |  | | | | | | |  / ____ \\| | | | (_| | | |_| \\__ \\ \\__ \\" << endl;
  cout << "   |______|_| .__/  |_|  |_|_|_| |_|_| /_/    \\_\\_| |_|\\__,_|_|\\__, |___/_|___/" << endl;
  cout << "            | |                                                 __/ |          " << endl;
  cout << "            |_|                                                |___/           " << endl;
  cout << "                                                                              " << endl;
  cout << " #################################################################################" << endl;
  cout << "                                                                              " << endl;
  cout << "   LipMiniAnalysis " << std::dec << std::fixed << std::right << std::setw(5)
    << std::setprecision(2) << LipMiniAnalysisVersion << endl;
  cout << "                                                                              " << endl;
  cout << "                                              \\      | | / |  /    |   |  ||  " << endl;
  cout << "                                               '-___-+-+'  |_/     |   |  |   " << endl;
  cout << "   By Antonio Onofre                          `--___-++--^'|       |   | /|   " << endl;
  cout << "      (antonio.onofre@cern.ch)                       ||    |       |   |' |   " << endl;
  cout << "   Date: 23.Nov.2012                          --___  ||    |       | _/|  |   " << endl;
  cout << "                                              ==---:-++___ |  ___--+'  |  |   " << endl;
  cout << "                                                '--_'|'---'+--___  |   |  |   " << endl;
  cout << "                                                    '+-_   |     '-+__ |  |   " << endl;
  cout << "                                              ._.          ._.       ._____   " << endl;
  cout << "                                              | |          | |       | ___ \\  " << endl;
  cout << "                                              | |_.        | |       | .___/  " << endl;
  cout << "                                              |___|        |_|       |_|      " << endl;
  cout << "                                                                              " << endl;
  cout << "   Note: 	this code was built from LipCbrAnalysis developed by             " << endl;
  cout << "     		Nuno Castro   (nfcastro@lipc.fis.uc.pt)  and                     " << endl;
  cout << "         	Filipe Veloso (fveloso@lipc.fis.uc.pt)                           " << endl;
  cout << "                                                                              " << endl;
  cout << "                                                                              " << endl;
  cout << " #################################################################################" << endl;


  // Prints date, time and host name
  struct tm * timeinfo;
  char timestring[50];
  time(&rawtimeI);
  timeinfo = localtime(&rawtimeI);
  strftime(timestring, 100, "%a %d %b %Y %T", timeinfo);
  char hostnamestring[100];
  gethostname(hostnamestring, 100);

  cout << endl;
  cout << "  Date and Time:  " << timestring << endl;
  cout << "  Host name: " << hostnamestring << endl;
  cout << endl;
  cout << " #################################################################################" << endl;
  cout << endl;

  // Look again into the CommandLineOptions but printing this time
  CommandLineOptions(1);

  // Create output file
  CreateOutputFile();

  if (DoLike) ++MaxCuts;

  if (isData==0){
     TMonteCarlo mc_tot(0, -1, 1, 1, "Total MC", MaxCuts);
     MonteCarlo.push_back(mc_tot);
  }else{
     TMonteCarlo data_tot(0, -1, 1, 1, "Total data", MaxCuts);
     MonteCarlo.push_back(data_tot);
  }

  cout << endl;
  cout << " #################################################################################" << endl;
  cout << endl;

  DefineSamples();

  cout << endl;
  cout << " #################################################################################" << endl;
  cout << endl;

  cout << "Will run for the following ref/systematics:" << endl;
  for (Int_t i_syst=0; i_syst<Syst.size(); ++i_syst)
    cout << " " << Syst[i_syst] << " " << SystName[i_syst] << endl;

  cout << endl;
  cout << " #################################################################################" << endl;
  cout << endl;

  if(isData==1){
    cout << " GRL Information Directly from Mini Ntuple " << endl;
  }

  for (Int_t i_syst=0; i_syst<Syst.size(); ++i_syst) {
    // book histograms (just for total MC/data)
    RecoType=Syst[i_syst];
    OutputFile[i_syst]->cd();
    THisto histos;
    MonteCarlo[0].histo.push_back(histos);
    BookHistograms(MonteCarlo[0].histo[i_syst].histo);
    // prepare counters of events (all MC/data samples)
    for (Int_t mc=0; mc<MonteCarlo.size(); ++mc) {
      MonteCarlo[mc].p_nSelEvt.push_back(std::vector<Double_t>());  
      MonteCarlo[mc].p_nSelWeightedEvt.push_back(std::vector<Double_t>());  
      MonteCarlo[mc].p_nSelWeightedEvtErr.push_back(std::vector<Double_t>());  
      for (int i=0; i <= MaxCuts; ++i) {
        MonteCarlo[mc].p_nSelEvt[i_syst].push_back(0.);
        MonteCarlo[mc].p_nSelWeightedEvt[i_syst].push_back(0.);
        MonteCarlo[mc].p_nSelWeightedEvtErr[i_syst].push_back(0.);
      }
    }
  }
  OutputFile[0]->cd();

/* Comentado para performance; descomentar de a análise necessitar de pdf e likelihoods
  // book pdf histos
  TDirectory *dir = gDirectory;
  for (Int_t i_syst=0; i_syst<Syst.size(); ++i_syst) {
    OutputFile[i_syst]->cd();
    std::vector<TH1D*> vec;
    PdfVec.push_back(vec);
    BookPdfHistos(i_syst);
  }
  dir->cd();

  // create ntuple to store surviving events
  CreateOutputNtuple();

  // book likelihood histos
  ThisEventPdfValues.resize(PdfVec[0].size());
  if (DoLike) {
    DefineSignalPdfSamples();
    ReadPdf();

    cout << endl;
    cout << " #################################################################################" << endl;
    cout << endl;
  }

  if (DoLike) {
    BookLikeHistos();
    ProcessBookLikeHistos();
  }
*/

  // dummy canvas to suppress some root warnings
  TCanvas *mydummycanvas=new TCanvas();

  // loop over files
  for (int f=0; f < Input.size(); ++f) {

    // time
    time_t rawtimeD;
    time(&rawtimeD);

    cout << "  Next file: " << Input.Name(f) << endl;
//    cout << "    ntuple is of type " << Input.Type(f) << endl;

    // declare nTuple type
    if (Input.Type(f) == "MiniTTHReader") nTuple = new MiniTTHReader(isData);
    else {
      cout << "  nTuple is of unknown type: " << Input.Type(f) << " " << endl;
      exit(0); 
    }

    // time
    time_t rawtimeA;
    time(&rawtimeA);
    cout << "  Initializing time was " << difftime(rawtimeA, rawtimeD) << " s " << endl;
    cout << endl;

    // opening a file stored in a lustre system is very slow in root if lustre 
    // is overloaded, so check if it is possible to copy the root file to 
    // memory or temporary disc space.

    Int_t UseTempFile=0;
    char localfilename[400];

// comentar daqui para baixo
    TRandom3 r(0);
    long int insize;
    struct statvfs buf;
    char dirtmp[200];
    long int freespace;
    Int_t error;



#ifdef MEASURE_READFILE

#endif




    // open original file and check its size
    std::ifstream infile(Input.Name(f).c_str(), std::ios_base::binary|std::ios::ate);
    insize = (long int) infile.tellg();
    printf("Size of Original File: %ld\n", insize);

    // try to copy th original file to /dev/shm, $TMPDIR or /tmp
    for (Int_t i=0; i<3; ++i) {
      if (i==0) sprintf(dirtmp, "/dev/shm");
      else if (i==1) {
        char * envtmpdir = getenv("TMPDIR");
        // if TMPDIR isn't defined then go to next iteraction (use /tmp)
        if (envtmpdir==NULL) continue;
        sprintf(dirtmp, "%s", envtmpdir);
      } else if (i==2) sprintf(dirtmp, "/tmp");

      // check free space (it also checks for dirtmp existence)
      freespace = -1;
      if ((statvfs(dirtmp, &buf))==0) freespace = buf.f_bfree*buf.f_bsize;
      printf("Free Space in %s: %ld\n", dirtmp, freespace);
      if (freespace<insize) continue;

      // create local file
      sprintf(localfilename, "%s/LipMiniAnalysisTempFile%010.0f", dirtmp, r.Rndm()*1e10);
      std::ofstream outfile(localfilename, std::ios_base::binary);

      // return to begin of original file
      infile.seekg(0, std::ios::beg);

      // copy files
      outfile << infile.rdbuf();

      // check errors
      if (!outfile.good()) {
        printf("Copy File Error!\n");
        error=1;
      } else error=0;

      // close files
      infile.close();
      outfile.close();

      // check size of local file
      std::ifstream outtestfile(localfilename, std::ios_base::binary|std::ios::ate);
      long int outtestsize = (long int) outtestfile.tellg();
      printf("Copy Size : %ld\n", outtestsize);

      // close file
      outtestfile.close();

      // test file sizes
      // if original and local file have the same size then we are done
      // (we should compare them byte by byte, but that takes time...)
      if (outtestsize==insize && !error) {
        UseTempFile = 1;
        break;
      } else {
        printf("Removing Error Local File: %s\n", localfilename);
        std::remove(localfilename);
      }
    }

    // time
    time_t rawtimeC;
    time(&rawtimeC);
    cout << "  Read time was " << difftime(rawtimeC, rawtimeA) << " s " << endl;
    cout << "  Speed was " << insize/difftime(rawtimeC, rawtimeA) << " B/s " << endl;
// comentar daqui para cima e comentar o std::remove(localfilename); !!!

    // use original file in case of problems
    if (!UseTempFile) sprintf(localfilename, "%s", Input.Name(f).c_str());

    // open nTuple
    //nTuple->Input((char *)Input.Name(f).c_str());
    nTuple->Input(localfilename);
    nTuple->Init();
    DefaultValues();
    UserValues();
    CommandLineOptions(0);

    // loop over events
    Loop();

    // close file
    delete nTuple;

    // remove localfile
    printf("Used File: %s\n", localfilename);
    if (UseTempFile) {
      printf("Error File Removing: %s\n", localfilename);
      std::remove(localfilename);
    }

    // time
    time_t rawtimeB;
    time(&rawtimeB);
    cout << "  Partial time was " << difftime(rawtimeB, rawtimeA) << " s " << endl;
    cout << endl;
  }

  // delete dummy canvas
  delete mydummycanvas;

   //...Compute errors for total background
  for (Int_t i_syst=0; i_syst<Syst.size(); ++i_syst) {
    for(int j = 0; j <= MaxCuts; j++) {
      double err = 0.;
      for(int i=1; i < MonteCarlo.size(); i++) {
//        if (MonteCarlo[i].type()==1 && MonteCarlo[i].nSelEvt(0)>0) {
//
//        // just for backgrounds with events!!!
//        if (MonteCarlo[i].type()==1 && MonteCarlo[i].nSelEvt(j)>0) {
        // considering also backgrounds without events for the error evaluation:
        if (MonteCarlo[i].type()==1) {
          err += MonteCarlo[i].nSelWeightedEvtErr(i_syst, j)*MonteCarlo[i].nSelWeightedEvtErr(i_syst, j);
        }
      }
      MonteCarlo[0].SetnSelWeightedEvtErr(i_syst, j,sqrt(err));
    }
  }

  //...Compute total number of generated events for total background
  int TotGenEvt = 0;
  for(int i=1; i < MonteCarlo.size(); i++) {
    if (MonteCarlo[i].type()==1 && MonteCarlo[i].nSelEvt(0, 0)>0) {
      TotGenEvt += MonteCarlo[i].nGenEvt();
    }
  }
  MonteCarlo[0].SetnGenEvt(TotGenEvt);

  for (Int_t i_syst=0; i_syst<Syst.size(); ++i_syst)
    PrintSummary(3, i_syst);

/* Comentado para performance; descomentar se a análise necessitar de pdf e likelihoods
  ProcessPdf();
*/
  PostLoopCalculations();
  WriteHistograms();

  // Prints date, time
  time(&rawtimeF);
  timeinfo = localtime(&rawtimeF);
  strftime(timestring, 100, "%a %d %b %Y %T", timeinfo);
  cout << "  Now is " << timestring << endl;
  cout << "  Total time was " << difftime(rawtimeF, rawtimeI) << " s " << endl;

}

// #############################################################################
LipMiniAnalysis::~LipMiniAnalysis() {
// #############################################################################

  if (!fChain) return;
  delete fChain->GetCurrentFile();

  delete nTuple;
  OutputFile.clear();
}

// #############################################################################
void LipMiniAnalysis::LeptonIsolationDeltaR(){
// #############################################################################
//
//  purpose: calculation of isolation angle (in R) of leptons
//
//  authors: fveloso
//  first version: 10.may.2007
//
//  last change: 18.may.2007
//  by: nfcastro
//
// #############################################################################

  double IsoDeltaR;
  double tmp;

  for(int i = 0; i < LeptonVec.size(); ++i){
    IsoDeltaR = 999.;
    tmp = 999.;

    for(int j = 0; j < JetVec.size(); j++) {
      tmp = JetVec[j].DeltaR(LeptonVec[i]);
      if( tmp < IsoDeltaR ) IsoDeltaR = tmp;
    }

    for(int j = 0; j < LeptonVec.size(); j++) {
      tmp = LeptonVec[j].DeltaR(LeptonVec[i]);
      if( tmp < IsoDeltaR && j != i ) IsoDeltaR = tmp;
    }

    for(int j = 0; j < TauVec.size(); j++) {
      tmp = TauVec[j].DeltaR(LeptonVec[i]);
      if( tmp < IsoDeltaR ) IsoDeltaR = tmp;
    }

    for(int j = 0; j < PhotonVec.size(); j++) {
      tmp = PhotonVec[j].DeltaR(LeptonVec[i]);
      if( tmp < IsoDeltaR ) IsoDeltaR = tmp;
    }

    if(doGoodObjSelection==1) LeptonVec[i].SetIsoDeltaR(IsoDeltaR);
  }
}

// #############################################################################
void LipMiniAnalysis::PhotonIsolationDeltaR(){
// #############################################################################
//
//  purpose: calculation of isolation angle (in R) of photons
//
//  authors: fveloso
//  first version: 10.may.2007
//
//  last change: 24.apr.2010
//  by: nfcastro
//
// #############################################################################

  double IsoDeltaR;
  double tmp;

  for(int i = 0; i < PhotonVec.size(); ++i){
    IsoDeltaR = 999.;
    tmp = 999.;

    for(int j = 0; j < JetVec.size(); j++) {
      tmp = JetVec[j].DeltaR(PhotonVec[i]);
      if (tmp < IsoDeltaR) IsoDeltaR = tmp;
    }

    for(int j = 0; j < LeptonVec.size(); j++) {
      tmp = LeptonVec[j].DeltaR(PhotonVec[i]);
      if (tmp < IsoDeltaR) IsoDeltaR = tmp;
    }

    for(int j = 0; j < TauVec.size(); j++) {
      tmp = TauVec[j].DeltaR(PhotonVec[i]);
      if( tmp < IsoDeltaR ) IsoDeltaR = tmp;
    }

    for(int j = 0; j < PhotonVec.size(); j++) {
      tmp = PhotonVec[j].DeltaR(PhotonVec[i]);
      if (tmp < IsoDeltaR && j != i) IsoDeltaR = tmp;
    }

    PhotonVec[i].SetIsoDeltaR(IsoDeltaR);
  }
}

// #############################################################################
void LipMiniAnalysis::JetIsolationDeltaR(){
// #############################################################################
//
//  purpose: calculation of isolation angle (in R) of jets
//
//  authors: fveloso
//  first version: 10.may.2007
//
//  last change: 24.apr.2010
//  by: nfcastro
//
// #############################################################################

  double IsoDeltaR;
  double tmp;

  for (int i = 0; i < JetVec.size(); ++i){
    IsoDeltaR = 999.;
    tmp = 999.;

    for (int j = 0; j < JetVec.size(); j++) {
      tmp = JetVec[j].DeltaR(JetVec[i]);
      if (tmp < IsoDeltaR && j != i) IsoDeltaR = tmp;
    }

    for (int j = 0; j < LeptonVec.size(); j++) {
      tmp = LeptonVec[j].DeltaR(JetVec[i]);
      if (tmp < IsoDeltaR) IsoDeltaR = tmp;
    }

    for(int j = 0; j < TauVec.size(); j++) {
      tmp = TauVec[j].DeltaR(JetVec[i]);
      if( tmp < IsoDeltaR ) IsoDeltaR = tmp;
    }

    for (int j = 0; j < PhotonVec.size(); j++) {
      tmp = PhotonVec[j].DeltaR(JetVec[i]);
      if (tmp < IsoDeltaR) IsoDeltaR = tmp;
    }

    JetVec[i].SetIsoDeltaR(IsoDeltaR);
  }
}

// #############################################################################
void LipMiniAnalysis::TauIsolationDeltaR(){
// #############################################################################
//
//  purpose: calculation of isolation angle (in R) of taus
//
//  authors: nfcastro
//  first version: 18.may.2007
//
//  last change: 18.may.2007
//  by: nfcastro
//
// #############################################################################

  double IsoDeltaR;
  double tmp;

  for(int i = 0; i < TauVec.size(); ++i){
    IsoDeltaR = 999.;
    tmp = 999.;

    for(int j = 0; j < JetVec.size(); j++) {
      tmp = JetVec[j].DeltaR(TauVec[i]);
      if (tmp < IsoDeltaR) IsoDeltaR = tmp;
    }

    for(int j = 0; j < TauVec.size(); j++) {
      tmp = TauVec[j].DeltaR(TauVec[i]);
      if( tmp < IsoDeltaR && j != i ) IsoDeltaR = tmp;
    }

    for(int j = 0; j < PhotonVec.size(); j++) {
      tmp = PhotonVec[j].DeltaR(TauVec[i]);
      if (tmp < IsoDeltaR) IsoDeltaR = tmp;
    }

    for(int j = 0; j < LeptonVec.size(); j++) {
      tmp = LeptonVec[j].DeltaR(TauVec[i]);
      if (tmp < IsoDeltaR) IsoDeltaR = tmp;
    }

    TauVec[i].SetIsoDeltaR(IsoDeltaR);
  }
}

// #############################################################################
void LipMiniAnalysis::FillElectronVec(){
// #############################################################################
//
//  purpose: fill vectors of TLorentzVectors-with-index for all electrons
//
//  authors: fveloso
//  first version: 10.may.2007
//
//  last change: 17.Nov.2012
//  by: A.Onofre
//
// #############################################################################

  ElectronVec.clear();

  for(Int_t i = 0; i < LeptonVec.size(); ++i){
    if( abs(LeptonVec[i].isb) == 11 ) {
      ElectronVec.push_back(LeptonVec[i]);
    }
  }

  // AO 17 Nov 2012 ===== sort Electrons inside the vectors
  sort(ElectronVec.begin(), ElectronVec.end(), LorentzVecComp);
  // AO 17 Nov 2012 ===== sort Electrons inside the vectors

}

// #############################################################################
void LipMiniAnalysis::FillMuonVec(){
// #############################################################################
//
//  purpose: fill vectors of TLorentzVectors-with-index for all muons
//
//  authors: fveloso
//  first version: 10.may.2007
//
//  last change: 17.Nov.2012
//  by: A.Onofre
//
// #############################################################################

  MuonVec.clear();

  for(Int_t i = 0; i < LeptonVec.size(); ++i){
    if( abs(LeptonVec[i].isb) == 13 ) {
      MuonVec.push_back(LeptonVec[i]);
    }
  }

  // AO 17 Nov 2012 ===== sort Muons inside the vectors
  sort(MuonVec.begin(), MuonVec.end(), LorentzVecComp);
  // AO 17 Nov 2012 ===== sort Muons inside the vectors

}

// #############################################################################
void LipMiniAnalysis::FillBTaggedJetVec(){
// #############################################################################
//
//  purpose: fill vectors of TLorentzVectors-with-index for all b jets
//
//  authors: fveloso
//  first version: 10.may.2007
//
//  last change: 17.Nov.2012
//  by: A.Onofre
//
// #############################################################################

  BTaggedJetVec.clear();

  for(Int_t i = 0; i < JetVec.size(); ++i){
    if( abs(JetVec[i].isb) == 5 ) {
      BTaggedJetVec.push_back(JetVec[i]);
    }
  }

  // AO 17 Nov 2012 ===== sort b-jets inside the vectors
  sort(BTaggedJetVec.begin(), BTaggedJetVec.end(), LorentzVecComp);
  // AO 17 Nov 2012 ===== sort b-jets inside the vectors

}

// #############################################################################
void LipMiniAnalysis::FillNonTaggedJetVec(){
// #############################################################################
//
//  purpose: fill vectors of TLorentzVectors-with-index for all non b jets
//
//  authors: fveloso
//  first version: 10.may.2007
//
//  last change: 12.Nov.2012
//  by: A.Onofre
//
// #############################################################################

  NonTaggedJetVec.clear();

  for(Int_t i = 0; i < JetVec.size(); ++i){
    if( abs(JetVec[i].isb) != 5 ) {
      NonTaggedJetVec.push_back(JetVec[i]);
    }
  }

  // AO 17 Nov 2012 ===== sort Non b-jets inside the vectors
  sort(NonTaggedJetVec.begin(), NonTaggedJetVec.end(), LorentzVecComp);
  // AO 17 Nov 2012 ===== sort Non b-jets inside the vectors

}

// #############################################################################
void LipMiniAnalysis::FillAllVectors(){
// #############################################################################
//
//  purpose: main function that calls the routines that fill Vectors of
//  TLorentzVectors-with-index for reconstructed detector objects
//  (muons, jets, electrons, etc...)
//
//  authors: nfcastro
//  first version: ??.???.??
//
//  last change: 17.Nov.2012
//  by: A.Onofre
//
// #############################################################################

  // fill vectors with particles
  nTuple->FillVectors(	TruthVec, 
			TauVec, 
			LeptonVec, 
			PhotonVec, 
			JetVec, 
			Vtx,
    			LightJetCalib, 
			BJetCalib, 
			BTagCut, 
			JetOverlap, 
			PtCutJet, 
			PtCutEle, 
			PtCutMuo, 
			RecoType);

  // sort particles inside the vectors
  sort(LeptonVec.begin(), LeptonVec.end(), LorentzVecComp);
  sort(PhotonVec.begin(), PhotonVec.end(), LorentzVecComp);
  sort(JetVec.begin(), JetVec.end(), LorentzVecComp);

  // isolation angles
  LeptonIsolationDeltaR();
  PhotonIsolationDeltaR();
  JetIsolationDeltaR();
  TauIsolationDeltaR();

  // AO 17 Nov 2012 ========== fill vectors with subtypes of particles and sort particles inside the vectors
  FillMuonVec();
  FillElectronVec();
  FillBTaggedJetVec();
  FillNonTaggedJetVec();
  // AO 17 Nov 2012 ========== fill vectors with subtypes of particles and sort particles inside the vectors

  // other stuff to be filled only after nTuple->FillVectors is called
  Isub = nTuple->Isub();
  LumiBlock = nTuple->LumiBlock();
  RunNumber = nTuple->RunNumber();
  EveNumber = nTuple->EveNumber();
  TruthEleNumber = nTuple->TruthEleNumber();
  TruthMuonNumber = nTuple->TruthMuonNumber();
  ElectronTrigger = nTuple->ElectronTrigger();
  MuonTrigger = nTuple->MuonTrigger();
  Cosmic = nTuple->Cosmic();
  HforFlag = nTuple->HforFlag();
  Ht_Mini = nTuple->Ht_Mini();
  massInv_LL_Mini = nTuple->massInv_LL_Mini();
  jet_n_Mini = nTuple->jet_n_Mini();
  EleMuoOverlap = nTuple->EleMuoOverlap();
  GoodRL = nTuple->GoodRL();
  // MissPt() must be called before MissPx() and MissPy()
  MissPt = nTuple->MissPt();
  MissPx = nTuple->MissPx();
  MissPy = nTuple->MissPy();
  Weight = nTuple->Weight();
  Sphericity = nTuple->Sphericity();
  Aplanarity = nTuple->Aplanarity();
  Planarity  = nTuple->Planarity();

}

// #############################################################################
void LipMiniAnalysis::DefineSamples(){
// #############################################################################
//
//  purpose: the user should declare here the ROOT samples that will be
//  analysed
//
//  notes: should be replaced by the user!
//
//  authors: fveloso
//  first version: 18.nov.2006
//
//  last change: 18.nov.2006
//  by: fveloso
//
// #############################################################################

  //cout << "LipMiniAnalysis::DefineSamples()" << endl;

}

// #############################################################################
void LipMiniAnalysis::ReadFiles(const char* ntu, const char * file, int type, double run, double lum, int nGenEvt, std::string title, int MaxCuts) {
// #############################################################################
//
//  purpose: to be used inside DefineSamples()
//
//  author: fveloso
//  first version: 06.jan.2011
//
//  last change:07.jan.2011
//  by: fveloso
//
// #############################################################################

  ifstream in;
  char s[500];
  //sprintf(s, "../lists/%s", file);
  sprintf(s, "%s", file);
  in.open(s);
  while (1) {
    std::string rootfile;
    in >> rootfile;
    if (!in.good()) break;
    Input.File(ntu, rootfile);
  }
  in.close();
  TMonteCarlo mc(type, run, lum, nGenEvt, title, MaxCuts);
  MonteCarlo.push_back(mc);
}


// #############################################################################
void LipMiniAnalysis::DefineSignalPdfSamples(){
// #############################################################################
//
//  purpose: the user should declare here the ROOT samples that will be
//  used as probability density functions
//
//  notes: should be replaced by the user!
//
//  authors: fveloso
//  first version: 18.jun.2007
//
// #############################################################################

  SignalPdfFiles.push_back("signal_pdf.root");

}

// #############################################################################
void LipMiniAnalysis::BookHistograms(THistos &){
// #############################################################################
//
//  purpose: to book histograms before the event's loop
//
//  notes: should be replaced by the user!
//
//  authors: fveloso
//  first version: 15.nov.2006
//
//  last change: 15.nov.2006
//  by: fveloso
//
// #############################################################################

}

// #############################################################################
void LipMiniAnalysis::BookPdfHistos(Int_t i_syst) {
// #############################################################################
//
//  purpose: to book histograms to be used as probability density functions  
//
//  notes: should be replaced by the user!
//
//  authors: fveloso
//  first version: 15.nov.2006
//
//  last change: 15.nov.2006
//  by: fveloso
//
// #############################################################################

}

// #############################################################################
void LipMiniAnalysis::BookLikeHistos(){
// #############################################################################
//
//  purpose: to book histograms to be used as discriminant variables
//
//  notes: should be replaced by the user!
//
//  authors: fveloso
//  first version: 15.nov.2006
//
//  last change: 15.nov.2006
//  by: fveloso
//
// #############################################################################

  LogELikeSHist = new TH1D("LogELikeS", "Ln(L_S)", 60, -30, 0);
  LogELikeSHist->Sumw2();
  LogELikeBHist = new TH1D("LogELikeB", "Ln(L_B)", 60, -30, 0);
  LogELikeBHist->Sumw2();
  LogELikelihoodHist = new TH1D("LogELikelihood", "Ln(L_S/L_B)", 100, -10, 10);
  LogELikelihoodHist->Sumw2();

  Log10LikeSHist = new TH1D("Log10LikeS", "Log(L_S)", 50, -5, 0);
  Log10LikeSHist->Sumw2();
  Log10LikeBHist = new TH1D("Log10LikeB", "Log(L_B)", 50, -5, 0);
  Log10LikeBHist->Sumw2();
  Log10LikelihoodHist = new TH1D("Log10Likelihood", "Log(L_S/L_B)", 45, -3, 3);
  Log10LikelihoodHist->Sumw2();

}

// #############################################################################
void LipMiniAnalysis::ProcessBookLikeHistos(){
// #############################################################################
//
//  purpose: to copie booked histograms to be used as discriminant variables
//
//  authors: fveloso
//  first version: 18.jun.2007
//
// #############################################################################

  TDirectory *dir = gDirectory;
  for (int i=0; i<SignalPdfFiles.size(); ++i){
    char t[10];
    sprintf(t, "LikesSig%02i", i);
    TDirectory *likedir = new TDirectory(t,t);
    likedir->cd();

    LogELikeSHistVec[i] = (TH1D*)LogELikeSHist->Clone("LogELikeS");
    LogELikeBHistVec[i] = (TH1D*)LogELikeBHist->Clone("LogELikeB");
    LogELikelihoodHistVec[i] = (TH1D*)LogELikelihoodHist->Clone("LogELikelihood");

    Log10LikeSHistVec[i] = (TH1D*)Log10LikeSHist->Clone("Log10LikeS");
    Log10LikeBHistVec[i] = (TH1D*)Log10LikeBHist->Clone("Log10LikeB");
    Log10LikelihoodHistVec[i] = (TH1D*)Log10LikelihoodHist->Clone("Log10Likelihood");
    dir->cd();
  }
}

// #############################################################################
void LipMiniAnalysis::FillHistograms(THistos &){
// #############################################################################
//
//  purpose: to fill histograms inside the event's loop
//
//  notes: should be replaced by the user!
//
//  authors: fveloso
//  first version: 15.nov.2006
//
//  last change: 15.nov.2006
//  by: fveloso
//
// #############################################################################

}

// #############################################################################
void LipMiniAnalysis::PdfValues(){
// #############################################################################
//
//  purpose: to define the pdf values of each event
//
//  notes: should be replaced by the user!
//
//  authors: fveloso
//  first version: 20.nov.2006
//
//  last change: 21.nov.2006
//  by: fveloso
//
// #############################################################################

}

// #############################################################################
void LipMiniAnalysis::FillPdf(Int_t i_syst){
// #############################################################################
//
//  purpose: to fill pdf inside the event's loop
//
//  authors: fveloso
//  first version: 21.nov.2006
//
//  last change: 21.nov.2006
//  by: fveloso
//
// #############################################################################

  if (PdfVec[i_syst].size() > 0) {
    for (int j=0; j<PdfVec[i_syst].size(); ++j){
      PdfVec[i_syst][j]->Fill(ThisEventPdfValues[j],Weight);
    }
  }
}

// #############################################################################
void LipMiniAnalysis::CreateOutputNtuple(){
// #############################################################################
//
//  purpose: create a ntuple to store surviving events
//
//  authors: fveloso
//  first version: 27.jun.2007
//
// #############################################################################

  TDirectory *dir = gDirectory;

  for (Int_t i_syst=0; i_syst<Syst.size(); ++i_syst) {
    char ntupledef[1000] = "Isub:Weight";
    for (int j=0; j<PdfVec[i_syst].size(); ++j){
      sprintf(ntupledef, "%s:pdf%02i", ntupledef, j);
    }
    OutputFile[i_syst]->cd();
    TNtuple* pOutputNtuple = new TNtuple("SelectedEvents","Events after sequential selection",ntupledef);
    OutputNtuple.push_back(pOutputNtuple);
  }

  dir->cd();
}

// #############################################################################
void LipMiniAnalysis::FillOutputNtuple(Int_t i_syst) {
// #############################################################################
//
//  purpose: to fill the ntuple with the surviving events
//
//  last change: 30.nov.2009
//  by: nfcastro
//
// #############################################################################

  //cout << " FillOutputNtuple() i_syst=" << i_syst << " PdfVec.size()=" << PdfVec.size() << " PdfVec[i_syst].size()=" << PdfVec[i_syst].size() << endl;

       if (PdfVec[i_syst].size() ==  0) OutputNtuple[i_syst]->Fill(Isub,Weight);
  else if (PdfVec[i_syst].size() ==  1) OutputNtuple[i_syst]->Fill(Isub,Weight,ThisEventPdfValues[0]);
  else if (PdfVec[i_syst].size() ==  2) OutputNtuple[i_syst]->Fill(Isub,Weight,ThisEventPdfValues[0],ThisEventPdfValues[1]);
  else if (PdfVec[i_syst].size() ==  3) OutputNtuple[i_syst]->Fill(Isub,Weight,ThisEventPdfValues[0],ThisEventPdfValues[1],ThisEventPdfValues[2]);
  else if (PdfVec[i_syst].size() ==  4) OutputNtuple[i_syst]->Fill(Isub,Weight,ThisEventPdfValues[0],ThisEventPdfValues[1],ThisEventPdfValues[2],ThisEventPdfValues[3]);
  else if (PdfVec[i_syst].size() ==  5) OutputNtuple[i_syst]->Fill(Isub,Weight,ThisEventPdfValues[0],ThisEventPdfValues[1],ThisEventPdfValues[2],ThisEventPdfValues[3],ThisEventPdfValues[4]);
  else if (PdfVec[i_syst].size() ==  6) OutputNtuple[i_syst]->Fill(Isub,Weight,ThisEventPdfValues[0],ThisEventPdfValues[1],ThisEventPdfValues[2],ThisEventPdfValues[3],ThisEventPdfValues[4],ThisEventPdfValues[5]);
  else if (PdfVec[i_syst].size() ==  7) OutputNtuple[i_syst]->Fill(Isub,Weight,ThisEventPdfValues[0],ThisEventPdfValues[1],ThisEventPdfValues[2],ThisEventPdfValues[3],ThisEventPdfValues[4],ThisEventPdfValues[5],ThisEventPdfValues[6]);
  else if (PdfVec[i_syst].size() ==  8) OutputNtuple[i_syst]->Fill(Isub,Weight,ThisEventPdfValues[0],ThisEventPdfValues[1],ThisEventPdfValues[2],ThisEventPdfValues[3],ThisEventPdfValues[4],ThisEventPdfValues[5],ThisEventPdfValues[6],ThisEventPdfValues[7]);
  else if (PdfVec[i_syst].size() ==  9) OutputNtuple[i_syst]->Fill(Isub,Weight,ThisEventPdfValues[0],ThisEventPdfValues[1],ThisEventPdfValues[2],ThisEventPdfValues[3],ThisEventPdfValues[4],ThisEventPdfValues[5],ThisEventPdfValues[6],ThisEventPdfValues[7],ThisEventPdfValues[8]);
  else if (PdfVec[i_syst].size() == 10) OutputNtuple[i_syst]->Fill(Isub,Weight,ThisEventPdfValues[0],ThisEventPdfValues[1],ThisEventPdfValues[2],ThisEventPdfValues[3],ThisEventPdfValues[4],ThisEventPdfValues[5],ThisEventPdfValues[6],ThisEventPdfValues[7],ThisEventPdfValues[8],ThisEventPdfValues[9]);
  else if (PdfVec[i_syst].size() == 11) OutputNtuple[i_syst]->Fill(Isub,Weight,ThisEventPdfValues[0],ThisEventPdfValues[1],ThisEventPdfValues[2],ThisEventPdfValues[3],ThisEventPdfValues[4],ThisEventPdfValues[5],ThisEventPdfValues[6],ThisEventPdfValues[7],ThisEventPdfValues[8],ThisEventPdfValues[9],ThisEventPdfValues[10]);
  else if (PdfVec[i_syst].size() == 12) OutputNtuple[i_syst]->Fill(Isub,Weight,ThisEventPdfValues[0],ThisEventPdfValues[1],ThisEventPdfValues[2],ThisEventPdfValues[3],ThisEventPdfValues[4],ThisEventPdfValues[5],ThisEventPdfValues[6],ThisEventPdfValues[7],ThisEventPdfValues[8],ThisEventPdfValues[9],ThisEventPdfValues[10],ThisEventPdfValues[11]);
  else if (PdfVec[i_syst].size()  > 12) cout << "ERROR: number of p.d.f.s excedeed 12: " << PdfVec[i_syst].size() << endl;

}

// #############################################################################
void LipMiniAnalysis::ProcessPdf(){
// #############################################################################
//
//  purpose: to normalize and smoth the pdf histograms
//
//  authors: fveloso
//  first version: 20.nov.2006
//
//  last change: 20.nov.2006
//  by: fveloso
//
// #############################################################################

  // loop over systematics
  for (Int_t i_syst=0; i_syst<Syst.size(); ++i_syst) {
    if (PdfVec[i_syst].size() > 0) {
      cout << endl;
      cout << " #################################################################################" << endl;
      cout << endl;
      cout << "  Probability density functions defined:" << endl;
      cout << endl;
      for (int j=0; j<PdfVec[i_syst].size(); ++j){
        // protection: if bin content is less or equal 0, change to 1e-30
        int err = 0;
        int nbins = PdfVec[i_syst][j]->GetXaxis()->GetNbins(); 
        for (int k = 1; k <= nbins; ++k) {
          if (PdfVec[i_syst][j]->GetBinContent(k) <= 0.) {
            PdfVec[i_syst][j]->SetBinContent(k, 1e-30);
            err = 1;
          }
        }

        // print info
        cout << "  " << j << " : \"" << PdfVec[i_syst][j]->GetTitle() << "\" ";
        if (err == 1) cout << " WITH NON POSITIVE BINS (CHANGED TO 1E-30)!!!" ;
        cout << endl;

        // smooth
        PdfVec[i_syst][j]->Smooth(PdfSmooth[j]);
        cout << " p.d.f. " << j << " smoothed with factor " << PdfSmooth[j] << endl;

        // normalization
        Double_t scale = 1/PdfVec[i_syst][j]->Integral();
        PdfVec[i_syst][j]->Scale(scale);
      }
      cout << endl;
   }
  }
}

// #############################################################################
void LipMiniAnalysis::ReadPdf(){
// #############################################################################
//
//  purpose: to read the pdf histograms made previously
//
//  authors: fveloso
//  first version: 20.nov.2006
//
//  last change: 18.jun.2007
//  by: fveloso
//
// #############################################################################

  if (PdfVec[0].size() > 0) {
    cout << endl;
    cout << "################################################################################" << endl;
    cout << endl;
    cout << "  Loading probability density functions:" << endl;
    cout << endl;

    TDirectory *dir = gDirectory;
    char t[5];

    cout << "  Background file \"" << BackgroundPdfFile.c_str() << "\":" << endl;
    TFile* BackgroundPdf = new TFile(BackgroundPdfFile.c_str());
    for (int j=0; j<PdfVec[0].size(); ++j){
      TH1D *pdfb = (TH1D*)BackgroundPdf->Get(pdf(t,j));
      BackgroundPdfVec.push_back(pdfb);
      cout << "  " << j << " : \"" << BackgroundPdfVec[j]->GetTitle() << "\" " << endl;
    }

    for (int i=0; i<SignalPdfFiles.size(); ++i){
      cout << endl;
      cout << "  Signal file \""<< SignalPdfFiles[i] << "\":" << endl;
      TFile* SignalPdf = new TFile(SignalPdfFiles[i]);
      for (int j=0; j<PdfVec[0].size(); ++j){
        TH1D *pdfs = (TH1D*)SignalPdf->Get(pdf(t,j));
        SignalPdfVec[i].push_back(pdfs);
       cout << "  " << j << " : \"" << SignalPdfVec[i][j]->GetTitle() << "\" " << endl;
      }
    }

    dir->cd();
  }
}

// #############################################################################
void LipMiniAnalysis::ComputeAndFillLikelihood(){
// #############################################################################
//
//  purpose: to compute the discriminant variables and to fill the
//  corresponding histograms
//
//  authors: fveloso
//  first version: 21.nov.2006
//
//  last change: 18.jun.2007
//  by: fveloso
//
// #############################################################################

  if (PdfVec[0].size() > 0) {
    for (int i=0; i<SignalPdfFiles.size(); ++i){
      double LikeSValue = 1;
      double LikeBValue = 1;
      for (int j=0; j<PdfVec[0].size(); ++j){
        LikeSValue *= SignalPdfVec[i][j]->GetBinContent(SignalPdfVec[i][j]->FindBin(ThisEventPdfValues[j]));
        LikeBValue *= BackgroundPdfVec[j]->GetBinContent(BackgroundPdfVec[j]->FindBin(ThisEventPdfValues[j]));
      }

      LogELikeSValue = log(LikeSValue);
      LogELikeBValue = log(LikeBValue);
      LogELikeLOverLikeBValue = log(LikeSValue/LikeBValue);

      Log10LikeSValue = log10(LikeSValue);
      Log10LikeBValue = log10(LikeBValue);
      Log10LikeLOverLikeBValue = log10(LikeSValue/LikeBValue);

      LogELikeSHistVec[i]->Fill(LogELikeSValue,Weight);   
      LogELikeBHistVec[i]->Fill(LogELikeBValue,Weight);   
      LogELikelihoodHistVec[i]->Fill(LogELikeLOverLikeBValue,Weight);   

      Log10LikeSHistVec[i]->Fill(Log10LikeSValue,Weight);
      Log10LikeBHistVec[i]->Fill(Log10LikeBValue,Weight);
      Log10LikelihoodHistVec[i]->Fill(Log10LikeLOverLikeBValue,Weight);

      if (i == 0) {
        LogELikeSHist->Fill(LogELikeSValue,Weight);   
        LogELikeBHist->Fill(LogELikeBValue,Weight);   
        LogELikelihoodHist->Fill(LogELikeLOverLikeBValue,Weight);   

        Log10LikeSHist->Fill(Log10LikeSValue,Weight);
        Log10LikeBHist->Fill(Log10LikeBValue,Weight);
        Log10LikelihoodHist->Fill(Log10LikeLOverLikeBValue,Weight);
      }
    }
  }
}

// #############################################################################
void LipMiniAnalysis::DoCuts(){
// #############################################################################
//
//  purpose: to do the selection cuts
//
//  notes: should be replaced by the user!
//
//  authors: fveloso
//  first version: 15.nov.2006
//
//  last change: 16.nov.2006
//  by: fveloso
//
// #############################################################################

}

// #############################################################################
void LipMiniAnalysis::CreateOutputFile(){
// #############################################################################
//
//  purpose: create a file for output
//
//  authors: fveloso
//  first version: 26.fev.2007
//
//  last change: 18.jun.2007
//  by: fveloso
//
// #############################################################################

  // loop over systematics
  for (Int_t i_syst=0; i_syst<Syst.size(); ++i_syst) {
    char FileName[1000];
    sprintf(FileName, "%s%s.root", SystName[i_syst].c_str(), OutputFileName.c_str());
    TFile* pOutputFile = new TFile(FileName,"recreate");
    OutputFile.push_back(pOutputFile);
  }

  OutputFile[0]->cd();

}

// #############################################################################
void LipMiniAnalysis::WriteHistograms(){
// #############################################################################
//
//  purpose: write all histograms to the output file
//
//  authors: fveloso
//  first version: 19.nov.2006
//
//  last change: 26.fev.2007
//  by: fveloso
//
// #############################################################################

  for (Int_t i_syst=0; i_syst<Syst.size(); ++i_syst) {
    OutputFile[i_syst]->Write();
  }
}

// #############################################################################
void LipMiniAnalysis::PrintSummary(int OutputLevel, Int_t i_syst){
// #############################################################################
//
//  purpose: print a summary of the analysis
//
//  authors: fveloso
//  first version: 19.nov.2006
//
//  last change: 19.nov.2006
//  by: fveloso
//
// #############################################################################

  // Define the output precision:
  cout.precision(5);

  cout << endl;
  cout << " ################################################################################" << endl;
  cout << endl;
  cout << "Summary for " << Syst[i_syst] << " " << SystName[i_syst] << endl;
  cout << endl;
  cout << "  Luminosity set to: " << Luminosity << " fb^-1" << std::endl;

  // remember number of events selected on each cut and cuts description

  if (OutputLevel==1 || OutputLevel==0) {

    cout << endl;
    if(isData==0) cout << "  Total number of analysed events (unweighted): " << MonteCarlo[0].nSelEvt(i_syst, 0) << std::endl;
    else          cout << "  Total number of analysed events: "              << MonteCarlo[0].nSelEvt(i_syst, 0) << std::endl;
    if(isData==0) cout << "  Total number of analysed events (weighted)  : " << MonteCarlo[0].nSelWeightedEvt(i_syst, 0) << std::endl;
    cout << endl;

    for(int i = 0; i <= MaxCuts; ++i) {
//      Double_t sel_eff = (Sumweights_total > 0.) ? 100.*(npassl[i]/Sumweights_total) : 0.;
//      Double_t err_eff = (Sumweights_total > 0.) ? 100.*(sqrt(npassl[i])/Sumweights_total) : 0.;
//      cout << "LEVEL " << i << ": " << npassl[i] << " out of "<< Sumweights_total << " events were selected -> eff = "
//        << sel_eff << " +- " << err_eff << " %" <<endl;
      cout << "  Level " << i << ": " << MonteCarlo[0].nSelWeightedEvt(i_syst, i) << " out of "<< MonteCarlo[0].nSelWeightedEvt(i_syst, 0) << " events were selected -> eff = "
        << MonteCarlo[0].nSelEvt(i_syst, i)/MonteCarlo[0].nSelEvt(i_syst, 0)*100.
        << " +- " << sqrt(MonteCarlo[0].nSelEvt(i_syst, i))/MonteCarlo[0].nSelEvt(i_syst, 0)*100 << " %" <<endl;
    }
    cout << endl;
    std::cout << "  Events description (after final selection):"<< std::endl;
    for(int i = 0; i < MonteCarlo.size(); i++) {
      if (MonteCarlo[i].nSelEvt(i_syst, 0)>0) {
        std::cout << "      " << MonteCarlo[i].title() << " (run " << MonteCarlo[i].run() << "): \t" << MonteCarlo[i].nSelWeightedEvt(i_syst, MaxCuts) << " +- " << MonteCarlo[i].nSelWeightedEvtErr(i_syst, MaxCuts)
                  << "\t (" << MonteCarlo[i].nSelEvt(i_syst, MaxCuts) << " events with weight " << MonteCarlo[i].nSelWeightedEvt(i_syst, MaxCuts)/MonteCarlo[i].nSelEvt(i_syst, MaxCuts) << ")" <<  endl;
        std::cout << "       \t \t --> overall efficiency (selected/generated) = "
                  << MonteCarlo[i].nSelEvt(i_syst, MaxCuts)/MonteCarlo[i].nGenEvt()*100 << " +- " << MonteCarlo[i].nSelWeightedEvtErr(i_syst, MaxCuts)/MonteCarlo[i].nGenEvt()*100 << " %" << endl;
      }
    }
    std::cout << "    total: " << MonteCarlo[0].nSelWeightedEvt(i_syst, MaxCuts) << " +- " << MonteCarlo[0].nSelWeightedEvtErr(i_syst, MaxCuts) << endl;

  }

  if (OutputLevel==2 || OutputLevel==0) {
/*
    cout << std::left << std::setw(40) << " ";
    for(int i = 0; i <= MaxCuts; i++) {
      cout << std::right << std::setw(39) << "LEVEL " << i << "  ";
    }
    cout << endl;
    for(int i=0; i < MonteCarlo.size(); i++) {
      if (MonteCarlo[i].nSelEvt(0)>0) {
        cout << std::left << std::setw(40) << MonteCarlo[i].title() << ":";
        for(int j = 0; j <= MaxCuts; j++) {
          cout << std::dec << std::fixed << std::right << std::setw(10) << std::setprecision(1) << std::showpoint << MonteCarlo[i].nSelWeightedEvt(j)       << " +- "
               << std::dec << std::fixed << std::right << std::setw( 7) << std::setprecision(1) << std::showpoint << MonteCarlo[i].nSelWeightedEvtErr(j);
// NOTE: verify code!!!
//     << " ("
//               << std::dec << std::fixed << std::right << std::setw( 5) << std::setprecision(1) << std::showpoint << MonteCarlo[i].nSelWeightedEvtEff(j)    << " +- "
//               << std::dec << std::fixed << std::right << std::setw( 5) << std::setprecision(3) << std::showpoint << MonteCarlo[i].nSelWeightedEvtEffErr(j) << " %) ";
        }
      cout << endl;
      }
    }
    cout << endl;
*/

/*
    cout << std::left << std::setw(40) << MonteCarlo[0].title() << ":";
    for(int j = 0; j <= MaxCuts; j++) {
      double err = 0.;
      for(int i=1; i < MonteCarlo.size(); i++) {
        if (MonteCarlo[i].nSelEvt(0)>0) {
          err+= MonteCarlo[i].nSelWeightedEvtErr(j)*MonteCarlo[i].nSelWeightedEvtErr(j);
        }
      }
      err=sqrt(err);
      cout << std::dec << std::fixed << std::right << std::setw(10) << std::setprecision(1) << std::showpoint << MonteCarlo[0].nSelWeightedEvt(j)       << " +- "
           << std::dec << std::fixed << std::right << std::setw( 7) << std::setprecision(1) << std::showpoint << err                                    << " ("
           << std::dec << std::fixed << std::right << std::setw( 5) << std::setprecision(1) << std::showpoint << MonteCarlo[0].nSelWeightedEvtEff(j)    << " +- "
           << std::dec << std::fixed << std::right << std::setw( 5) << std::setprecision(3) << std::showpoint << MonteCarlo[0].nSelWeightedEvtEffErr(j) << " %) ";
    }
    cout << endl;
*/
  }

  if (OutputLevel==3 || OutputLevel==0) {
    for(int i = 0; i < MonteCarlo.size(); i++) {
      if (MonteCarlo[i].nSelEvt(i_syst, 0)>0) {

        cout << endl;
        cout << "  ----------------------------------------------------------------------------" << endl;
        cout << "  "
          << MonteCarlo[i].title()
          << " (RunNumber = "
          << std::dec << std::fixed << std::right << std::setw(0) << std::setprecision(0)
          << MonteCarlo[i].run()
          << "):"
          << endl;
        cout << "  ----------------------------------------------------------------------------" << endl;
        cout << endl;
       if (isData != 1){ 
          cout << "    Number of generated events (unweighted) = "
            << std::dec << std::fixed << std::right << std::setw(10) << std::setprecision(0)
            << MonteCarlo[i].nGenEvt() << endl;
       }
        cout << "    Number of analysed events (unweighted)  = "
          << std::dec << std::fixed << std::right << std::setw(10) << std::setprecision(0)
          << MonteCarlo[i].nSelEvt(i_syst, 0) << endl;
        cout << "    Number of analysed events (weighted)    = "
          << std::dec << std::fixed << std::right << std::setw(12) << std::setprecision(1)
          << MonteCarlo[i].nSelWeightedEvt(i_syst, 0) << endl;
        cout << "    Weight of analysed events               = "
          << std::dec << std::fixed << std::right << std::setw(16) << std::setprecision(5)
          << MonteCarlo[i].nSelWeightedEvt(i_syst, 0)/MonteCarlo[i].nSelEvt(i_syst, 0) << endl;
        cout << "    Luminosity of the sample                = "
          << std::dec << std::fixed << std::right << std::setw(16) << std::setprecision(5)
          << MonteCarlo[i].lum() << endl;
        cout << endl;

        cout << "    Level         Weighted             (Unweighted)" << endl;

        for(int l = 0; l <= MaxCuts; ++l) {
          cout
            // level
            << std::dec << std::fixed << std::right << std::setw(7) << std::setprecision(0)
            << l << ": "

            // number of selected events (weighted)
            << std::dec << std::fixed << std::right << std::setw(10) << std::setprecision(1)
            << MonteCarlo[i].nSelWeightedEvt(i_syst, l)

            // statistic error of the number of selected events (weighted)
            << " +- "
            << std::dec << std::fixed << std::right << std::setw(9) << std::setprecision(2)
            << MonteCarlo[i].nSelWeightedEvtErr(i_syst, l)

            // number of selected events (unweighted)
            << " ("
            << std::dec << std::fixed << std::right << std::setw(10) << std::setprecision(0)
            << MonteCarlo[i].nSelEvt(i_syst, l);

            if(isData != 1){
              // efficiency (of unweighted events)
              cout << " ; "
              << std::dec << std::fixed << std::right << std::setw(6) << std::setprecision(2)
              << MonteCarlo[i].nSelEvt(i_syst, l)/MonteCarlo[i].nGenEvt()*100

              // statistic error of the efficiency (of unweighted events)
              << " +- "
              << std::dec << std::fixed << std::right << std::setw(5) << std::setprecision(3)
              << sqrt(MonteCarlo[i].nSelEvt(i_syst, l))/MonteCarlo[i].nGenEvt()*100 << " % "
              << ")"        
              << endl;
            }else{
              // efficiency (of unweighted events)
              cout << " ; "
              << std::dec << std::fixed << std::right << std::setw(6) << std::setprecision(2)
              << MonteCarlo[i].nSelEvt(i_syst, l)/MonteCarlo[i].nSelEvt(i_syst, 0)*100

              // statistic error of the efficiency (of unweighted events)
              << " +- "
              << std::dec << std::fixed << std::right << std::setw(5) << std::setprecision(3)
              << sqrt(MonteCarlo[i].nSelEvt(i_syst, l))/MonteCarlo[i].nSelEvt(i_syst, 0)*100 << " % "
              << ")"        
              << endl;
            }
        }
      }
    }
  }
}

// #############################################################################
void LipMiniAnalysis::PostLoopCalculations(){
// #############################################################################
//
//  purpose: to compute things after the events loop (eg. asymmetries)
//
//  notes: should be replaced by the user!
//
//  authors: nfcastro, fveloso
//  first version: ??.???.????
//
//  last change: 19.nov.2006
//  by: fveloso
//
// #############################################################################

}

// #############################################################################
void LipMiniAnalysis::Loop() {
// #############################################################################

#ifdef MEASURE_READFILE
  long long totaltime = 0;
#endif

  if (nTuple->fChain == 0) return;
  Int_t nentries = Int_t(nTuple->fChain->GetEntriesFast());

  // start loop over all events
  for (Long64_t i_event = 0; nentries; ++i_event) {

#ifdef MEASURE_READFILE
    long long int time;
    timeval t;

    gettimeofday(&t, NULL);
    time = t.tv_sec * 1000000.0 + t.tv_usec;
#endif

    // standard stuff to get the event in memory
    Int_t ientry = nTuple->LoadTree(i_event);
    if (ientry < 0) break;
    nTuple->fChain->GetEntry(i_event);

    // loop over systematics
    for (Int_t i_syst=0; i_syst<Syst.size(); ++i_syst) {

      // Convert contents of memory buffer into vectors of TLorentzVectors for ease use
      RecoType=Syst[i_syst];
      FillAllVectors();
      //cout << " Loop() i_syst=" << i_syst << " RecoType=" << RecoType << endl;

      // Remember total number of events
      int mc_process=-999;
      for (int i = 1; i < MonteCarlo.size(); i++) {
        if(Isub == MonteCarlo[i].run()) mc_process=i;
      }
      if (mc_process==-999) {
        std::cout << "ERROR: MonteCarlo process " << Isub << " not defined" << endl;
        exit(0);
      }

#ifdef MEASURE_READFILE
      timeval t;
    long long int end;

    gettimeofday(&t, NULL);

    end = t.tv_sec * 1000000.0 + t.tv_usec;
    end -= init;

    totaltime += end;
#endif
      // Normalize weight to luminosity
      if(isData==0 || isQCDantie==1 || isQCDmmm==1) Weight=Weight*Luminosity/MonteCarlo[mc_process].lum();
      //else Weight=1;

//    int mc_process=1;
//    Weight=1.;

      // Do selection cuts
      LastCut = 0;
      DoCuts();

      // Check if the value of 'MaxCuts' is enough for the cuts made by 'DoCuts()'
      if (LastCut > MaxCuts - DoLike){
        cout << endl;
        cout << "Oops! The number of cuts in DoCuts() exceeded MaxCuts:" << endl;
        cout << endl;
        cout << "  MaxCuts is " << MaxCuts << " but you required " << LastCut << " cuts in DoCuts()..." << endl;
        cout << "  Try again with a higher value for MaxCuts!" << endl;
        cout << endl;
        exit(0);
      }

/* Comentado para performance; descomentar de a análise necessitar de pdf e likelihoods
      // probability density functions
      if (LastCut == MaxCuts - DoLike) {
        // fill pdf
        PdfValues();
        // compute likelihood
        if (DoLike) {
          ComputeAndFillLikelihood();
          if (LogELikeLOverLikeBValue >= LogELikeMinCut && Log10LikeLOverLikeBValue >= Log10LikeMinCut
           && LogELikeLOverLikeBValue <= LogELikeMaxCut && Log10LikeLOverLikeBValue <= Log10LikeMaxCut) LastCut++;
        }
      }

      if (LastCut == MaxCuts) {
        // fill pdf
        FillPdf(i_syst);
        // fill output ntuple
        FillOutputNtuple(i_syst);
      }
*/
      // Count events
      for(int i=0; i <= LastCut; ++i) {
        MonteCarlo[mc_process].AddSelEvt(i_syst, i);
        MonteCarlo[mc_process].AddSelWeightedEvt(i_syst, i,Weight);
        // total background
        if (MonteCarlo[mc_process].type()==1) {
          MonteCarlo[0].AddSelEvt(i_syst, i);
          MonteCarlo[0].AddSelWeightedEvt(i_syst, i, Weight);
        }
      }

      // Fill histograms
      //if (MonteCarlo[mc_process].type()==1) {
      //  FillHistograms(MonteCarlo[mc_process].histo[i_syst]);
      //}
      FillHistograms(MonteCarlo[0].histo[i_syst].histo);
      //if (MonteCarlo[mc_process].type()==1) FillHistograms(MonteCarlo[0].histo);

    // end loop over systematics
    }

  // end loop over events
  }

  #ifdef MEASURE_READFILE
  cout << "BADALHOQUICES " << totaltime << endl;
  #endif

}

// #############################################################################
THisto::THisto() {
// #############################################################################

  for (int i = 0; i < THistoX; ++i) {
    for (int j = 0; j < THistoY; ++j) {
      histo[i][j] = NULL;
    }
  }
}

// #############################################################################
TMonteCarlo::TMonteCarlo(int type, double run, double lum, int nGenEvt, std::string title, int MaxCuts) {
// #############################################################################
  if (type != 0 && type != 1 && type != 2 ) {
    std::cout << "ERROR: unknown MC type for process " << title << endl;
    exit(0);
  }
  if (lum <= 0.) {
    std::cout << "ERROR: luminosity is less than or equal to 0 for process " << title << endl;
    exit(0);
  }
  p_type=type;
  p_run=run;
  p_lum=lum;
  p_nGenEvt=nGenEvt;
  p_title=title;
  p_CrossSection=nGenEvt/lum;
/*
  for (int i=0; i <= MaxCuts; ++i) {
    p_nSelEvt.push_back(0.);
    p_nSelWeightedEvt.push_back(0.);
    p_nSelWeightedEvtErr.push_back(0.);
  }
*/
//  gDirectory->mkdir(title.c_str(), title.c_str());
//  gDirectory->mkdir("ola", "olatitulo");
  if (p_type == 1) {
    dir = new TDirectory(title.c_str(), title.c_str());
  } else if (p_type == 0) {
    dir = gDirectory;
  }

}

// #############################################################################
TMonteCarlo::~TMonteCarlo() {
// #############################################################################

}

// #############################################################################
void TMonteCarlo::SetnGenEvt(int nGenEvt) {
// #############################################################################

  if (p_type == 0 ) {
    p_nGenEvt=nGenEvt;
  } else {
    std::cout << "ERROR: Cannot set manually the statistical error for process " << p_title << endl;
    exit(0);
  }
}

// #############################################################################
void TMonteCarlo::SetnSelWeightedEvtErr(Int_t i_syst, int level, double err) {
// #############################################################################

  if (p_type == 0 ) {
    p_nSelWeightedEvtErr[i_syst][level]=err;
  } else {
    std::cout << "ERROR: Cannot set manually the statistical error for process " << p_title << endl;
    exit(0);
  }
}

// #############################################################################
double TMonteCarlo::nSelWeightedEvtErr(Int_t i_syst,int level) {
// #############################################################################

  if (p_type != 0) {
    if (nSelEvt(i_syst, level) != 0) {
      return nSelWeightedEvt(i_syst, level)/sqrt(nSelEvt(i_syst, level));
    } else {
      return nSelWeightedEvt(i_syst, 0)/nSelEvt(i_syst, 0);
    }
  } else {
    return p_nSelWeightedEvtErr[i_syst][level];
  }

}

// #############################################################################
TInput::TInput() {
// #############################################################################
}

// #############################################################################
TInput::~TInput() {
// #############################################################################
  p_Name.clear();
  p_Type.clear();
}
