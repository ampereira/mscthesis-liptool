if (!isData) {

  cout <<      "\n   ATLAS LIP Mini (8TeV) ttH (dilep) Analysis" << endl;
  cout <<      "\n   Mini TopTools Monte Carlo (MiniTTHReader)"<< std::endl;
  const char *  ntu = "MiniTTHReader";


  // Initialize all Samples (make sure you only read what you want)
  Int_t Pythia8_AU2CTEQ6L1_ttH125_dilepbb = 0; 	// ttH ( mH = 125 GeV )
  Int_t TTbarLF_Alpgen_8TeV = 0;		// ttbar+LF Alpgen @ 8 TeV

  // Define Sample to Run
  // ttH ( mH = 125 GeV )
  if (Sample == 901) {cout << "Pythia8_AU2CTEQ6L1_ttH125_dilepbb" << endl;
    Pythia8_AU2CTEQ6L1_ttH125_dilepbb = 1;
  }
  // ttbar+LF Alpgen
  if (Sample == 801) {cout << "TTbarLF_Alpgen_8TeV" << endl;
    TTbarLF_Alpgen_8TeV  = 1;
  }



  // ----------------------------------------------------------------------------------------------------------------------------------------------- 
  // 								R E A D    N T U P L E S 
  // ----------------------------------------------------------------------------------------------------------------------------------------------- 
  // TMonteCarlo    Name  (1, Run Number, Luminosity, Events,            Label for Output               , MaxCuts); MonteCarlo.push_back(   Name  );
  // ----------------------------------------------------------------------------------------------------------------------------------------------- 
  // ttbar+LF Alpgen
  if ( TTbarLF_Alpgen_8TeV ) 
    Input.File(ntu, "../Links/Data_MC_MarkOwen_8TeV/MC_Alpgen_tt-LF_8TeV.root");
    double MCatNLOlum = Luminosity;   // Luminosity of sample
    TMonteCarlo mc_164440(1,     164440, MCatNLOlum,  106585, "195847_MC_Alpgen_tt-LF_8TeV", MaxCuts); MonteCarlo.push_back(mc_164440);

  // ttH ( mH = 125 GeV )
  if ( Pythia8_AU2CTEQ6L1_ttH125_dilepbb ) 
    Input.File(ntu, "../Links/mini_Pythia8_AU2CTEQ6L1_ttH125_dilepbb.root");
    MCatNLOlum = Luminosity;   // Luminosity of sample
    TMonteCarlo mc_161871(1,     161871, MCatNLOlum,   5703, "161871_Pythia8_AU2CTEQ6L1_ttH125_dilepbb", MaxCuts); MonteCarlo.push_back(mc_161871);

}
