if (isData) {

  cout <<      "\n   ATLAS LIP Mini (8TeV) ttH (dilep) Analysis" << endl;
  cout <<      "\n   Mini TopTools Data (MiniTTHReader)"<< std::endl;
  const char *  ntu = "MiniTTHReader";

  if(alldata_ele){
	// input data root files from Egamma stream
        cout << " Define Samples Selected:    define_samples_data_MarkOwen.cxx " << endl;
	Input.File(ntu, "../Links/Data_MC_MarkOwen_8TeV/Data_Egamma_8TeV.root");
	// input TMonteCarlo file from Egamma stream
        cout << "               Runs List:    ../Links/Data_MC_MarkOwen_8TeV/Data_Egamma_8TeV_Make-TMonteCarlo.txt" << endl;
	#include "../Links/Data_MC_MarkOwen_8TeV/Data_Egamma_8TeV_Make-TMonteCarlo.txt"
  }


  if(alldata_muo){
	// input data root files from Egamma stream
        cout << " Define Samples Selected:    define_samples_data_MarkOwen.cxx " << endl;
	Input.File(ntu, "../Links/Data_MC_MarkOwen_8TeV/Data_Muons_8TeV.root");
	// input TMonteCarlo file from Egamma stream
        cout << "               Runs List:    ../Links/Data_MC_MarkOwen_8TeV/Data_Muons_8TeV_Make-TMonteCarlo.txt" << endl;
	#include "../Links/Data_MC_MarkOwen_8TeV/Data_Muons_8TeV_Make-TMonteCarlo.txt"
  }

}
