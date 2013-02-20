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
//   TLorentzVectorWFlags                       \      | | / |  /    |   |  ||  
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
//   Note:      this code was built from TLorentzVectorWIdx developed by 
//              Nuno Castro   (nfcastro@lipc.fis.uc.pt)  and 
//              Filipe Veloso (fveloso@lipc.fis.uc.pt)             
//
//
// #################################################################################
//
//  file: TLorentzVectorWFlags.cxx
//
//  purpose: 	to build objects as TLorentz vectors with:
//	     	- idx(=k),  	index of particle in event listing
//	     	- isb(=kf), 	code  of particle (kf=13 for muons, etc.)
//		- IsoDeltaR, 	particle isolation 
//		- itruthMatch, 	particle truth   match information (integer=0,1)  
//		- itrigMatch, 	particle trigger match information (integer=0,1)  
//
//  authors: A.Onofre
//  first version: 25.Nov.2012
//
//  last change: 25.Nov.2012
//  by: A.Onofre
//
// #############################################################################

#define TLorentzVectorWFlags_cxx

#include "TLorentzVectorWFlags.h"

// #############################################################################
TLorentzVectorWFlags::TLorentzVectorWFlags(): TLorentzVector(), idx(-1), isb(-1), IsoDeltaR(999), itruthMatch(-1), itrigMatch(-1){
// #############################################################################

}

// #############################################################################
TLorentzVectorWFlags::TLorentzVectorWFlags(Double_t px, Double_t py, Double_t pz, Double_t E, Int_t index, Int_t index2, double p_IsoDeltaR, Int_t index3, Int_t index4) :
  TLorentzVector(px,py,pz,E), idx(index), isb(index2), IsoDeltaR(p_IsoDeltaR), itruthMatch(index3), itrigMatch(index4) {
// #############################################################################

}

// #############################################################################
TLorentzVectorWFlags::TLorentzVectorWFlags(TLorentzVector v, int index, int index2, double p_IsoDeltaR, int index3, int index4) :
  TLorentzVector(v), idx(index), isb(index2), IsoDeltaR(p_IsoDeltaR), itruthMatch(index3), itrigMatch(index4) {
// #############################################################################

}

// #############################################################################
TLorentzVectorWFlags::TLorentzVectorWFlags(const TLorentzVectorWFlags& other) :
  TLorentzVector(other), idx(other.idx), isb(other.isb), IsoDeltaR(other.IsoDeltaR), itruthMatch(other.itruthMatch), itrigMatch(other.itrigMatch) {
// #############################################################################

}

// #############################################################################
TLorentzVectorWFlags& TLorentzVectorWFlags::operator=(const TLorentzVectorWFlags& other) {
// #############################################################################

  if (&other==this) {
    return *this ;
  }
  TLorentzVector::operator=(other) ;
  idx = other.idx ;
  isb = other.isb ;
  IsoDeltaR = other.IsoDeltaR;
  itruthMatch = other.itruthMatch ;
  itrigMatch  = other.itrigMatch ;
  return *this ;
}

// #############################################################################
TLorentzVectorWFlags::~TLorentzVectorWFlags() {
// #############################################################################

}
