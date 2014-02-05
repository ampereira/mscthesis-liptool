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
//  purpose:    to build objects as TLorentz vectors with:
//              - idx(=k),      index of particle in event listing
//              - isb(=kf),     code  of particle (kf=13 for muons, etc.)
//              - IsoDeltaR,    particle isolation 
//              - itruthMatch,  particle truth   match information (integer=0,1)  
//              - itrigMatch,   particle trigger match information (integer=0,1)  
//
//  authors: A.Onofre
//  first version: 25.Nov.2012
//
//  last change: 25.Nov.2012
//  by: A.Onofre
//
// #############################################################################    

#ifndef TLorentzVectorWFlags_h
#define TLorentzVectorWFlags_h

#include <TLorentzVector.h>
#include <vector>

// #############################################################################
class TLorentzVectorWFlags : public TLorentzVector {
// #############################################################################
// 
//  purpose: TlorentzVector with additional relevant info
//
//  authors: A.Onofre
//  first version: 25.Nov.2012
//
//  last change: 25.Nov.2012
//  by: A.Onofre
//
// #############################################################################

 public:
  // 
  TLorentzVectorWFlags() ;
  TLorentzVectorWFlags(Double_t px, Double_t py, Double_t pz, Double_t E, int idx, int isb, double IsoDeltaR, int itruthMatch, int itrigMatch);
  TLorentzVectorWFlags(TLorentzVector v, int idx, int isb, double IsoDeltaR, int itruthMatch, int itrigMatch);
  TLorentzVectorWFlags(const TLorentzVectorWFlags& other);
  //
  ~TLorentzVectorWFlags();
  TLorentzVectorWFlags& operator=(const TLorentzVectorWFlags& other) ;

  inline void SetIsoDeltaR(double p_IsoDeltaR) {IsoDeltaR = p_IsoDeltaR;};

  int idx ;
  int isb ;
  double IsoDeltaR;
  int itruthMatch; 
  int itrigMatch;

 private:

  //  ClassDef(TLorentzVectorWFlags,0)
} ;

#endif
