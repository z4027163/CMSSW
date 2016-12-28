#ifndef MELA_SpinTwoqqMinimalMela_h
#define MELA_SpinTwoqqMinimalMela_h

/** \class SpinTwoqqMinimalMela
 *
 *  SpinTwoqqMinimalMELA discriminator 
 *
 *
 *  $Date: 2013/01/10 02:43:12 $
 *  $Revision: 1.1 $
 *  \author JHU
 */

#include <TLorentzVector.h>

class AngularPdfFactory;
class TensorPdfFactory;
class RooRealVar;


class SpinTwoqqMinimalMELA{

public:

  SpinTwoqqMinimalMELA();

  ~SpinTwoqqMinimalMELA();

  void computeKD(TLorentzVector Z1_lept1, int Z1_lept1Id,
		 TLorentzVector Z1_lept2, int Z1_lept2Id,
		 TLorentzVector Z2_lept1, int Z2_lept1Id,
		 TLorentzVector Z2_lept2, int Z2_lept2Id,
		 float& kd, 
		 float& psig,
		 float& pbkg);
  
  void computeKD(float zzmass, float z1mass, float z2mass, 
		 float costhetstar, 
		 float costheta1, 
		 float costheta2, 
		 float phi, 
		 float phistar1,
		 float& kd, 
		 float& psig, 
		 float& psigALT);

private:
  void checkZorder(float& z1mass, float& z2mass,
		   float& costhetastar, float& costheta1,
		   float& costheta2, float& phi,
		   float& phistar1);

  RooRealVar* z1mass_rrv;
  RooRealVar* z2mass_rrv;
  RooRealVar* costheta1_rrv;
  RooRealVar* costheta2_rrv;
  RooRealVar* phi_rrv;
  RooRealVar* costhetastar_rrv;
  RooRealVar* phistar1_rrv;
  RooRealVar* mzz_rrv;

  AngularPdfFactory *SMHiggs;
  TensorPdfFactory *minGrav;   

};

#endif
