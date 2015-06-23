//////Takes in the values of track quantites and returns an 
//////L1TReginalMuonCandidate with proper values and legths
//////
//////
//////
//////


#include "DataFormats/L1TMuon/interface/L1TRegionalMuonCandidate.h"
#include "DataFormats/L1TMuon/interface/L1TRegionalMuonCandidateFwd.h"



int GetPackedEta(float theta, int sector){

	float scale = 1/0.010875;
	
	float theta_angle = (theta*0.851562 + 8.5)*(3.14159265359/180);
	float eta = (-1)*log(tan(theta_angle/2));
	if(sector > 5)
		eta *= -1;

	int PackedEta = eta*scale;
	if(eta < 0)
		PackedEta -= 1;
		
	if(PackedEta > 239)
		PackedEta = 239;
	
	if(PackedEta < -240)
		PackedEta = -240;


	return PackedEta;

}


l1t::L1TRegionalMuonCandidate MakeRegionalCand(float pt, int phi, int theta, 
											   int sign, int quality, 
											   int trackaddress, int sector){
											   
	l1t::L1TRegionalMuonCandidate Cand;									   

	int iEta = GetPackedEta(theta,sector);
	
	l1t::tftype TFtype = l1t::tftype::emtf_neg;
	if(sector > 5){
		TFtype = l1t::tftype::emtf_pos;
		sector -= 6;
	}
		
	int iPt = pt*2;
	if(iPt > 511)
		iPt = 511;
	
	if(iPt < 0)
		iPt = 0;
		
	int iQual = quality/8;
		
	Cand.setHwPt(iPt);
	Cand.setHwEta(iEta);
  	Cand.setHwPhi(phi/4);//this is relative phi not global. Needs to be decided on still
  	Cand.setHwSign(1);
	Cand.setHwSignValid(0);
  	Cand.setHwQual(iQual);
  	Cand.setHwTrackAddress(1);
	Cand.setTFIdentifiers(sector,TFtype);
  	

	return Cand;

}
