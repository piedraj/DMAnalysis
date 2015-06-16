// -*- C++ -*-
//
// Package:    DMAnalysis/DMAnalysisTreeMaker
// Class:      DMAnalysisTreeMaker
// 
/**\class DMAnalysisTreeMaker DMAnalysisTreeMaker.cc DMAnalysis/DMAnalysisTreeMaker/plugins/DMAnalysisTreeMaker.cc

 Description: This class reads B2G EDM files and runs a ttbar + DM search in the dilepton channel

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Jonatan Piedra Gomez
//         Created:  Mon, 15 Jun 2015 11:47:44 GMT
//
//


// System include files
#include <memory>

// User include files
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "TH1F.h"
#include "DataFormats/Math/interface/LorentzVector.h"
////////////////////////#include "TLorentzVector.h"


class DMAnalysisTreeMaker : public edm::EDAnalyzer {
public:
  explicit DMAnalysisTreeMaker(const edm::ParameterSet&);
  ~DMAnalysisTreeMaker();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  

  // Data members
  TH1F* h_muEta;
  TH1F* h_muPt;
};


//
// constants, enums and typedefs
//

struct Lepton
{
  UInt_t         flavor;  // Muon, Electron
  float          charge;
  float          Pt;
  float		 Eta;
  float          Phi;
  float          E; 

  //Bool_t operator<(const Lepton& a) const
  //{
  //  return a.Pt() < a.Pt();
  //}
};

//
// static data member definitions
//


DMAnalysisTreeMaker::DMAnalysisTreeMaker(const edm::ParameterSet& iConfig)
{
}


DMAnalysisTreeMaker::~DMAnalysisTreeMaker()
{
}


void DMAnalysisTreeMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  printf(" Hello 00\n");

  // Electrons
  edm::Handle< std::vector<float> > handle_elCharge  ;
  edm::Handle< std::vector<float> > handle_elE       ;
  edm::Handle< std::vector<float> > handle_elEta     ;
  edm::Handle< std::vector<float> > handle_elIso03   ;
  edm::Handle< std::vector<float> > handle_elPhi     ;
  edm::Handle< std::vector<float> > handle_elPt      ;
  edm::Handle< std::vector<float> > handle_elisMedium;

  iEvent.getByLabel(edm::InputTag("electrons","elCharge"  ), handle_elCharge  );
  iEvent.getByLabel(edm::InputTag("electrons","elE"       ), handle_elE       );
  iEvent.getByLabel(edm::InputTag("electrons","elEta"     ), handle_elEta     );
  iEvent.getByLabel(edm::InputTag("electrons","elIso03"   ), handle_elIso03   );
  iEvent.getByLabel(edm::InputTag("electrons","elPhi"     ), handle_elPhi     );
  iEvent.getByLabel(edm::InputTag("electrons","elPt"      ), handle_elPt      );
  iEvent.getByLabel(edm::InputTag("electrons","elisMedium"), handle_elisMedium);

  if (!handle_elCharge.isValid())   return;
  if (!handle_elE.isValid())        return;
  if (!handle_elEta.isValid())      return;
  if (!handle_elIso03.isValid())    return;
  if (!handle_elPhi.isValid())      return;
  if (!handle_elPt.isValid())       return;
  if (!handle_elisMedium.isValid()) return;

  const std::vector<float> elCharge   = *(handle_elCharge.product()  );
  const std::vector<float> elE        = *(handle_elE.product()       );
  const std::vector<float> elEta      = *(handle_elEta.product()     );
  const std::vector<float> elIso03    = *(handle_elIso03.product()   );
  const std::vector<float> elPhi      = *(handle_elPhi.product()     );
  const std::vector<float> elPt       = *(handle_elPt.product()      );
  const std::vector<float> elisMedium = *(handle_elisMedium.product());


  // Muons
  edm::Handle< std::vector<float> > handle_muCharge     ;
  edm::Handle< std::vector<float> > handle_muE          ;
  edm::Handle< std::vector<float> > handle_muEta        ;
  edm::Handle< std::vector<float> > handle_muIsTightMuon;
  edm::Handle< std::vector<float> > handle_muIso04      ;
  edm::Handle< std::vector<float> > handle_muPhi        ;
  edm::Handle< std::vector<float> > handle_muPt         ;

  iEvent.getByLabel(edm::InputTag("muons","muCharge"     ), handle_muCharge     );
  iEvent.getByLabel(edm::InputTag("muons","muE"          ), handle_muE          );
  iEvent.getByLabel(edm::InputTag("muons","muEta"        ), handle_muEta        );
  iEvent.getByLabel(edm::InputTag("muons","muIsTightMuon"), handle_muIsTightMuon);
  iEvent.getByLabel(edm::InputTag("muons","muIso04"      ), handle_muIso04      );
  iEvent.getByLabel(edm::InputTag("muons","muPhi"        ), handle_muPhi        );
  iEvent.getByLabel(edm::InputTag("muons","muPt"         ), handle_muPt         );

  if (!handle_muCharge.isValid())      return;
  if (!handle_muE.isValid())           return;
  if (!handle_muEta.isValid())         return;
  if (!handle_muIsTightMuon.isValid()) return;
  if (!handle_muIso04.isValid())       return;
  if (!handle_muPhi.isValid())         return;
  if (!handle_muPt.isValid())          return;

  const std::vector<float> muCharge      = *(handle_muCharge.product()     );
  const std::vector<float> muE           = *(handle_muE.product()          );
  const std::vector<float> muEta         = *(handle_muEta.product()        );
  const std::vector<float> muIsTightMuon = *(handle_muIsTightMuon.product());
  const std::vector<float> muIso04       = *(handle_muIso04.product()      );
  const std::vector<float> muPhi         = *(handle_muPhi.product()        );
  const std::vector<float> muPt          = *(handle_muPt.product()         );


  // Jets

  edm::Handle< std::vector<float> > handle_jetAK4E  ;
  edm::Handle< std::vector<float> > handle_jetAK4Eta;
  edm::Handle< std::vector<float> > handle_jetAK4Phi;
  edm::Handle< std::vector<float> > handle_jetAK4Pt ;

  iEvent.getByLabel(edm::InputTag("jetAK4","jetAK4E"  ), handle_jetAK4E  );
  iEvent.getByLabel(edm::InputTag("jetAK4","jetAK4Eta"), handle_jetAK4Eta);
  iEvent.getByLabel(edm::InputTag("jetAK4","jetAK4Phi"), handle_jetAK4Phi);
  iEvent.getByLabel(edm::InputTag("jetAK4","jetAK4Pt" ), handle_jetAK4Pt );

  if (!handle_jetAK4E.isValid())   return;
  if (!handle_jetAK4Eta.isValid()) return;
  if (!handle_jetAK4Phi.isValid()) return;
  if (!handle_jetAK4Pt.isValid())  return;

  const std::vector<float> jetAK4E   = *(handle_jetAK4E.product()  );
  const std::vector<float> jetAK4Eta = *(handle_jetAK4Eta.product());
  const std::vector<float> jetAK4Phi = *(handle_jetAK4Phi.product());
  const std::vector<float> jetAK4Pt  = *(handle_jetAK4Pt.product() );


  printf(" Hello 01\n");


//---------------------------------------------------------------------------------- 


  std::vector<Lepton> AnalysisLeptons;
  //  std::vector<TLorentzVector> SelectedJets;

  for ( UInt_t i = 0; i < elPt.size(); i++ ) {

	if ( 1 == 2 ) continue;   // id

	if ( 1 == 2 ) continue;   // iso

 	Lepton elec;

	elec.flavor = 1;

	elec.charge = elCharge.at(i);

	elec.Pt     = elPt.at(i);

	elec.Eta    = elEta.at(i);

	elec.Phi    = elPhi.at(i);

	elec.E      = elE.at(i);

	AnalysisLeptons.push_back(elec);

  }


  printf(" Hello 02\n");


  for ( UInt_t i = 0; i < muPt.size(); i++ ) {

	if ( 1 == 2 ) continue;   // id 

	if ( 1 == 2 ) continue;   // iso

 	Lepton muon;

	muon.flavor = 2;

	muon.charge = muCharge.at(i);

	muon.Pt     = muPt.at(i);

	muon.Eta    = muEta.at(i);

	muon.Phi    = muPhi.at(i);

	muon.E      = muE.at(i);

	AnalysisLeptons.push_back(muon);

  }
  printf(" Hello 03\n");


// Require exactly two leptons
//--------------------------------------------------------------------------
//    if (AnalysisLeptons.size() != 2) continue;


// Sort
//--------------------------------------------------------------------------
// ... 


// Lepton Pair Selection
//--------------------------------------------------------------------------
//  if ( AnalysisLeptons.size() < 2 ) continue;

//  if ( AnalysisLeptons[0].charge * AnalysisLeptons[1].charge < 0 ) continue;

  Lepton lep1, lep2;

  lep1 = AnalysisLeptons[0];

  lep2 = AnalysisLeptons[1];

  math::PtEtaPhiELorentzVector lep1_tlv(lep1.Pt, lep1.Eta, lep1.Phi, lep1.E);
    //lep2_tlv; 

  //  lep1_tlv.SetPtEtaPhiE

  printf(" pt: %f \t eta: %f \t m: %f\n", lep1_tlv.Pt(), lep1_tlv.Eta(), lep1_tlv.M());

  //  lep2_tlv.SetPtEtaPhiE(lep2.Pt, lep2.Eta, lep2.Phi, lep2.E);

// if ( (lep1_tlv+lep2_tlv).M() > 20. ) continue;


// Jet Selection
//-------------------------------------------------------------------------- 
/*  for ( UInt_t i = 0; i < jetAK4Pt.size(); i++ ) {

	if ( fabs(jetAK4Eta.at(i)) >= 2.4 ) continue;

//	TLorentzVector Jet;
      
//      	Jet.SetPtEtaPhiE(jetAK4Pt.at(i), jetAK4Eta.at(i), jetAK4Phi.at(i), jetAK4E.at(i));

      	Bool_t thisJetIsLepton = false;
      
      	for ( UInt_t j = 0; j < AnalysisLeptons.size(); j++ ) {

//		if ( Jet.DeltaR(lep1_tlv) <= 0.4 && Jet.DeltaR(lep2_tlv) ) {

		      thisJetIsLepton = true;
	      
		      break;

		}

	}

      	if (thisJetIsLepton) continue;

	if (jetAK4Pt[i] > 30.) {

	//		SelectedJets.push_back(Jet);
	  
	  	//njet++;
	  
	  	//if (jettche[i] > 2.1 && jetid[i] >= 4) nbjet++;  // jetid 4 = MVA LOOSE
	
	}

}*/


// channel assignment
//---------------------------------------------------------------------------------- 
//  int ch; 
//
//  if ( lep1.flavor * lep2.flavor == 1 ) ch = 0;   // ee
//  if ( lep1.flavor * lep2.flavor == 4 ) ch = 1;   // mumu
//  if ( lep1.flavor * lep2.flavor == 2 ) ch = 2;   // emu


//---------------------------------------------------------------------------------- 

// int muPt_size = muPt.size();


    /*for (int i=0; i<muPt_size; i++)
    {      
      h_muEta->Fill(muIsTightMuon.at(i));
      h_muPt ->Fill(muIso04.at(i));
    }*/


}   // end DMAnalysisTreeMaker::analyze


void DMAnalysisTreeMaker::beginJob()
{
  edm::Service<TFileService> fs;

  h_muEta = fs->make<TH1F>("h_muEta", ";muon #eta",        100, -3,   3);
  h_muPt  = fs->make<TH1F>("h_muPt",  ";muon p_{T} [GeV]", 100,  0, 300);

}


void DMAnalysisTreeMaker::endJob()
{
}


// Fill 'descriptions' with the allowed parameters for the module
void DMAnalysisTreeMaker::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
  // The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


// Define this as a plug-in
DEFINE_FWK_MODULE(DMAnalysisTreeMaker);
