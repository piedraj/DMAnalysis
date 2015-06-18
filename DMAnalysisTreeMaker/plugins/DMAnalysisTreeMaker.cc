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
#include "Math/GenVector/VectorUtil.h"
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
  //TH1F* h_muEta;
  //TH1F* h_muPt;

  TH1F* h_MET_all_01           ;
  TH1F* h_JetPt_all_01         ; 
  TH1F* h_LeptonPt_all_01      ;  
  TH1F* h_LeptonDeltaPhi_all_01; 

  TH1F* h_MET_all_02           ; 
  TH1F* h_JetPt_all_02         ; 
  TH1F* h_LeptonPt_all_02      ;    
  TH1F* h_LeptonDeltaPhi_all_02; 

  TH1F* h_MET_all_03           ; 
  TH1F* h_JetPt_all_03         ;  
  TH1F* h_LeptonPt_all_03      ; 
  TH1F* h_LeptonDeltaPhi_all_03;  

  TH1F* h_MET_all_04           ;
  TH1F* h_JetPt_all_04         ;             
  TH1F* h_LeptonPt_all_04      ;     
  TH1F* h_LeptonDeltaPhi_all_04; 

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
  //float          Mass; 
  float          Phi;
  float          E; 

  Bool_t operator<(const Lepton& a) const {
	return Pt > a.Pt;
  }

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

  // Electrons
  edm::Handle< std::vector<float> > handle_elCharge  ;
  edm::Handle< std::vector<float> > handle_elE       ;
  edm::Handle< std::vector<float> > handle_elEta     ;
  edm::Handle< std::vector<float> > handle_elIso03   ;
  edm::Handle< std::vector<float> > handle_elMass    ;
  edm::Handle< std::vector<float> > handle_elPhi     ;
  edm::Handle< std::vector<float> > handle_elPt      ;
  edm::Handle< std::vector<float> > handle_elisMedium;

  iEvent.getByLabel(edm::InputTag("electrons","elCharge"  ), handle_elCharge  );
  iEvent.getByLabel(edm::InputTag("electrons","elE"       ), handle_elE       );
  iEvent.getByLabel(edm::InputTag("electrons","elEta"     ), handle_elEta     );
  iEvent.getByLabel(edm::InputTag("electrons","elIso03"   ), handle_elIso03   );
  iEvent.getByLabel(edm::InputTag("electrons","elMass"    ), handle_elMass    );
  iEvent.getByLabel(edm::InputTag("electrons","elPhi"     ), handle_elPhi     );
  iEvent.getByLabel(edm::InputTag("electrons","elPt"      ), handle_elPt      );
  iEvent.getByLabel(edm::InputTag("electrons","elisMedium"), handle_elisMedium);

  if (!handle_elCharge.isValid())   return;
  if (!handle_elE.isValid())        return;
  if (!handle_elEta.isValid())      return;
  if (!handle_elIso03.isValid())    return;
  if (!handle_elMass.isValid())    return;
  if (!handle_elPhi.isValid())      return;
  if (!handle_elPt.isValid())       return;
  if (!handle_elisMedium.isValid()) return;

  const std::vector<float> elCharge   = *(handle_elCharge.product()  );
  const std::vector<float> elE        = *(handle_elE.product()       );
  const std::vector<float> elEta      = *(handle_elEta.product()     );
  const std::vector<float> elIso03    = *(handle_elIso03.product()   );
  const std::vector<float> elMass     = *(handle_elMass.product()    );
  const std::vector<float> elPhi      = *(handle_elPhi.product()     );
  const std::vector<float> elPt       = *(handle_elPt.product()      );
  const std::vector<float> elisMedium = *(handle_elisMedium.product());


  // Muons
  edm::Handle< std::vector<float> > handle_muCharge     ;
  edm::Handle< std::vector<float> > handle_muE          ;
  edm::Handle< std::vector<float> > handle_muEta        ;
  edm::Handle< std::vector<float> > handle_muIsTightMuon;
  edm::Handle< std::vector<float> > handle_muIso04      ;
  edm::Handle< std::vector<float> > handle_muMass       ;
  edm::Handle< std::vector<float> > handle_muPhi        ;
  edm::Handle< std::vector<float> > handle_muPt         ;

  iEvent.getByLabel(edm::InputTag("muons","muCharge"     ), handle_muCharge     );
  iEvent.getByLabel(edm::InputTag("muons","muE"          ), handle_muE          );
  iEvent.getByLabel(edm::InputTag("muons","muEta"        ), handle_muEta        );
  iEvent.getByLabel(edm::InputTag("muons","muIsTightMuon"), handle_muIsTightMuon);
  iEvent.getByLabel(edm::InputTag("muons","muIso04"      ), handle_muIso04      );
  iEvent.getByLabel(edm::InputTag("muons","muMass"       ), handle_muMass       );
  iEvent.getByLabel(edm::InputTag("muons","muPhi"        ), handle_muPhi        );
  iEvent.getByLabel(edm::InputTag("muons","muPt"         ), handle_muPt         );

  if (!handle_muCharge.isValid())      return;
  if (!handle_muE.isValid())           return;
  if (!handle_muEta.isValid())         return;
  if (!handle_muIsTightMuon.isValid()) return;
  if (!handle_muIso04.isValid())       return;
  if (!handle_muMass.isValid())        return;
  if (!handle_muPhi.isValid())         return;
  if (!handle_muPt.isValid())          return;

  const std::vector<float> muCharge      = *(handle_muCharge.product()     );
  const std::vector<float> muE           = *(handle_muE.product()          );
  const std::vector<float> muEta         = *(handle_muEta.product()        );
  const std::vector<float> muIsTightMuon = *(handle_muIsTightMuon.product());
  const std::vector<float> muIso04       = *(handle_muIso04.product()      );
  const std::vector<float> muMass        = *(handle_muMass.product()       );
  const std::vector<float> muPhi         = *(handle_muPhi.product()        );
  const std::vector<float> muPt          = *(handle_muPt.product()         );


  // Jets

  edm::Handle< std::vector<float> > handle_jetAK4E  ;
  edm::Handle< std::vector<float> > handle_jetAK4Eta;
  edm::Handle< std::vector<float> > handle_jetAK4Phi;
  edm::Handle< std::vector<float> > handle_jetAK4Pt ;

  iEvent.getByLabel(edm::InputTag("jetsAK4","jetAK4E"  ), handle_jetAK4E  );
  iEvent.getByLabel(edm::InputTag("jetsAK4","jetAK4Eta"), handle_jetAK4Eta);
  iEvent.getByLabel(edm::InputTag("jetsAK4","jetAK4Phi"), handle_jetAK4Phi);
  iEvent.getByLabel(edm::InputTag("jetsAK4","jetAK4Pt" ), handle_jetAK4Pt );

  if (!handle_jetAK4E.isValid())   return;
  if (!handle_jetAK4Eta.isValid()) return;
  if (!handle_jetAK4Phi.isValid()) return;
  if (!handle_jetAK4Pt.isValid())  return;

  const std::vector<float> jetAK4E   = *(handle_jetAK4E.product()  );
  const std::vector<float> jetAK4Eta = *(handle_jetAK4Eta.product());
  const std::vector<float> jetAK4Phi = *(handle_jetAK4Phi.product());
  const std::vector<float> jetAK4Pt  = *(handle_jetAK4Pt.product() );

  // MET

  edm::Handle< std::vector<float> > handle_metPt;
      
  iEvent.getByLabel(edm::InputTag("met","metPt"  ), handle_metPt);

  if (!handle_metPt.isValid()) return;

  const std::vector<float> metPt = *(handle_metPt.product()  );

//---------------------------------------------------------------------------------- 


  std::vector<Lepton> AnalysisLeptons;

  for ( UInt_t i = 0; i < elPt.size(); i++ ) {

	if ( elisMedium.at(i) != 1.0 ) continue;   // id

	if ( elIso03.at(i) > 0.11 ) continue;   // iso

	//printf("elIso03 = %f \t elPt = %f \n", elIso03.at(i), elPt.at(i));

 	Lepton elec;

	elec.flavor = 1;

	elec.charge = elCharge.at(i);

	elec.Pt     = elPt.at(i);

	elec.Eta    = elEta.at(i);

	//elec.Mass   = elMass.at(i);

	elec.Phi    = elPhi.at(i);

	elec.E      = elE.at(i);

	AnalysisLeptons.push_back(elec);

  }



  for ( UInt_t i = 0; i < muPt.size(); i++ ) {

	if ( muIsTightMuon.at(i) != 1.0 ) continue;   // id 

	if ( muIso04.at(i) > 0.12 ) continue;   // iso

 	Lepton muon;

	muon.flavor = 2;

	muon.charge = muCharge.at(i);

	muon.Pt     = muPt.at(i);

	muon.Eta    = muEta.at(i);

	//muon.Mass   = muMass.at(i);

	muon.Phi    = muPhi.at(i);

	muon.E      = muE.at(i);

	AnalysisLeptons.push_back(muon);

  }

  int n_leptons = AnalysisLeptons.size(); 
  
// Require exactly two leptons
//--------------------------------------------------------------------------
//    if (AnalysisLeptons.size() != 2) return;


// Sort
//--------------------------------------------------------------------------

//printf(" pt before: ");
//for(int i = 0; i < n_leptons; i++) printf( " %f \t", AnalysisLeptons[i].Pt );
//printf("\n pt after: ");

std::sort(AnalysisLeptons.begin(), AnalysisLeptons.end());

//for(int i = 0; i < n_leptons; i++) printf( " %f \t", AnalysisLeptons[i].Pt );
//printf("\n .................... \n");

// Lepton Pair Selection
//--------------------------------------------------------------------------
  if ( n_leptons < 2 ) return;

  //printf("number of leptons = %i \n", n_leptons); 

  bool  opposite_charge = false;	

  opposite_charge = ( AnalysisLeptons[0].charge * AnalysisLeptons[1].charge < 0 );

  if ( !opposite_charge ) return;

  //printf("pair charge = %f \n", AnalysisLeptons[0].charge * AnalysisLeptons[1].charge);
 
  Lepton lep1, lep2;

  lep1 = AnalysisLeptons[0];
  lep2 = AnalysisLeptons[1];

  math::PtEtaPhiELorentzVector lep1_tlv(lep1.Pt, lep1.Eta, lep1.Phi, lep1.E);
  math::PtEtaPhiELorentzVector lep2_tlv(lep2.Pt, lep2.Eta, lep2.Phi, lep2.E);

  //printf(" pt: %f \t eta: %f \t m: %f  \n", lep1_tlv.Pt(), lep1_tlv.Eta(), lep1_tlv.M()); 

  float m_ll = (lep1_tlv + lep2_tlv).M(); 

  if ( m_ll < 20. ) return;

  //printf("invariant mass: %f \n", m_ll );

  //printf("delta R =    %f \n", ROOT::Math::VectorUtil::DeltaR(lep1_tlv, lep2_tlv) * ROOT::Math::VectorUtil::DeltaR(lep1_tlv, lep2_tlv));
  //printf("my delta R = %f \n", (lep1.Eta-lep2.Eta)*(lep1.Eta-lep2.Eta) + (lep1.Phi-lep2.Phi)*(lep1.Phi-lep2.Phi));
  //printf("........\n");

  //int metSize = metPt.size(); 
  //printf("met vector size is: %i \n", metSize );

// Jet Selection
//-------------------------------------------------------------------------- 
  std::vector<math::PtEtaPhiELorentzVector> SelectedJets; 

  for ( UInt_t i = 0; i < jetAK4Pt.size(); i++ ) {

	if ( jetAK4Pt.at(i) < 30 ) continue; 

	if ( fabs(jetAK4Eta.at(i)) > 2.4 ) continue; 
 
        math::PtEtaPhiELorentzVector jet( jetAK4Pt.at(i), jetAK4Eta.at(i), jetAK4Phi.at(i), jetAK4E.at(i) );

        //printf(" pt: %f \t eta: %f \t phi: %f \t E: %f \t mass: %f \n", jet.Pt(), jet.Eta(), jet.Phi(), jet.E(), jet.M() ); 

	double deltaR1 = ROOT::Math::VectorUtil::DeltaR(jet, lep1_tlv);
	double deltaR2 = ROOT::Math::VectorUtil::DeltaR(jet, lep2_tlv);

	if ( deltaR1 < 0.4 || deltaR2 < 0.4 ) continue; 

	//printf( "deltaR1 = %f \t deltaR2 = %f \n", deltaR1, deltaR2 ); 

	SelectedJets.push_back(jet);
	   
	//if (jettche[i] > 2.1 && jetid[i] >= 4) nbjet++;  // jetid 4 = MVA LOOSE
	
  }   // end for 

  int n_jet = SelectedJets.size(); 


// Sort SelectedJets !!
//---------------------------------------------------------------------------------- 
// ...


  if( n_jet < 2) return;

//  printf("number of jets = %i \n", n_jet); 

  math::PtEtaPhiELorentzVector jet1_tlv = SelectedJets[0];
  math::PtEtaPhiELorentzVector jet2_tlv = SelectedJets[1];


// channel assignment
//---------------------------------------------------------------------------------- 
   int ch; 

    if ( lep1.flavor * lep2.flavor == 1 ) ch = 0;   // ee
    if ( lep1.flavor * lep2.flavor == 4 ) ch = 1;   // mumu
    if ( lep1.flavor * lep2.flavor == 2 ) ch = 2;   // emu


//---------------------------------------------------------------------------------- 

// int muPt_size = muPt.size();

// for (int i=0; i<muPt_size; i++)
//    {      
//      h_muEta->Fill(muIsTightMuon.at(i));
//      h_muPt ->Fill(muIso04.at(i));
//    }


float Zradius      = 15. ;

float lim_MET      = 320.;

float lim_jetPt    = 400.;

float lim_lepPt    = 120.;

float lim_deltaPhi = 2.  ; 


bool Zveto    = false;  if ( m_ll > 91. - Zradius   &&   m_ll < 91. + Zradius )                                                             Zveto    = true;

bool pre      = false;	if ( n_jet  < 2  ||  (   ( ch == 0 || ch == 1 )  &&  ( Zveto == true )  ) )                                         pre      = true;

bool MET      = false;  if ( metPt.at(0) < lim_MET )                                                                                        MET      = true;

bool jetPt    = false;  if ( jet1_tlv.Pt() + jet2_tlv.Pt() > lim_jetPt )                                                                    jetPt    = true;

bool lepPt    = false;  if ( lep1.Pt + lep2.Pt < lim_lepPt )                                                                                lepPt    = true; 

bool deltaPhi = false;  if ( fabs( ROOT::Math::VectorUtil::DeltaPhi(lep1_tlv, lep2_tlv) ) > lim_deltaPhi )                                  deltaPhi = true; 


if (  pre == true  ) return;

h_MET_all_01           ->Fill(metPt.at(0));
h_JetPt_all_01         ->Fill(jet1_tlv.Pt() + jet2_tlv.Pt()); 
h_LeptonPt_all_01      ->Fill(lep1.Pt + lep2.Pt);  
h_LeptonDeltaPhi_all_01->Fill(fabs( ROOT::Math::VectorUtil::DeltaPhi(lep1_tlv, lep2_tlv) )); 


if ( metPt.at(0) < 80  ) return;
 
h_MET_all_02           ->Fill(metPt.at(0)); 
h_JetPt_all_02         ->Fill(jet1_tlv.Pt() + jet2_tlv.Pt()); 
h_LeptonPt_all_02      ->Fill(lep1.Pt + lep2.Pt);    
h_LeptonDeltaPhi_all_02->Fill(fabs( ROOT::Math::VectorUtil::DeltaPhi(lep1_tlv, lep2_tlv) )); 


if (  jetPt == true  ||  lepPt == true  ||  deltaPhi == true  ) return; 

h_MET_all_03           ->Fill(metPt.at(0)); 
h_JetPt_all_03         ->Fill(jet1_tlv.Pt() + jet2_tlv.Pt());  
h_LeptonPt_all_03      ->Fill(lep1.Pt + lep2.Pt); 
h_LeptonDeltaPhi_all_03->Fill(fabs( ROOT::Math::VectorUtil::DeltaPhi(lep1_tlv, lep2_tlv) ));  


if (  MET == true  ) return; 

h_MET_all_04           ->Fill(metPt.at(0));
h_JetPt_all_04         ->Fill(jet1_tlv.Pt() + jet2_tlv.Pt());             
h_LeptonPt_all_04      ->Fill(lep1.Pt + lep2.Pt);     
h_LeptonDeltaPhi_all_04->Fill(fabs( ROOT::Math::VectorUtil::DeltaPhi(lep1_tlv, lep2_tlv) )); 


}   // end DMAnalysisTreeMaker::analyze




void DMAnalysisTreeMaker::beginJob()
{
  edm::Service<TFileService> fs;

  //h_muEta = fs->make<TH1F>("h_muEta", ";muon #eta",        100, -3,   3);
  //h_muPt  = fs->make<TH1F>("h_muPt",  ";muon p_{T} [GeV]", 100,  0, 300);

//histo definition at Mkhistogr2.c:  h[i][l][m] = new TH1F( id[i][l][m],    ntu[i][1],    blup[i][0], blup[i][1], blup[i][2] )     with  id[i][l][m] = ntu[i][0] + "_" + namech[l] + "_" + namecut[m]

 h_MET_all_01            = fs->make<TH1F>( "MET_all_01"           ,  "E_{T}^{miss} [GeV]"               ,  15, 0 , 600  );
 h_JetPt_all_01          = fs->make<TH1F>( "JetPt_all_01"         ,  "sum of p_{T} of two jets [GeV]"   ,  12, 0 , 1200 );    
 h_LeptonPt_all_01       = fs->make<TH1F>( "LeptonPt_all_01"      ,  "sum of p_{T} of two leptons [GeV]",  14, 40, 600  );   
 h_LeptonDeltaPhi_all_01 = fs->make<TH1F>( "LeptonDeltaPhi_all_01",  "#Delta#phi_{ll}"                  ,  16, 0 , 3.2  ); 

 h_MET_all_02            = fs->make<TH1F>( "MET_all_02"           ,  "E_{T}^{miss} [GeV]"               ,  15, 0 , 600  );
 h_JetPt_all_02          = fs->make<TH1F>( "JetPt_all_02"         ,  "sum of p_{T} of two jets [GeV]"   ,  12, 0 , 1200 );    
 h_LeptonPt_all_02       = fs->make<TH1F>( "LeptonPt_all_02"      ,  "sum of p_{T} of two leptons [GeV]",  14, 40, 600  );   
 h_LeptonDeltaPhi_all_02 = fs->make<TH1F>( "LeptonDeltaPhi_all_02",  "#Delta#phi_{ll}"                  ,  16, 0 , 3.2  ); 

 h_MET_all_03            = fs->make<TH1F>( "MET_all_03"           ,  "E_{T}^{miss} [GeV]"               ,  15, 0 , 600  );
 h_JetPt_all_03          = fs->make<TH1F>( "JetPt_all_03"         ,  "sum of p_{T} of two jets [GeV]"   ,  12, 0 , 1200 );    
 h_LeptonPt_all_03       = fs->make<TH1F>( "LeptonPt_all_03"      ,  "sum of p_{T} of two leptons [GeV]",  14, 40, 600  );   
 h_LeptonDeltaPhi_all_03 = fs->make<TH1F>( "LeptonDeltaPhi_all_03",  "#Delta#phi_{ll}"                  ,  16, 0 , 3.2  );   

 h_MET_all_04            = fs->make<TH1F>( "MET_all_04"           ,  "E_{T}^{miss} [GeV]"               ,  15, 0 , 600  );
 h_JetPt_all_04          = fs->make<TH1F>( "JetPt_all_04"         ,  "sum of p_{T} of two jets [GeV]"   ,  12, 0 , 1200 );    
 h_LeptonPt_all_04       = fs->make<TH1F>( "LeptonPt_all_04"      ,  "sum of p_{T} of two leptons [GeV]",  14, 40, 600  );   
 h_LeptonDeltaPhi_all_04 = fs->make<TH1F>( "LeptonDeltaPhi_all_04",  "#Delta#phi_{ll}"                  ,  16, 0 , 3.2  ); 

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
