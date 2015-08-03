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
#include "TTree.h"
#include "DataFormats/Math/interface/LorentzVector.h" 
#include "Math/GenVector/VectorUtil.h"


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

 bool readGen_;

//-----   branches   ----------------------------------------

  TTree* DMTree = new TTree();

  float t_metPt;

  // leptons
  float t_ch; 

  float t_lep1Pt ;
  float t_lep1Eta;
  float t_lep1Phi;
  float t_lep1E  ;

  float t_lep2Pt ;
  float t_lep2Eta;
  float t_lep2Phi;
  float t_lep2E  ;

  float t_dilepInvMass;

  // jets
  float t_jet1Pt ;
  float t_jet1Eta;
  float t_jet1Phi;
  float t_jet1E  ;

  float t_jet2Pt ;
  float t_jet2Eta;
  float t_jet2Phi;
  float t_jet2E  ;

  // event shape definitions
  float t_C_a          ;
  float t_C_b          ;
  float t_D_a          ;
  float t_D_b          ;
  float t_aplanarity_a ;
  float t_aplanarity_b ;
  float t_circularity_a;
  float t_circularity_b;
  float t_isotropy_a   ;
  float t_isotropy_b   ;
  float t_sphericity_a ;
  float t_sphericity_b ;
  float t_thrust_a     ;
  float t_thrust_b     ;
  float t_centrality   ;

  //genPart
  float t_genTopPt     ;	
  float t_genTopPhi    ;	
  float t_genTopEta    ;	
  float t_genTopE      ;
	
  float t_genAntiTopPt ;	
  float t_genAntiTopPhi;	
  float t_genAntiTopEta;	
  float t_genAntiTopE  ;

  /*std::vector<float> *t_genPartCharge;   
  std::vector<float> *t_genPartE     ;
  std::vector<float> *t_genPartEta   ;
  std::vector<float> *t_genPartId    ;
  std::vector<float> *t_genPartMass  ;
  std::vector<float> *t_genPartMomID ;
  std::vector<float> *t_genPartPhi   ;
  std::vector<float> *t_genPartPt    ;
  std::vector<float> *t_genPartStatus;
  std::vector<float> *t_genPartY     ;*/

  //trigger
  std::vector<float>       *t_triggerBitTree;
  std::vector<int>         *t_triggerPrescaleTree;
  std::vector<std::string> *t_triggerNameTree;

//-----   histos   ------------------------------------------

/*  TH1F* h_MET_all_01           ;
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
  TH1F* h_LeptonDeltaPhi_all_04;*/
};


//
// constants, enums and typedefs
//
enum {Electron=1, Muon=2};

struct Lepton
{
  UInt_t Flavor;
  float  Charge;
  float  Pt;
  float	 Eta;
  float  Mass; 
  float  Phi;
  float  E; 

  Bool_t operator<(const Lepton& a) const {
    return Pt > a.Pt;
  }
};

bool myfunction (math::PtEtaPhiELorentzVector v1, math::PtEtaPhiELorentzVector v2) { 
	return ( v1.Pt() > v2.Pt() ); 
} 


//
// static data member definitions
//


DMAnalysisTreeMaker::DMAnalysisTreeMaker(const edm::ParameterSet& iConfig) :
readGen_ (iConfig.getUntrackedParameter<bool>("readGen", false))
{

}


DMAnalysisTreeMaker::~DMAnalysisTreeMaker()
{
}


void DMAnalysisTreeMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{



  //-- modelo -----
  // edm::Handle< <column1> > handle_metPt; 
  // iEvent.getByLabel(edm::InputTag("<column2>","<column3>"), handle_metPt);
  // if (!handle_metPt.isValid()) return;
  // const std::vector<float> metPt = *(handle_metPt.product());

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
  if (!handle_elMass.isValid())     return;
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
      
  iEvent.getByLabel(edm::InputTag("met","metPt"), handle_metPt);

  if (!handle_metPt.isValid()) return;

  const std::vector<float> metPt = *(handle_metPt.product());
  

  // event shape definitions
  edm::Handle< double > handle_C_a          ;
  edm::Handle< double > handle_C_b          ;
  edm::Handle< double > handle_D_a          ;
  edm::Handle< double > handle_D_b          ;
  edm::Handle< double > handle_aplanarity_a ;
  edm::Handle< double > handle_aplanarity_b ;
  edm::Handle< double > handle_circularity_a;
  edm::Handle< double > handle_circularity_b;
  edm::Handle< double > handle_isotropy_a   ;
  edm::Handle< double > handle_isotropy_b   ;
  edm::Handle< double > handle_sphericity_a ;
  edm::Handle< double > handle_sphericity_b ;
  edm::Handle< double > handle_thrust_a     ;
  edm::Handle< double > handle_thrust_b     ;
  edm::Handle< double > handle_centrality   ;
 
  iEvent.getByLabel(edm::InputTag("eventShapePFJetVars", "C"          ), handle_C_a          );
  iEvent.getByLabel(edm::InputTag("eventShapePFVars"   , "C"          ), handle_C_b          );
  iEvent.getByLabel(edm::InputTag("eventShapePFJetVars", "D"          ), handle_D_a          );
  iEvent.getByLabel(edm::InputTag("eventShapePFVars"   , "D"          ), handle_D_b          );
  iEvent.getByLabel(edm::InputTag("eventShapePFJetVars", "aplanarity" ), handle_aplanarity_a );
  iEvent.getByLabel(edm::InputTag("eventShapePFVars"   , "aplanarity" ), handle_aplanarity_b );
  iEvent.getByLabel(edm::InputTag("eventShapePFJetVars", "circularity"), handle_circularity_a);
  iEvent.getByLabel(edm::InputTag("eventShapePFVars"   , "circularity"), handle_circularity_b);
  iEvent.getByLabel(edm::InputTag("eventShapePFJetVars", "isotropy"   ), handle_isotropy_a   );
  iEvent.getByLabel(edm::InputTag("eventShapePFVars"   , "isotropy"   ), handle_isotropy_b   );
  iEvent.getByLabel(edm::InputTag("eventShapePFJetVars", "sphericity" ), handle_sphericity_a );
  iEvent.getByLabel(edm::InputTag("eventShapePFVars"   , "sphericity" ), handle_sphericity_b );
  iEvent.getByLabel(edm::InputTag("eventShapePFJetVars", "thrust"     ), handle_thrust_a     );
  iEvent.getByLabel(edm::InputTag("eventShapePFVars"   , "thrust"     ), handle_thrust_b     );
  iEvent.getByLabel(edm::InputTag("centrality"         , "centrality" ), handle_centrality   );

  if (!handle_C_a.isValid())           return;
  if (!handle_C_b.isValid())           return;
  if (!handle_D_a.isValid())           return;
  if (!handle_D_b.isValid())           return;
  if (!handle_aplanarity_a.isValid())  return;
  if (!handle_aplanarity_b.isValid())  return;
  if (!handle_circularity_a.isValid()) return;
  if (!handle_circularity_b.isValid()) return;
  if (!handle_isotropy_a.isValid())    return;
  if (!handle_isotropy_b.isValid())    return;
  if (!handle_sphericity_a.isValid())  return;
  if (!handle_sphericity_b.isValid())  return;
  if (!handle_thrust_a.isValid())      return;
  if (!handle_thrust_b.isValid())      return;
  if (!handle_centrality.isValid())    return;

  const double C_a           = *(handle_C_a.product()          );
  const double C_b           = *(handle_C_b.product()          );
  const double D_a           = *(handle_D_a.product()          );
  const double D_b           = *(handle_D_b.product()          );
  const double aplanarity_a  = *(handle_aplanarity_a.product() );
  const double aplanarity_b  = *(handle_aplanarity_b.product() );
  const double circularity_a = *(handle_circularity_a.product());
  const double circularity_b = *(handle_circularity_b.product());
  const double isotropy_a    = *(handle_isotropy_a.product()   );
  const double isotropy_b    = *(handle_isotropy_b.product()   );
  const double sphericity_a  = *(handle_sphericity_a.product() );
  const double sphericity_b  = *(handle_sphericity_b.product() );
  const double thrust_a      = *(handle_thrust_a.product()     );
  const double thrust_b      = *(handle_thrust_b.product()     );
  const double centrality    = *(handle_centrality.product()   );


  //trigger
  edm::Handle< std::vector<float> >       handle_triggerBitTree     ;
  edm::Handle< std::vector<int> >         handle_triggerPrescaleTree;
  edm::Handle< std::vector<std::string> > handle_triggerNameTree    ;

  iEvent.getByLabel(edm::InputTag("TriggerUserData","triggerBitTree"     ), handle_triggerBitTree     );
  iEvent.getByLabel(edm::InputTag("TriggerUserData","triggerPrescaleTree"), handle_triggerPrescaleTree);
  iEvent.getByLabel(edm::InputTag("TriggerUserData","triggerNameTree"    ), handle_triggerNameTree    );

  if (!handle_triggerBitTree.isValid())      return;
  if (!handle_triggerPrescaleTree.isValid()) return;
  if (!handle_triggerNameTree.isValid())     return;

  const std::vector<float>       triggerBitTree      = *(handle_triggerBitTree.product()     );
  const std::vector<int>         triggerPrescaleTree = *(handle_triggerPrescaleTree.product());
  const std::vector<std::string> triggerNameTree     = *(handle_triggerNameTree.product()    );



  // Tree variables   // vectorial branch

  t_triggerBitTree      = new std::vector<float>      ;
  t_triggerPrescaleTree = new std::vector<int>        ;
  t_triggerNameTree     = new std::vector<std::string>;



//---------------------------------------------------------------------------------- 


  std::vector<Lepton> AnalysisLeptons;

  for ( UInt_t i = 0; i < elPt.size(); i++ ) {

	if ( elisMedium.at(i) != 1.0 ) continue;   // id

	if ( elIso03.at(i) > 0.11 ) continue;   // iso

	//printf("elIso03 = %f \t elPt = %f \n", elIso03.at(i), elPt.at(i));

 	Lepton elec;

	elec.Flavor = Electron;

	elec.Charge = elCharge.at(i);

	elec.Pt     = elPt.at(i);

	elec.Eta    = elEta.at(i);

	//elec.Mass   = elMass.at(i);

	elec.Phi    = elPhi.at(i);

	elec.E      = elE.at(i);

	AnalysisLeptons.push_back(elec);

  }



  for ( UInt_t i = 0; i < muPt.size(); i++ ) {

    //t_muPt->push_back(muPt.at(i));    // vectorial branch

	if ( muIsTightMuon.at(i) != 1.0 ) continue;   // id 

	if ( muIso04.at(i) > 0.12 ) continue;   // iso

 	Lepton muon;

	muon.Flavor = 2;

	muon.Charge = muCharge.at(i);

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
  if ( n_leptons < 2  ) return;

  bool opposite_charge = false;	

  opposite_charge = ( AnalysisLeptons[0].Charge * AnalysisLeptons[1].Charge < 0 );

  if ( !opposite_charge ) return;

  Lepton lep1, lep2;

  lep1 = AnalysisLeptons[0];
  lep2 = AnalysisLeptons[1];

  math::PtEtaPhiELorentzVector lep1_tlv(lep1.Pt, lep1.Eta, lep1.Phi, lep1.E);
  math::PtEtaPhiELorentzVector lep2_tlv(lep2.Pt, lep2.Eta, lep2.Phi, lep2.E);

  //printf(" pt: %f \t eta: %f \t m: %f  \n", lep1_tlv.Pt(), lep1_tlv.Eta(), lep1_tlv.M()); 

  float m_ll = (lep1_tlv + lep2_tlv).M(); 

  if ( m_ll < 20. ) return; 

  printf(" -- 3 \n");

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
  //printf(" pt before: ");
  //for(int i = 0; i < n_jet; i++) printf( " %f \t", SelectedJets[i].Pt() );
  //printf("\n pt after: ");

  std::sort(SelectedJets.begin(), SelectedJets.end(), myfunction);

  //for(int i = 0; i < n_jet; i++) printf( " %f \t", SelectedJets[i].Pt() );
  // printf("\n .................... \n");


  if( n_jet < 2) return;

  //  printf("number of jets = %i \n", n_jet); 

  math::PtEtaPhiELorentzVector jet1_tlv = SelectedJets[0];   // comented at 8.vii
  math::PtEtaPhiELorentzVector jet2_tlv = SelectedJets[1];   // comented at 8.vii




// top generation variables
// ------------------------

  // genPart
  edm::Handle< std::vector<float> > handle_genPartCharge;      
  edm::Handle< std::vector<float> > handle_genPartE     ;       
  edm::Handle< std::vector<float> > handle_genPartEta   ;      
  edm::Handle< std::vector<float> > handle_genPartID    ;      
  edm::Handle< std::vector<float> > handle_genPartMass  ;
  edm::Handle< std::vector<float> > handle_genPartMomID ;      
  edm::Handle< std::vector<float> > handle_genPartPhi   ;      
  edm::Handle< std::vector<float> > handle_genPartPt    ;      
  edm::Handle< std::vector<float> > handle_genPartStatus;      
  edm::Handle< std::vector<float> > handle_genPartY     ;      

if (readGen_) {

  iEvent.getByLabel(edm::InputTag("genPart","genPartCharge"), handle_genPartCharge);
  iEvent.getByLabel(edm::InputTag("genPart","genPartE"     ), handle_genPartE     );
  iEvent.getByLabel(edm::InputTag("genPart","genPartEta"   ), handle_genPartEta   );
  iEvent.getByLabel(edm::InputTag("genPart","genPartID"    ), handle_genPartID    );
  iEvent.getByLabel(edm::InputTag("genPart","genPartMass"  ), handle_genPartMass  );
  iEvent.getByLabel(edm::InputTag("genPart","genPartMomID" ), handle_genPartMomID );
  iEvent.getByLabel(edm::InputTag("genPart","genPartPhi"   ), handle_genPartPhi   );
  iEvent.getByLabel(edm::InputTag("genPart","genPartPt"    ), handle_genPartPt    );
  iEvent.getByLabel(edm::InputTag("genPart","genPartStatus"), handle_genPartStatus);
  iEvent.getByLabel(edm::InputTag("genPart","genPartY"     ), handle_genPartY     );

  if (!handle_genPartCharge.isValid()) return;
  if (!handle_genPartE.isValid())      return;
  if (!handle_genPartEta.isValid())    return;
  if (!handle_genPartID.isValid())     return;
  if (!handle_genPartMass.isValid())   return;
  if (!handle_genPartMomID.isValid())  return;
  if (!handle_genPartPhi.isValid())    return;
  if (!handle_genPartPt.isValid())     return;
  if (!handle_genPartStatus.isValid()) return;
  if (!handle_genPartY.isValid())      return;

  const std::vector<float> genPartCharge = *(handle_genPartCharge.product());
  const std::vector<float> genPartE      = *(handle_genPartE.product()     );
  const std::vector<float> genPartEta    = *(handle_genPartEta.product()   );
  const std::vector<float> genPartID     = *(handle_genPartID.product()    );
  const std::vector<float> genPartMass   = *(handle_genPartMass.product()  );
  const std::vector<float> genPartMomID  = *(handle_genPartMomID.product() );
  const std::vector<float> genPartPhi    = *(handle_genPartPhi.product()   );
  const std::vector<float> genPartPt     = *(handle_genPartPt.product()    );
  const std::vector<float> genPartStatus = *(handle_genPartStatus.product());
  const std::vector<float> genPartY      = *(handle_genPartY.product()     );



  float genTopPt  = 0 ;	
  float genTopPhi = 0 ;	
  float genTopEta = 0 ;	
  float genTopE   = 0 ;
	
  float genAntiTopPt  = 0;	
  float genAntiTopPhi = 0;	
  float genAntiTopEta = 0;	
  float genAntiTopE   = 0;


  int genSize = genPartCharge.size();

  //for ( int i = 0; i < genSize; i++) {

  //	if( genPartID.at(i) == 6.0 ){

  //		printf("%i \t %f \t %f  \t %f \t %f \t %f \t %f \n", i, genPartID.at(i), genPartStatus.at(i), genPartMass.at(i), genPartPt.at(i), genPartEta.at(i), genPartPhi.at(i)); 

  //	}

  //}

  //printf("-----------\n");
	

  for ( int i = 0; i < genSize; i++) {

	if ( genPartID.at(i) ==  6.0 ) { 

		genTopPt  = genPartPt.at(i) ; 
		genTopPhi = genPartPhi.at(i);
		genTopEta = genPartEta.at(i);  	
  		genTopE   = genPartE.at(i)  ;

		//printf("i-top = %i \n", i );

		break;

	}

  }


  for ( int i = 0; i < genSize; i++) {
 
        if ( genPartID.at(i) == -6.0 ) { 

		genAntiTopPt  = genPartPt.at(i) ; 
		genAntiTopPhi = genPartPhi.at(i);
		genAntiTopEta = genPartEta.at(i);  	
  		genAntiTopE   = genPartE.at(i)  ;

		//printf("i-antitop = %i \n", i );

		break;

	}
	
  }


  t_genTopPt      = genTopPt     ;	
  t_genTopPhi     = genTopPhi    ;	
  t_genTopEta     = genTopEta    ;	
  t_genTopE       = genTopE      ;	
  t_genAntiTopPt  = genAntiTopPt ;	
  t_genAntiTopPhi = genAntiTopPhi;	
  t_genAntiTopEta = genAntiTopEta;	
  t_genAntiTopE   = genAntiTopE  ;


}   // end readGen

  int triggerSize = triggerBitTree.size();

  for ( int i = 0; i < triggerSize; i++ ) {

	t_triggerBitTree      -> push_back(triggerBitTree.at(i)     );
	t_triggerPrescaleTree -> push_back(triggerPrescaleTree.at(i));
	t_triggerNameTree     -> push_back(triggerNameTree.at(i)    );

  }



// channel assignment
//---------------------------------------------------------------------------------- 
    int ch; 

    if ( lep1.Flavor * lep2.Flavor == 1 ) ch = 0;   // ee
    if ( lep1.Flavor * lep2.Flavor == 4 ) ch = 1;   // mumu
    if ( lep1.Flavor * lep2.Flavor == 2 ) ch = 2;   // emu 


//---------------------------------------------------------------------------------- 


  // Fill the tree and delete the pointers
  // Replace (1) by the skim selection
  t_metPt   = metPt.at(0);

  t_ch = ch;

  t_lep1Pt  = lep1.Pt ;
  t_lep1Eta = lep1.Eta;
  t_lep1Phi = lep1.Phi;
  t_lep1E   = lep1.E  ;

  t_lep2Pt  = lep2.Pt ;
  t_lep2Eta = lep2.Eta;
  t_lep2Phi = lep2.Phi;
  t_lep2E   = lep2.E  ;

  t_dilepInvMass = m_ll;

  t_jet1Pt  = jet1_tlv.Pt() ;
  t_jet1Eta = jet1_tlv.Eta();
  t_jet1Phi = jet1_tlv.Phi();
  t_jet1E   = jet1_tlv.E()  ;

  t_jet2Pt  = jet2_tlv.Pt() ;
  t_jet2Eta = jet2_tlv.Eta();
  t_jet2Phi = jet2_tlv.Phi();
  t_jet2E   = jet2_tlv.E()  ;

  t_C_a           = C_a          ;
  t_C_b           = C_b          ;
  t_D_a           = D_a          ;
  t_D_b           = D_b          ;
  t_aplanarity_a  = aplanarity_a ;
  t_aplanarity_b  = aplanarity_b ;
  t_circularity_a = circularity_a;
  t_circularity_b = circularity_b;
  t_isotropy_a    = isotropy_a   ;
  t_isotropy_b    = isotropy_b   ;
  t_sphericity_a  = sphericity_a ;
  t_sphericity_b  = sphericity_b ;
  t_thrust_a      = thrust_a     ;
  t_thrust_b      = thrust_b     ;
  t_centrality    = centrality   ;


  if (1)
    {
      DMTree->Fill();
    }

  delete t_triggerBitTree;
  delete t_triggerPrescaleTree;
  delete t_triggerNameTree;

//---------------------------------------------------------------------------------- 




/*float Zradius      = 15. ;

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
h_LeptonDeltaPhi_all_04->Fill(fabs( ROOT::Math::VectorUtil::DeltaPhi(lep1_tlv, lep2_tlv) )); */ 





}   // end DMAnalysisTreeMaker::analyze




void DMAnalysisTreeMaker::beginJob()
{
  edm::Service<TFileService> fs;

  DMTree = fs->make<TTree>("DMTree", "DMTree", 0);

  DMTree->Branch("t_metPt", &t_metPt, "t_metPt/F");

  DMTree->Branch("t_ch", &t_ch, "t_ch");
  
  DMTree->Branch("t_lep1Pt" , &t_lep1Pt , "t_lep1Pt/F" );
  DMTree->Branch("t_lep1Eta", &t_lep1Eta, "t_lep1Eta/F"); 
  DMTree->Branch("t_lep1Phi", &t_lep1Phi, "t_lep1Phi/F"); 
  DMTree->Branch("t_lep1E"  , &t_lep1E  , "t_lep1E/F"  ); 

  DMTree->Branch("t_lep2Pt" , &t_lep2Pt , "t_lep2Pt/F" );  
  DMTree->Branch("t_lep2Eta", &t_lep2Eta, "t_lep2Eta/F");
  DMTree->Branch("t_lep2Phi", &t_lep2Phi, "t_lep2Phi/F");
  DMTree->Branch("t_lep2E"  , &t_lep2E  , "t_lep2E/F"  );  

  DMTree->Branch("t_dilepInvMass", &t_dilepInvMass , "t_dilepInvMass/F" ); 

  DMTree->Branch("t_jet1Pt" , &t_jet1Pt , "t_jet1Pt/F" ); 
  DMTree->Branch("t_jet1Eta", &t_jet1Eta, "t_jet1Eta/F");
  DMTree->Branch("t_jet1Phi", &t_jet1Phi, "t_jet1Phi/F");
  DMTree->Branch("t_jet1E"  , &t_jet1E  , "t_jet1E/F"  ); 
  DMTree->Branch("t_jet2Pt" , &t_jet2Pt , "t_jet2Pt/F" );
  DMTree->Branch("t_jet2Eta", &t_jet2Eta, "t_jet2Eta/F");
  DMTree->Branch("t_jet2Phi", &t_jet2Phi, "t_jet2Phi/F");
  DMTree->Branch("t_jet2E"  , &t_jet2E  , "t_jet2E/F"  );  

  DMTree->Branch("t_C_a"          , &t_C_a          , "t_C_a/F"          );     
  DMTree->Branch("t_C_b"          , &t_C_b          , "t_C_b/F"          );          
  DMTree->Branch("t_D_a"          , &t_D_a          , "t_D_a/F"          );          
  DMTree->Branch("t_D_b"          , &t_D_b          , "t_D_b/F"          );          
  DMTree->Branch("t_aplanarity_a" , &t_aplanarity_a , "t_aplanarity_a/F" ); 
  DMTree->Branch("t_aplanarity_b" , &t_aplanarity_b , "t_aplanarity_b/F" ); 
  DMTree->Branch("t_circularity_a", &t_circularity_a, "t_circularity_a/F");
  DMTree->Branch("t_circularity_b", &t_circularity_b, "t_circularity_b/F");
  DMTree->Branch("t_isotropy_a"   , &t_isotropy_a   , "t_isotropy_a/F"   );   
  DMTree->Branch("t_isotropy_b"   , &t_isotropy_b   , "t_isotropy_b/F"   );  
  DMTree->Branch("t_sphericity_a" , &t_sphericity_a , "t_sphericity_a/F" );
  DMTree->Branch("t_sphericity_b" , &t_sphericity_b , "t_sphericity_b/F" );
  DMTree->Branch("t_thrust_a"     , &t_thrust_a     , "t_thrust_a/F"     );     
  DMTree->Branch("t_thrust_b"     , &t_thrust_b     , "t_thrust_b/F"     );     
  DMTree->Branch("t_centrality"   , &t_centrality   , "t_centrality/F"   );   


if (readGen_) {  

  DMTree->Branch("t_genTopPt"     , &t_genTopPt     , "t_genTopPt/F"      ); 
  DMTree->Branch("t_genTopPhi"    , &t_genTopPhi    , "t_genTopPhi/F"     ); 
  DMTree->Branch("t_genTopEta"    , &t_genTopEta    , "t_genTopEta/F"     ); 
  DMTree->Branch("t_genTopE"      , &t_genTopE      , "t_genTopE/F"       ); 
  DMTree->Branch("t_genAntiTopPt" , &t_genAntiTopPt , "t_genAntiTopPt/F"  ); 
  DMTree->Branch("t_genAntiTopPhi", &t_genAntiTopPhi, "t_genAntiTopPhi/F" ); 
  DMTree->Branch("t_genAntiTopEta", &t_genAntiTopEta, "t_genAntiTopEta/F" ); 
  DMTree->Branch("t_genAntiTopE"  , &t_genAntiTopE  , "t_genAntiTopE/F"   ); 

}

  DMTree->Branch( "t_triggerBitTree",      "std::vector<float>",       &t_triggerBitTree      );
  DMTree->Branch( "t_triggerPrescaleTree", "std::vector<int>",         &t_triggerPrescaleTree );
  DMTree->Branch( "t_triggerNameTree",     "std::vector<std::string>", &t_triggerNameTree     );

  //DMTree->Branch("t_muPt",  "std::vector<float>", &t_muPt);   // vectorial branch
 



//histo definition at Mkhistogr2.c:  h[i][l][m] = new TH1F( id[i][l][m],    ntu[i][1],    blup[i][0], blup[i][1], blup[i][2] )     with  id[i][l][m] = ntu[i][0] + "_" + namech[l] + "_" + namecut[m]

/* h_MET_all_01            = fs->make<TH1F>( "MET_all_01"           ,  "E_{T}^{miss} [GeV]"               ,  15, 0 , 600  );
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
 h_LeptonDeltaPhi_all_04 = fs->make<TH1F>( "LeptonDeltaPhi_all_04",  "#Delta#phi_{ll}"                  ,  16, 0 , 3.2  ); */

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
