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
  TH1F* h_muPt;
};


//
// constants, enums and typedefs
//

//
// static data member definitions
//


DMAnalysisTreeMaker::DMAnalysisTreeMaker(const edm::ParameterSet& iConfig) {}

DMAnalysisTreeMaker::~DMAnalysisTreeMaker() {}

void DMAnalysisTreeMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle< std::vector<float> > handle_muPt;
  iEvent.getByLabel(edm::InputTag("muons","muPt"), handle_muPt);

  if (!handle_muPt.isValid()) return;

  const std::vector<float> muPt = *(handle_muPt.product());

  int muPt_size = muPt.size();

  for (int i=0; i<muPt_size; i++) h_muPt->Fill(muPt.at(i));
}


void DMAnalysisTreeMaker::beginJob()
{
  edm::Service<TFileService> fs;

  h_muPt = fs->make<TH1F>("h_muPt", "", 200, 0, 200);
}

void DMAnalysisTreeMaker::endJob() {}


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
