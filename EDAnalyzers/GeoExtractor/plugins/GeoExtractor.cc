// -*- C++ -*-
//
// Package:    EDAnalyzers/GeoExtractor
// Class:      GeoExtractor
//
/**\class GeoExtractor GeoExtractor.cc EDAnalyzers/GeoExtractor/plugins/GeoExtractor.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:
//         Created:  Sat, 11 May 2019 13:14:55 GMT
//
//

// system include files
#include <memory>

// user include files

#include "CommonTools/UtilAlgos/interface/TFileService.h"
// #include "DataFormats/CaloTowers/interface/CaloTowerDefs.h"
// #include "DataFormats/Common/interface/MapOfVectors.h"
// #include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
// #include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/FWLite/interface/ESHandle.h"
#include "DataFormats/ForwardDetId/interface/HGCEEDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
// #include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
// #include "DataFormats/HGCalReco/interface/Trackster.h"
// #include "DataFormats/HGCRecHit/interface/HGCRecHit.h"
// #include "DataFormats/HepMCCandidate/interface/GenParticle.h"
// #include "DataFormats/JetReco/interface/PFJet.h"
// #include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/ParticleFlowReco/interface/HGCalMultiCluster.h"
// #include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
// #include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
// #include "DataFormats/TrackReco/interface/Track.h"
// #include "DataFormats/TrackReco/interface/TrackFwd.h"
// #include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "Geometry/CaloTopology/interface/HGCalTopology.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/Records/interface/HGCalGeometryRecord.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
// #include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
// #include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"
// #include "SimDataFormats/CaloHit/interface/PCaloHit.h"
// #include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
// #include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "EDAnalyzers/GeoExtractor/interface/TreeOutputInfo.h"
#include "EDAnalyzers/GeoExtractor/interface/DetO.h"

#include <Compression.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMatrixD.h>
#include <TTree.h>
#include <TVector2.h>
#include <TVectorD.h>

#include <fstream>
// #include <iostream>
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

class GeoExtractor : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
public:
  explicit GeoExtractor(const edm::ParameterSet &);
  ~GeoExtractor();

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event &, const edm::EventSetup &) override;
  virtual void endJob() override;

  edm::ESHandle<CaloGeometry> geom;
  hgcal::RecHitTools recHitTools;
  std::ofstream myfile;
  edm::Service<TFileService> fs;
  TreeOutputInfo::TreeOutput *treeOutput;
  //Setup the map, that contains all the Classures
  //Det -> SubDet -> Layer -> Wafer -> Cell
  DetColl detcol;
};

GeoExtractor::GeoExtractor(const edm::ParameterSet &iConfig)
{
  myfile.open("geometry.yaml");
  myfile.clear();
  usesResource("TFileService");
  treeOutput = new TreeOutputInfo::TreeOutput("tree", fs);
}
GeoExtractor::~GeoExtractor() { myfile.close(); }

//
// member functions
//


// ------------ method called for each event  ------------
void GeoExtractor::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
  iSetup.get<CaloGeometryRecord>().get(geom);
  recHitTools.setGeometry(*(geom.product()));

  std::map<int, edm::ESHandle<HGCalTopology> > m_topo;
  m_topo[DetId::HGCalEE];
  iSetup.get<IdealGeometryRecord>().get("HGCalEESensitive", m_topo[DetId::HGCalEE]);
  m_topo[DetId::HGCalHSi];
  iSetup.get<IdealGeometryRecord>().get("HGCalHESiliconSensitive", m_topo[DetId::HGCalHSi]);
  m_topo[DetId::HGCalHSc];
  iSetup.get<IdealGeometryRecord>().get("HGCalHEScintillatorSensitive", m_topo[DetId::HGCalHSc]);

  int n_printed = 0;
  //get all valid cells in the geometry, will be filtered later
  const std::vector<DetId> v_allCellIds = geom->getValidDetIds();

  printf("#All cell ids %i\n", (int)v_allCellIds.size());
  // Filter the Ids
  std::vector<DetId> v_detId;
  int rejected_det = 0;
  int rejected_pos = 0;
  for (int i = 0; i < (int)v_allCellIds.size(); i++)
  {
    // Skip IDs from other detector parts
    if (v_allCellIds[i].det() != DetId::HGCalEE && v_allCellIds[i].det() != DetId::HGCalHSi &&
        v_allCellIds[i].det() != DetId::HGCalHSc)
    {
      rejected_det++;
      continue;
    }
    // Todo Position check
    auto x = recHitTools.getPosition(v_allCellIds[i]).x();
    if (x < -50 || x > 50)
    {
      rejected_pos++;
      continue;
    }

    auto y = recHitTools.getPosition(v_allCellIds[i]).y();
    if (y < -150 || y > -50)
    {
      rejected_pos++;
      continue;
    }
    // If all checks are pass, add the detId to the list of valid Ids.
    v_detId.push_back(v_allCellIds[i]);
  }
  printf("Rejected by detector: %i\n", rejected_det);
  printf("Rejected by position: %i\n", rejected_pos);
  printf("Cells left: %i\n", (int)v_detId.size());


  for (int i = 0; i < (int)v_detId.size(); i++)
  {
    DetId cID = v_detId[i];
    // check angle, forward

    // Logmessage
    if (i % 100 == 0)
    {
      printf("Processing %i\n", i);
    }
    // Setup the geometry
    edm::ESHandle<HGCalTopology> &handle_topo_HGCal = m_topo[cID.det()];

    if (!handle_topo_HGCal.isValid())
    {
      printf("Error: Invalid HGCal topology. \n");
      exit(EXIT_FAILURE);
    }

    // Setup the detector
    int detectorid = cID.det();
    // if detector not in map, initialize it;
    if (detcol.detectors.find(detectorid) == detcol.detectors.end())
    {
      detcol.detectors[detectorid];
    }
    Det &detector = detcol.detectors[detectorid];

    //Setup the subdetector
    int subdetid = cID.subdetId();
    if (detector.subdetectors.find(subdetid) == detector.subdetectors.end())
    {
      detector.subdetectors[subdetid];
    }
    Subdet &subdet = detector.subdetectors[subdetid];

    //Setup the layer
    int layerid = recHitTools.getLayer(cID);
    if (subdet.layers.find(layerid) == subdet.layers.end())
    {
      subdet.layers[layerid];
      subdet.layers[layerid].z = recHitTools.getPosition(cID).z();
    }
    Layer &layer = subdet.layers[layerid];

    //Setup the Wafer
    std::pair<int, int> waferid = recHitTools.getWafer(cID);
    if (layer.wafers.find(waferid) == layer.wafers.end())
    {
      layer.wafers[waferid];
      // NA
      // layer.wafers[waferid].middle_x= recHitTools.getPosition(cID)[0];
      // layer.wafers[waferid].middle_y= recHitTools.getPosition(cID)[1];
      layer.wafers[waferid].si_thickness = recHitTools.getSiThickIndex(cID);
    }
    Wafer &wafer = layer.wafers[waferid];

    //Setup the Cell
    wafer.cells[cID];
    Cell &cell = wafer.cells[cID];

    cell.Id = cID;
    cell.x = recHitTools.getPosition(cID).x();
    cell.y = recHitTools.getPosition(cID).y();
    cell.neighbors = handle_topo_HGCal->neighbors(cID);

    treeOutput->cellid.push_back(cell.Id);
    treeOutput->detectorid.push_back(detectorid);
    treeOutput->subdetid.push_back(subdetid);
    treeOutput->layerid.push_back(layerid);
    treeOutput->waferid.push_back(waferid);

    // myfile << "\tDet " << cID.det();
    // myfile << "\tSubdet " << subdetid;
    // myfile << "\tLayer " << layerid;
    // myfile << "\tWafer (" << waferid.first << "," << waferid.second << ")";
    // myfile << "\tCellID " << cell.Id.rawId();
    // myfile << "\tx " << cell.x;
    // myfile << "\ty " << cell.y;
    // myfile << "\tneighbors ";
    // for (DetId elem : cell.neighbors)
    // {
    //   myfile << elem.rawId() << ", ";
    // }
    // myfile << "\n";
    n_printed++;
    if (n_printed > 1000)
    {
      detcol.toyaml(myfile, 0);
      treeOutput->fill();

      return;
    }
  }
}

// ------------ method called once each job just before starting event loop  ------------
void GeoExtractor::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void GeoExtractor::endJob() {}

//define this as a plug-in
DEFINE_FWK_MODULE(GeoExtractor);

// for (int i = 0; i < (int)v_detId.size(); i++) {

//     // Skip IDs from other detector parts
//     if (v_detId[i].det() < 8) {
//       continue;
//     }
//     // Logmessage
//     if (i%100 == 0 ) {
//       printf("Processing %i", i);
//     }
//     DetId curId = v_detId[i];
//     myfile << i << ": Id: " << curId.det();
//     myfile << "\tSubdetector " << curId.subdetId();
//     myfile << "\tPosition" << recHitTools.getPosition(curId);
//     myfile << "\tWafer Id" << recHitTools.getWafer(curId).first << "," << recHitTools.getWafer(curId).second << "\n";
//     n_printed++;
//     if (n_printed > 1000) {
//       return;
//     }
//   }