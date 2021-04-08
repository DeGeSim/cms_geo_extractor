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
#include "DataFormats/ForwardDetId/interface/HGCScintillatorDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCSiliconDetId.h"

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
#include "EDAnalyzers/GeoExtractor/interface/addOutput.h"
#include "EDAnalyzers/GeoExtractor/interface/DetO.h"
#include "EDAnalyzers/GeoExtractor/interface/MakePrintable.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <Compression.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMatrixD.h>
#include <TTree.h>
#include <TVector2.h>
#include <TVectorD.h>

#include <fstream>
#include <stdexcept>
#include <stdlib.h>

// setup  the logger outside
#include "EDAnalyzers/GeoExtractor/interface/Log.h"
structlog LOGCFG = {};

// class declaration
#include "GeoExtractor.h"
#include "utils.h"
#include "filters.h"
#include "search.h"
#include "gapfixer.h"
//
// member functions
//

void GeoExtractor::instanciateMapForCell(DetId &iterId)
{
  // Setup the detector
  unsigned int detectorid = iterId.det();
  // if detector not in map, initialize it;
  if (detcol.detectors.find(detectorid) == detcol.detectors.end())
  {
    detcol.detectors[detectorid];
  }
  Det &detector = detcol.detectors[detectorid];

  //Setup the subdetector
  unsigned int subdetid = iterId.subdetId();
  if (detector.subdetectors.find(subdetid) == detector.subdetectors.end())
  {
    detector.subdetectors[subdetid];
  }
  Subdet &subdet = detector.subdetectors[subdetid];

  //Setup the layer
  unsigned int layerid = recHitTools.getLayer(iterId);
  if (subdet.layers.find(layerid) == subdet.layers.end())
  {
    subdet.layers[layerid];
    subdet.layers[layerid].z = recHitTools.getPosition(iterId).z();
  }
  Layer &layer = subdet.layers[layerid];

  std::pair<int, int> cellid;
  std::pair<int, int> waferid;
  // setup the tiles or wafer depending on the detector
  if (isSiliconDet(detectorid))
  {
    //Setup the Wafer
    waferid = recHitTools.getWafer(iterId);
    if (layer.wafers.find(waferid) == layer.wafers.end())
    {
      layer.wafers[waferid];
    }
    Wafer &wafer = layer.wafers[waferid];

    cellid = recHitTools.getCell(iterId);
    wafer.cells[cellid];
    // LOG(DEBUG) << "Adding cell:" << wafer.cells[cellid]<<"\n";
  }
  else
  {
    HGCScintillatorDetId scid = HGCScintillatorDetId(iterId);
    cellid = scid.ietaphi();
    if (layer.tiles.find(cellid) == layer.tiles.end())
    {
      layer.tiles[cellid];
    }
    // LOG(DEBUG) << "Adding cell:" << layer.tiles[cellid]<<"\n";
  }
}

// ------------ method called for each event  ------------
void GeoExtractor::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
  iSetup.get<CaloGeometryRecord>().get(geom);
  recHitTools.setGeometry(*(geom.product()));

  iSetup.get<IdealGeometryRecord>().get("HGCalEESensitive", m_geom[DetId::HGCalEE]);
  iSetup.get<IdealGeometryRecord>().get("HGCalHESiliconSensitive", m_geom[DetId::HGCalHSi]);
  iSetup.get<IdealGeometryRecord>().get("HGCalHEScintillatorSensitive", m_geom[DetId::HGCalHSc]);

  iSetup.get<IdealGeometryRecord>().get("HGCalEESensitive", m_topo[DetId::HGCalEE]);
  iSetup.get<IdealGeometryRecord>().get("HGCalHESiliconSensitive", m_topo[DetId::HGCalHSi]);
  iSetup.get<IdealGeometryRecord>().get("HGCalHEScintillatorSensitive", m_topo[DetId::HGCalHSc]);

  if (!m_geom[DetId::HGCalEE].isValid() or !m_geom[DetId::HGCalHSi].isValid() or !m_geom[DetId::HGCalHSc].isValid())
  {
    LOG(ERROR) << "Error: Invalid HGCal geometry."
               << "\n";
    exit(EXIT_FAILURE);
  }
  int n_printed = 0;
  //get all valid cells in the geometry, will be filtered later
  const std::vector<DetId> v_allCellIds = geom->getValidDetIds();

  // Filter the Ids
  LOG(INFO) << "Filter all Cells for HGCal cells."
            << "\n";
  v_validHGCalIds = filterCellIds(v_allCellIds);

  LOG(INFO) << "Filling the detector structure."
            << "\n";
  for (DetId &iterId : v_validHGCalIds)
  {
    edm::ESHandle<HGCalTopology> &handle_topo_HGCal = m_topo[iterId.det()];

    if (!handle_topo_HGCal.isValid())
    {
      LOG(ERROR) << "Error: Invalid HGCal topology."
                 << "\n";
      exit(EXIT_FAILURE);
    }
    instanciateMapForCell(iterId);
    Cell *cellptr = getCellPtr(iterId);
    auto [detectorid, subdetid, layerid, waferortileid, cellid] = getCellHash(iterId);
    cellptr->globalid = iterId;
    cellptr->x = recHitTools.getPosition(iterId).x();
    cellptr->y = recHitTools.getPosition(iterId).y();

    if (isSiliconDet(detectorid))
    {
      cellptr->issilicon = true;
      cellptr->type = recHitTools.getSiThickIndex(iterId);
    }
    else
    {
      HGCScintillatorDetId scid = HGCScintillatorDetId(iterId);

      // cellptr->type = scid.type();
      cellptr->type = scid.sipm();
      cellptr->issilicon = false;
    }
    // only assign neighbors with valid ids
    for (auto &neighbor : handle_topo_HGCal->neighbors(iterId))
    {
      if (std::find(v_validHGCalIds.begin(), v_validHGCalIds.end(), neighbor) != v_validHGCalIds.end())
      {
        cellptr->neighbors.insert(neighbor);
      }
    }

    treeOutput->globalid.push_back(cellptr->globalid.rawId());
    treeOutput->detectorid.push_back(detectorid);
    treeOutput->subdetid.push_back(subdetid);
    treeOutput->layerid.push_back(layerid);
    treeOutput->waferortileid.push_back(waferortileid);
    treeOutput->cellid.push_back(cellid);

    n_printed++;
  }
  LOG(INFO) << "Assing the Z neighbors."
            << "\n";
  assignZNeighbors(v_validHGCalIds);
  LOG(INFO) << "Done.\n";

  LOG(INFO) << "Start fixing the bounderies.\n";
  fixGap(v_validHGCalIds);
  LOG(INFO) << "Done.\n";

  LOG(INFO) << "Filling TBranches.\n";

  for (DetId &iterId : v_validHGCalIds)
  {
    Cell *cellptr = getCellPtr(iterId);

    treeOutput->x.push_back(cellptr->x);
    treeOutput->y.push_back(cellptr->y);
    treeOutput->celltype.push_back(cellptr->type);
    treeOutput->issilicon.push_back(cellptr->issilicon);
    treeOutput->next.push_back(cellptr->next);
    treeOutput->previous.push_back(cellptr->previous);

    int neighborsSize = (int)cellptr->neighbors.size();
    treeOutput->nneighbors.push_back(neighborsSize);
    int gapneighborsSize = (int)cellptr->gapneighbors.size();
    treeOutput->ngapneighbors.push_back(gapneighborsSize);

    // Put the pointer to the tbranches in a vector so we can iterate over them more easily.

    std::vector<std::vector<unsigned int> *> v_neighborTreePtrs;
    v_neighborTreePtrs.push_back(&treeOutput->n0);
    v_neighborTreePtrs.push_back(&treeOutput->n1);
    v_neighborTreePtrs.push_back(&treeOutput->n2);
    v_neighborTreePtrs.push_back(&treeOutput->n3);
    v_neighborTreePtrs.push_back(&treeOutput->n4);
    v_neighborTreePtrs.push_back(&treeOutput->n5);
    v_neighborTreePtrs.push_back(&treeOutput->n6);
    v_neighborTreePtrs.push_back(&treeOutput->n7);

    LOG(DEBUG) << iterId.rawId();
    LOG(DEBUG) << " neighborsSize:" << neighborsSize;
    LOG(DEBUG) << " gapneighborsSize:" << gapneighborsSize << "\n";

    for (int ineighbor = 0; ineighbor < 8; ineighbor++)
    {

      LOG(DEBUG) << "\t" << ineighbor << "\n";
      if (ineighbor < neighborsSize)
      {
        DetId neighbortoadd = *std::next(cellptr->neighbors.begin(), ineighbor);
        v_neighborTreePtrs[ineighbor]->push_back(neighbortoadd);

        LOG(DEBUG) << "\t" << iterId.rawId() << " ineighbor (" << ineighbor << "): adding neighbor " << neighbortoadd.rawId() << "\n";
      }

      else if (ineighbor < gapneighborsSize + neighborsSize)
      {
        DetId neighbortoadd = *std::next(cellptr->gapneighbors.begin(), ineighbor - neighborsSize);
        v_neighborTreePtrs[ineighbor]->push_back(neighbortoadd);

        LOG(DEBUG) << "\t" << iterId.rawId() << " ineighbor (" << ineighbor << "): adding gapneighbor " << neighbortoadd.rawId() << "\n";
      }
      else
      {
        LOG(DEBUG) << "\t" << iterId.rawId() << " ineighbor (" << ineighbor << "): adding 0"
                   << "\n";
        v_neighborTreePtrs[ineighbor]->push_back(0);
      }
    }
    LOG(DEBUG) << "\n";
  }
}

// ------------ method called once each job just before starting event loop  ------------
void GeoExtractor::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void GeoExtractor::endJob() {}

//define this as a plug-in
DEFINE_FWK_MODULE(GeoExtractor);