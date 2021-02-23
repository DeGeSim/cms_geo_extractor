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
#include <stdexcept>
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
  //default edanalyzer fundionst
  virtual void beginJob() override;
  virtual void analyze(const edm::Event &, const edm::EventSetup &) override;
  virtual void endJob() override;

  //Tools
  edm::ESHandle<CaloGeometry> geom;
  hgcal::RecHitTools recHitTools;

  // container for the topologyies
  std::map<int, edm::ESHandle<HGCalTopology> > m_topo;

  //// IO
  // file to write the yaml structured detector information out
  std::ofstream myfile;
  edm::Service<TFileService> fs;
  // output tree
  TreeOutputInfo::TreeOutput *treeOutput;

  std::vector<DetId> filterCellIds(const std::vector<DetId> v_allCellIds);
  //Stuff or seaching the members

  void assignZneighbors(std::vector<DetId> &v_validHGCalIds);
  std::pair<DetId, float> findCellCloseToXYpos(DetId cellId, unsigned int detectorid, unsigned int subdetid, unsigned int layerid);
  DetId getcell(unsigned int detectorid, unsigned int subdetid, unsigned int layerid, std::pair<int, int> wafer, std::pair<int, int> cell);
  DetId findNextCell(DetId cellId);
  std::string printcell(unsigned int detectorid, unsigned int subdetid, unsigned int layerid, std::pair<int, int> wafer, std::pair<int, int> cell);

  //The map, that contains the detector structure
  //Det -> SubDet -> Layer -> Wafer -> Cell
  DetColl detcol;

  //vector with the numbers of the detector part of the hgcal
  std::vector<int> v_HGCalDets;

  //map that saves which cell are rejected in which step
  std::map<int, std::map<std::string, int> > m_rej;
};

GeoExtractor::GeoExtractor(const edm::ParameterSet &iConfig)
{
  usesResource("TFileService");
  treeOutput = new TreeOutputInfo::TreeOutput("tree", fs);
  m_topo[DetId::HGCalEE];
  m_topo[DetId::HGCalHSi];
  m_topo[DetId::HGCalHSc];

  v_HGCalDets.push_back(DetId::HGCalEE);
  v_HGCalDets.push_back(DetId::HGCalHSi);
  v_HGCalDets.push_back(DetId::HGCalHSc);
}
GeoExtractor::~GeoExtractor()
{
  treeOutput->fill();
  myfile.open("output/geometry.yaml");
  myfile.clear();
  detcol.toyaml(myfile, 0);
  myfile.close();
}

//
// member functions
//

std::vector<DetId> GeoExtractor::filterCellIds(const std::vector<DetId> v_allCellIds)
{
  std::vector<DetId> v_validHGCalIds;
  printf("#All cell ids %i\n", (int)v_allCellIds.size());
  for (auto detectorid : v_HGCalDets)
  {
    m_rej[detectorid];
    m_rej[detectorid]["det"] = 0;
    m_rej[detectorid]["x"] = 0;
    m_rej[detectorid]["y"] = 0;
    m_rej[detectorid]["z"] = 0;
  }
  for (int i = 0; i < (int)v_allCellIds.size(); i++)
  {
    auto detectorid = v_allCellIds[i].det();
    // Skip IDs from other detector parts
    if (detectorid != DetId::HGCalEE && detectorid != DetId::HGCalHSi &&
        detectorid != DetId::HGCalHSc)
    {
      continue;
    }
    m_rej[detectorid]["det"]++;
    // Todo Position check

    auto x = recHitTools.getPosition(v_allCellIds[i]).x();
    if (x < -50 || x > 50)
    {
      continue;
    }
    m_rej[detectorid]["x"]++;

    auto y = recHitTools.getPosition(v_allCellIds[i]).y();
    if (y < -150 || y > -50)
    {
      continue;
    }
    m_rej[detectorid]["y"]++;

    auto z = recHitTools.getPosition(v_allCellIds[i]).z();
    if (z < 0)
    {
      continue;
    }
    m_rej[detectorid]["z"]++;
    // If all checks are pass, add the detId to the list of valid Ids.
    v_validHGCalIds.push_back(v_allCellIds[i]);
  }
  std::cout << "\tdet \t";
  std::cout << "x \t";
  std::cout << "y \t";
  std::cout << "z \n";
  for (auto detectorid : v_HGCalDets)
  {
    std::cout << detectorid << "\t";
    std::cout << m_rej[detectorid]["det"] << "\t";
    std::cout << m_rej[detectorid]["x"] << "\t";
    std::cout << m_rej[detectorid]["y"] << "\t";
    std::cout << m_rej[detectorid]["z"] << "\n";
  }
  printf("Cells left: %i\n", (int)v_validHGCalIds.size());
  return v_validHGCalIds;
}

// ------------ method called for each event  ------------
void GeoExtractor::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
  iSetup.get<CaloGeometryRecord>().get(geom);
  recHitTools.setGeometry(*(geom.product()));

  iSetup.get<IdealGeometryRecord>().get("HGCalEESensitive", m_topo[DetId::HGCalEE]);
  iSetup.get<IdealGeometryRecord>().get("HGCalHESiliconSensitive", m_topo[DetId::HGCalHSi]);
  iSetup.get<IdealGeometryRecord>().get("HGCalHEScintillatorSensitive", m_topo[DetId::HGCalHSc]);

  int n_printed = 0;
  //get all valid cells in the geometry, will be filtered later
  const std::vector<DetId> v_allCellIds = geom->getValidDetIds();

  // Filter the Ids
  //vector to save the ids
  std::vector<DetId> v_validHGCalIds;
  v_validHGCalIds = filterCellIds(v_allCellIds);

  for (int i = 0; i < (int)v_validHGCalIds.size(); i++)
  {
    DetId cID = v_validHGCalIds[i];
    // Logmessage
    // if (i % 1000 == 0)
    // {
    //   printf("Processing %i\n", i);
    // }
    // Setup the geometry
    edm::ESHandle<HGCalTopology> &handle_topo_HGCal = m_topo[cID.det()];

    if (!handle_topo_HGCal.isValid())
    {
      printf("Error: Invalid HGCal topology. \n");
      exit(EXIT_FAILURE);
    }

    // Setup the detector
    unsigned int detectorid = cID.det();
    // if detector not in map, initialize it;
    if (detcol.detectors.find(detectorid) == detcol.detectors.end())
    {
      detcol.detectors[detectorid];
    }
    Det &detector = detcol.detectors[detectorid];

    //Setup the subdetector
    unsigned int subdetid = cID.subdetId();
    if (detector.subdetectors.find(subdetid) == detector.subdetectors.end())
    {
      detector.subdetectors[subdetid];
    }
    Subdet &subdet = detector.subdetectors[subdetid];

    //Setup the layer
    unsigned int layerid = recHitTools.getLayer(cID);
    if (subdet.layers.find(layerid) == subdet.layers.end())
    {
      subdet.layers[layerid];
      subdet.layers[layerid].z = recHitTools.getPosition(cID).z();
    }
    Layer &layer = subdet.layers[layerid];

    std::pair<int, int> cellid;
    std::pair<int, int> waferid;
    Cell *cellptr;
    // setup the tiles or wafer depending on the detector
    if (detectorid == DetId::HGCalEE || detectorid == DetId::HGCalHSi)
    {
      //Setup the Wafer
      waferid = recHitTools.getWafer(cID);
      if (layer.wafers.find(waferid) == layer.wafers.end())
      {
        layer.wafers[waferid];
        layer.wafers[waferid].si_thickness = recHitTools.getSiThickIndex(cID);
      }
      Wafer &wafer = layer.wafers[waferid];

      cellid = recHitTools.getCell(cID);
      wafer.cells[cellid];
      cellptr = &wafer.cells[cellid];
      // std::cout << "Adding cell:" << printcell(detectorid, subdetid, layerid, waferid, cellid) << "\n";
    }
    else if (detectorid == DetId::HGCalHSc)
    {
      HGCScintillatorDetId scid = HGCScintillatorDetId(cID);
      cellid = scid.ietaphi();
      if (layer.tiles.find(cellid) == layer.tiles.end())
      {
        layer.tiles[cellid];
      }
      cellptr = &layer.tiles[cellid];
      // std::cout << "Adding cell:" << printcell(detectorid, subdetid, layerid, cellid, std::pair<int, int>(0,0)) << "\n";
    }
    Cell &cell = *cellptr;
    if (detectorid != DetId::HGCalHSc)
    {
      continue;
    }

    cell.globalid = cID;
    cell.x = recHitTools.getPosition(cID).x();
    cell.y = recHitTools.getPosition(cID).y();

    if (detectorid == DetId::HGCalHSc)
    {
      HGCScintillatorDetId scid = HGCScintillatorDetId(cID);
      cell.type = scid.type();
    }

    cell.neighbors = handle_topo_HGCal->neighbors(cID);

    treeOutput->globalid.push_back(cell.globalid.rawId());
    treeOutput->detectorid.push_back(detectorid);
    treeOutput->subdetid.push_back(subdetid);
    treeOutput->layerid.push_back(layerid);
    if (detectorid == DetId::HGCalHSc)
    {
      treeOutput->waferortileid.push_back(cellid);
      treeOutput->cellid.push_back(std::pair<int, int>(0,0));
    }
    else
    {
      treeOutput->waferortileid.push_back(waferid);
      treeOutput->cellid.push_back(cellid);
    }

    n_printed++;
  }
  assignZneighbors(v_validHGCalIds);
}

// wraps findNextCell and loops over the ids
void GeoExtractor::assignZneighbors(std::vector<DetId> &v_validHGCalIds)
{
  for (int i = 0; i < (int)v_validHGCalIds.size(); i++)
  {
    // if (i <253000){
    //   continue;
    // }
    DetId cID = v_validHGCalIds[i];
    if (i % 1000 == 0)
    {
      printf("Assinging z neighbors %i\n", i);
    }

    //Setup the detector
    unsigned int detectorid = cID.det();
    Det &detector = detcol.detectors[detectorid];

    //Setup the subdetector
    unsigned int subdetid = cID.subdetId();
    Subdet &subdet = detector.subdetectors[subdetid];

    //Setup the layer
    unsigned int layerid = recHitTools.getLayer(cID);

    Cell *cellptr;
    Layer &layer = subdet.layers[layerid];

    if (detectorid == DetId::HGCalEE || detectorid == DetId::HGCalHSi)
    {
      // get the cellid
      std::pair<int, int> cellid = recHitTools.getCell(cID);
      //Setup the Wafer
      std::pair<int, int> waferid = recHitTools.getWafer(cID);
      cellptr = &layer.wafers[waferid].cells[cellid];
    }
    else if (detectorid == DetId::HGCalHSi)
    {
      //Setup the Wafer
      HGCScintillatorDetId scid = HGCScintillatorDetId(cID);
      std::pair<int, int> cellid = scid.ietaphi();
      cellptr = &layer.tiles[cellid];
    }

    Cell &cell = *cellptr;
    cell.next = findNextCell(cID);
  }
}

// This function decices in which detector(s) to search for the neigbor.
// The heavy lifting is done by the findCellCloseToXYpos
DetId GeoExtractor::findNextCell(DetId cellId)
{
  unsigned int detectorid = cellId.det();
  unsigned int subdetid = cellId.subdetId();
  unsigned int layerid = recHitTools.getLayer(cellId);

  //
  // HGCalEE = 8, layer 1-28
  // HGCalHSi = 9, layer 1-22
  // HGCalHSc = 10, layer 9-22
  // HGCalTrigger = 11 X

  // For the ee cal we can easily move forward
  if (detectorid == DetId::HGCalEE)
  {
    if (layerid < 22)
    {
      printf("A");
      return findCellCloseToXYpos(cellId, detectorid, subdetid, layerid + 1).first;
    }
    // to the the next detector
    else
    {
      printf("B");
      return findCellCloseToXYpos(cellId, DetId::HGCalHSi, subdetid, 1).first;
    }
  }
  if (detectorid == DetId::HGCalHSi || detectorid == DetId::HGCalHSc)
  {
    // for layer <8 all cells we can just search in the Si part
    if (layerid < 8 && detectorid == DetId::HGCalHSi)
    {
      printf("C");
      return findCellCloseToXYpos(cellId, DetId::HGCalHSi, subdetid, layerid + 1).first;
    }
    else if (layerid == 22)
    {
      return DetId(0);
    }
    else
    {
      // printf("D");
      auto [sicanid, deltasi] = findCellCloseToXYpos(cellId, DetId::HGCalHSi, subdetid, layerid + 1);
      auto [sccanid, deltasc] = findCellCloseToXYpos(cellId, DetId::HGCalHSc, subdetid, layerid + 1);
      if (deltasi < deltasc)
      {
        return sicanid;
      }
      else
      {
        return sccanid;
      }
    }
  }
  printf("Det %i, Subdet %i, Layer %i", detectorid, subdetid, layerid);
  throw std::invalid_argument("Wont find neighbor. This part should never be reached.\n");
  return DetId(0);
}

// This is the function that does the acutal search.
// the coordinates give the layer, that is to be searched
std::pair<DetId, float> GeoExtractor::findCellCloseToXYpos(
    DetId cellId,
    unsigned int detectorid,
    unsigned int subdetid,
    unsigned int layerid)
{
  // Get the topo of the detector
  edm::ESHandle<HGCalTopology> &handle_topo_HGCal = m_topo[detectorid];
  //Get xy from the origin cell
  float x = recHitTools.getPosition(cellId).x();
  float y = recHitTools.getPosition(cellId).y();
  std::pair<int, int> wafer = recHitTools.getWafer(cellId);
  std::pair<int, int> cell = recHitTools.getCell(cellId);

  // The cell that is the closest, to be replace by close cells
  // if getcell can't find the wafer / cell it starts with the (0,0) coordinates
  DetId closest_cell = getcell(detectorid, subdetid, layerid, wafer, cell);

  printf("3\n");
  printf("closest_cell %u, det %u, subdet %u\n", closest_cell.rawId(),closest_cell.det(), closest_cell.subdetId());
  printf("layer %u\n", recHitTools.getLayer(closest_cell));

  float x_cur = recHitTools.getPosition(closest_cell).x();
  float y_cur = recHitTools.getPosition(closest_cell).y();
  float d_cur = (x_cur - x) * (x_cur - x) + (y_cur - y) * (y_cur - y);
  bool improvement = true;
  printf("3.5 ");
  while (improvement)
  {
    improvement = false;

    for (DetId &neighbor : handle_topo_HGCal->neighbors(closest_cell))
    {
      printf("4");
      float xn = recHitTools.getPosition(closest_cell).x();
      float yn = recHitTools.getPosition(closest_cell).y();
      float dn = (xn - x) * (xn - x) + (yn - y) * (yn - y);
      if (dn < d_cur)
      {
        improvement = true;
        closest_cell = neighbor;
        d_cur = dn;
        break;
      }
    }
  }
  return std::make_pair(closest_cell, d_cur);
}

DetId GeoExtractor::getcell(unsigned int detectorid, unsigned int subdetid, unsigned int layerid, std::pair<int, int> waferortileid, std::pair<int, int> cellid)
{
  if (detcol.detectors.find(detectorid) == detcol.detectors.end())
  {
    throw std::invalid_argument("No such detector\n" + printcell(detectorid, subdetid, layerid, waferortileid, cellid));
  }
  Det &detector = detcol.detectors[detectorid];

  if (detector.subdetectors.find(subdetid) == detector.subdetectors.end())
  {
    throw std::invalid_argument("No such subdetector\n" + printcell(detectorid, subdetid, layerid, waferortileid, cellid));
  }
  Subdet &subdet = detector.subdetectors[subdetid];

  if (subdet.layers.find(layerid) == subdet.layers.end())
  {
    throw std::invalid_argument("No such layer\n" + printcell(detectorid, subdetid, layerid, waferortileid, cellid));
  }
  Layer &layer = subdet.layers[layerid];

  if (detectorid == DetId::HGCalEE || detectorid == DetId::HGCalHSi)
  {
    if (layer.wafers.find(waferortileid) == layer.wafers.end())
    {
      waferortileid = layer.wafers.begin()->first;
    }
    Wafer &wafer = layer.wafers[waferortileid];

    if (wafer.cells.find(cellid) == wafer.cells.end())
    {
      cellid = wafer.cells.begin()->first;
    }
    Cell &cell = wafer.cells[cellid];
    return cell.globalid;
  }
  else if (detectorid == DetId::HGCalHSi)
  {
    if (layer.tiles.find(waferortileid) == layer.tiles.end())
    {
      waferortileid = layer.tiles.begin()->first;
    }
    return layer.tiles[waferortileid].globalid;
  }
  else
  {
    throw std::invalid_argument("");
  }
  return DetId(0);
}

std::string GeoExtractor::printcell(unsigned int detectorid, unsigned int subdetid, unsigned int layerid, std::pair<int, int> waferid, std::pair<int, int> cellid)
{
  std::ostringstream stringStream;
  stringStream << " Det " << detectorid;
  stringStream << " Subdet " << subdetid;
  stringStream << " Layer " << layerid;
  stringStream << " Wafer (" << waferid.first << "," << waferid.second << ")";
  stringStream << " Cell (" << cellid.first << "," << cellid.second << ")";
  std::string copyOfStr = stringStream.str();

  return copyOfStr;
}

// ------------ method called once each job just before starting event loop  ------------
void GeoExtractor::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void GeoExtractor::endJob() {}

//define this as a plug-in
DEFINE_FWK_MODULE(GeoExtractor);