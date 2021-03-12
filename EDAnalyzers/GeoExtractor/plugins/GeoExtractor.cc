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
  DetId getstartcell(unsigned int detectorid, unsigned int subdetid, unsigned int layerid, std::pair<int, int> wafer, std::pair<int, int> cell);
  DetId findNextCell(DetId cellId);
  std::string printcell(unsigned int detectorid, unsigned int subdetid, unsigned int layerid, std::pair<int, int> wafer, std::pair<int, int> cell);
  CellHash getCellHashKeys(DetId &iterId);
  bool isSiliconDet(int ndet);

  void instanciateMapForCell(DetId &iterId);
  Cell *getCellptr(DetId &iterId);
  bool validId(DetId id);
  void validateId(DetId id);

  //The map, that contains the detector structure
  //Det -> SubDet -> Layer -> Wafer -> Cell
  DetColl detcol;

  //vector with the numbers of the detector part of the hgcal
  std::vector<int> v_HGCalDets;

  //map that saves which cell are rejected in which step
  std::map<int, std::map<std::string, int> > m_rej;
  //vector to save the ids
  std::vector<DetId> v_validHGCalIds;
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

  //set the loglevel here DEBUG < INFO < WARN < ERROR
  LOGCFG.headers = false;
  LOGCFG.level = INFO;
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
  LOG(DEBUG) << "#All cell ids " << v_allCellIds.size() << "\n";
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
  LOG(INFO) << "\tdet    ";
  LOG(INFO) << "x    ";
  LOG(INFO) << "y    ";
  LOG(INFO) << "z \n";
  for (auto detectorid : v_HGCalDets)
  {
    LOG(INFO) << detectorid << "\t";
    LOG(INFO) << m_rej[detectorid]["det"] << "\t";
    LOG(INFO) << m_rej[detectorid]["x"] << "\t";
    LOG(INFO) << m_rej[detectorid]["y"] << "\t";
    LOG(INFO) << m_rej[detectorid]["z"] << "\n";
  }
  LOG(INFO) << "Cells left: " << v_validHGCalIds.size() << "\n";
  return v_validHGCalIds;
}

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
      layer.wafers[waferid].si_thickness = recHitTools.getSiThickIndex(iterId);
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

Cell *GeoExtractor::getCellptr(DetId &iterId)
{
  auto [detectorid, subdetid, layerid, waferortileid, cellid] = getCellHashKeys(iterId);
  // Setup the detector
  Det &detector = detcol.detectors[detectorid];

  //Setup the subdetector
  Subdet &subdet = detector.subdetectors[subdetid];

  //Setup the layer
  Layer &layer = subdet.layers[layerid];

  Cell *cellptr;
  // setup the tiles or wafer depending on the detector
  if (isSiliconDet(detectorid))
  {
    //Setup the Wafer
    Wafer &wafer = layer.wafers[waferortileid];
    cellptr = &wafer.cells[cellid];
  }
  else
  {
    cellptr = &layer.tiles[waferortileid];
  }
  return cellptr;
}

CellHash GeoExtractor::getCellHashKeys(DetId &iterId)
{
  // Setup the detector
  unsigned int detectorid = iterId.det();
  Det &detector = detcol.detectors[detectorid];

  //Setup the subdetector
  unsigned int subdetid = iterId.subdetId();
  Subdet &subdet = detector.subdetectors[subdetid];

  //Setup the layer
  unsigned int layerid = recHitTools.getLayer(iterId);
  Layer &layer = subdet.layers[layerid];

  std::pair<int, int> waferortileid;
  std::pair<int, int> cellid;
  // setup the tiles or wafer depending on the detector
  if (isSiliconDet(detectorid))
  {
    //Setup the Wafer
    waferortileid = recHitTools.getWafer(iterId);
    cellid = recHitTools.getCell(iterId);
  }
  else
  {
    HGCScintillatorDetId scid = HGCScintillatorDetId(iterId);
    waferortileid = scid.ietaphi();
    cellid = std::make_pair(0, 0);
  }
  CellHash retuple = std::make_tuple(detectorid, subdetid, layerid, waferortileid, cellid);
  return retuple;
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
  LOG(INFO) << "Fillter Cells"<<"\n";
  v_validHGCalIds = filterCellIds(v_allCellIds);

  LOG(INFO) << "Build the map"<<"\n";
  for (int i = 0; i < (int)v_validHGCalIds.size(); i++)
  {
    DetId iterId = v_validHGCalIds[i];
    edm::ESHandle<HGCalTopology> &handle_topo_HGCal = m_topo[iterId.det()];

    if (!handle_topo_HGCal.isValid())
    {
      LOG(ERROR) << "Error: Invalid HGCal topology."<<"\n";
      exit(EXIT_FAILURE);
    }
    instanciateMapForCell(iterId);
    Cell *cellptr = getCellptr(iterId);
    auto [detectorid, subdetid, layerid, waferortileid, cellid] = getCellHashKeys(iterId);
    cellptr->globalid = iterId;
    cellptr->x = recHitTools.getPosition(iterId).x();
    cellptr->y = recHitTools.getPosition(iterId).y();

    if (isSiliconDet(detectorid))
    {
      cellptr->issilicon = true;
    }
    else
    {
      HGCScintillatorDetId scid = HGCScintillatorDetId(iterId);
      cellptr->type = scid.type();
      cellptr->issilicon = false;
    }
    // only assign neighbors with valid ids
    for (auto &neighbor : handle_topo_HGCal->neighbors(iterId))
    {
      if (std::find(v_validHGCalIds.begin(), v_validHGCalIds.end(), neighbor) != v_validHGCalIds.end())
      {
        cellptr->neighbors.push_back(neighbor);
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
  LOG(INFO) << "Assing the Z neighbors."<<"\n";
  assignZneighbors(v_validHGCalIds);
}

// wraps findNextCell and loops over the ids
void GeoExtractor::assignZneighbors(std::vector<DetId> &v_validHGCalIds)
{
  for (int i = 0; i < (int)v_validHGCalIds.size(); i++)
  {
    DetId iterId = v_validHGCalIds[i];
    // if (iterId.rawId() != 2426030080)
    // {
    //   continue;
    // }
    // LOG(DEBUG) << "skipping everything but 2426030080"<<"\n";
    if (i % 1000 == 0)
    {
      LOG(DEBUG) << "Assinging z neighbors " << i <<"\n";
    }
    Cell *cellptr = getCellptr(iterId);
    validateId(cellptr->globalid);
    cellptr->next = findNextCell(cellptr->globalid);
  }
}

// This function decices in which detector(s) to search for the neigbor.
// The heavy lifting is done by the findCellCloseToXYpos
DetId GeoExtractor::findNextCell(DetId cellId)
{
  LOG(DEBUG) << "Start findNextCell"<<"\n";
  CellHash hash = getCellHashKeys(cellId);
  auto [detectorid, subdetid, layerid, waferortileid, cellid] = hash;
  LOG(DEBUG) << "find cell for id" << cellId.rawId()<<"\n";
  LOG(DEBUG) << getCellHashKeys(cellId)<<"\n";
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
      LOG(DEBUG) << "A"<<"\n";
      return findCellCloseToXYpos(cellId, detectorid, subdetid, layerid + 1).first;
    }
    // to the the next detector EE-> HSi
    else
    {
      LOG(DEBUG) << "B"<<"\n";
      // TODO check if there are nodes in layer 1 of the HGCalSc
      return findCellCloseToXYpos(cellId, DetId::HGCalHSi, subdetid, 1).first;
    }
  }
  if (detectorid == DetId::HGCalHSi || detectorid == DetId::HGCalHSc)
  {
    // for layer <8 all cells we can just search in the Si part
    if (layerid < 8 && detectorid == DetId::HGCalHSi)
    {
      LOG(DEBUG) << "C"<<"\n";
      return findCellCloseToXYpos(cellId, DetId::HGCalHSi, subdetid, layerid + 1).first;
    }
    // return 0 in the last layer
    else if (layerid == 22)
    {
      LOG(DEBUG) << "D"<<"\n";
      return DetId(0);
    }
    else
    {
      LOG(DEBUG) << "E"<<"\n";
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
  // LOG(ERROR)  << hash << ""<<"\n"; //doenst work ?
  std::cout << hash;
  LOG(ERROR) << "Wont find neighbor. This part should never be reached."<<"\n";
  exit(EXIT_FAILURE);
  return DetId(0);
}

// This is the function that does the acutal search.
// the coordinates give the layer, that is to be searched
std::pair<DetId, float> GeoExtractor::findCellCloseToXYpos(
    DetId originCellDetID,
    unsigned int targetdetectorid,
    unsigned int targetsubdetid,
    unsigned int targetlayerid)
{
  LOG(DEBUG) << "Start findCellCloseToXYpos for" << originCellDetID.rawId()<<"\n";

  //Get xy from the origin cell
  float x = recHitTools.getPosition(originCellDetID).x();
  float y = recHitTools.getPosition(originCellDetID).y();
  std::pair<int, int> targetwaferortileid;
  std::pair<int, int> targetcellid;

  // get the position of the cell, that we are searching the neighbor for.
  if (isSiliconDet(originCellDetID.det()))
  {
    LOG(DEBUG) << "findCellCloseToXYpos:SiliconPos"<<"\n";
    targetwaferortileid = recHitTools.getWafer(originCellDetID);
    targetcellid = recHitTools.getCell(originCellDetID);
  }
  else
  {
    LOG(DEBUG) << "findCellCloseToXYpos:Scintillator"<<"\n";
    HGCScintillatorDetId scid = HGCScintillatorDetId(originCellDetID);
    targetwaferortileid = scid.ietaphi();
    targetcellid = std::make_pair(0, 0);
  }

  LOG(DEBUG) << "search from cell"<<"\n";
  LOG(DEBUG) << *getCellptr(originCellDetID)<<"\n";

  // The cell that is the closest, to be replace by close cells
  // if getcell can't find the wafer / cell it starts with the (0,0) coordinates
  DetId closest_cellDetId = getstartcell(targetdetectorid, targetsubdetid, targetlayerid, targetwaferortileid, targetcellid);
  Cell *closest_cellptr = getCellptr(closest_cellDetId);

  LOG(DEBUG) << "3: closest_cell"<<"\n";

  float x_cur = closest_cellptr->x;
  float y_cur = closest_cellptr->y;
  float d_cur = (x_cur - x) * (x_cur - x) + (y_cur - y) * (y_cur - y);
  bool improvement = true;
  LOG(DEBUG) << "3.5"<<"\n";
  while (improvement)
  {
    improvement = false;

    for (DetId &neighbor : closest_cellptr->neighbors)
    {
      LOG(DEBUG) << "    4"<<"\n";
      Cell *nptr = getCellptr(neighbor);
      LOG(DEBUG) << "\t" << getCellHashKeys(neighbor)<<"\n";
      LOG(DEBUG) << "\t" << *nptr<<"\n";
      float xn = nptr->x;
      float yn = nptr->y;
      float dn = (xn - x) * (xn - x) + (yn - y) * (yn - y);
      if (dn < d_cur)
      {
        LOG(DEBUG) << "\t"
                   << "improvement: " << dn << "<" << d_cur;
        improvement = true;
        closest_cellptr = nptr;
        d_cur = dn;
        break;
      }
    }
  }
  LOG(DEBUG) << "no improvement, exiting"<<"\n";
  return std::make_pair(closest_cellptr->globalid, d_cur);
}

DetId GeoExtractor::getstartcell(unsigned int detectorid, unsigned int subdetid, unsigned int layerid, std::pair<int, int> waferortileid, std::pair<int, int> cellid)
{
  if (detcol.detectors.find(detectorid) == detcol.detectors.end())
  {
    LOG(ERROR) << "No such detector:\n"
               << printcell(detectorid, subdetid, layerid, waferortileid, cellid);
    exit(EXIT_FAILURE);
  }
  Det &detector = detcol.detectors[detectorid];

  if (detector.subdetectors.find(subdetid) == detector.subdetectors.end())
  {
    LOG(ERROR) << "No such subdetector:\n"
               << printcell(detectorid, subdetid, layerid, waferortileid, cellid);
    exit(EXIT_FAILURE);
  }
  Subdet &subdet = detector.subdetectors[subdetid];

  if (subdet.layers.find(layerid) == subdet.layers.end())
  {
    LOG(ERROR) << "No such layer:\n"
               << printcell(detectorid, subdetid, layerid, waferortileid, cellid);
    exit(EXIT_FAILURE);
  }
  Layer &layer = subdet.layers[layerid];

  // for the EE and HSi part the cells are accessed via the wafers
  if (isSiliconDet(detectorid))
  {
    if (layer.wafers.begin() == layer.wafers.end())
    {
      LOG(ERROR) << "Wafer is empty:\n"
                 << printcell(detectorid, subdetid, layerid, waferortileid, cellid);
      exit(EXIT_FAILURE);
    }

    //find the wafer, it may not exist in this layer, in this case search from the first in the list.
    if (layer.wafers.find(waferortileid) == layer.wafers.end())
    {
      LOG(DEBUG) << "Cant find wafer, using first wafer."<<"\n";
      waferortileid = layer.wafers.begin()->first;
    }
    Wafer &wafer = layer.wafers[waferortileid];

    //Find the cell
    if (wafer.cells.find(cellid) == wafer.cells.end())
    {
      LOG(DEBUG) << "Cant find cell, using first."<<"\n";
      cellid = wafer.cells.begin()->first;
    }
    Cell &cell = wafer.cells[cellid];

    LOG(DEBUG) << "search: " << detectorid << " " << subdetid << " " << layerid << " " << waferortileid << " " << cellid<<"\n";

    LOG(DEBUG) << cell<<"\n";
    HGCSiliconDetId siid = HGCSiliconDetId(cell.globalid);
    LOG(DEBUG) << "startcell: "
               << " " << cell.globalid.det() << " " << cell.globalid.subdetId() << " " << recHitTools.getLayer(cell.globalid) << " " << siid.waferUV() << " " << siid.cellUV()<<"\n";

    return cell.globalid;
  }
  else
  {
    if (layer.tiles.begin() == layer.tiles.end())
    {
      LOG(ERROR) << "Tiles is empty\n"
                 << printcell(detectorid, subdetid, layerid, waferortileid, cellid);
      exit(EXIT_FAILURE);
    }
    if (layer.tiles.find(waferortileid) == layer.tiles.end())
    {
      waferortileid = layer.tiles.begin()->first;
    }
    return layer.tiles[waferortileid].globalid;
  }
  LOG(ERROR) << "Dead End1"<<"\n";
  exit(EXIT_FAILURE);
  return DetId(0);
}

std::string GeoExtractor::printcell(unsigned int detectorid, unsigned int subdetid, unsigned int layerid, std::pair<int, int> waferid, std::pair<int, int> cellid)
{
  std::ostringstream stringStream;
  stringStream << " Det " << detectorid;
  stringStream << " Subdet " << subdetid;
  stringStream << " Layer " << layerid;
  stringStream << " Wafer (" << waferid.first << "," << waferid.second << ")"<<"\n";
  stringStream << " Cell (" << cellid.first << "," << cellid.second << ")"<<"\n";
  std::string copyOfStr = stringStream.str();

  return copyOfStr;
}

bool GeoExtractor::validId(DetId id)
{
  auto [detectorid, subdetid, layerid, waferortileid, cellid] = getCellHashKeys(id);
  // unsigned int detectorid = id.det();
  // unsigned int subdetid = id.subdetId();
  // unsigned int layerid = recHitTools.getLayer(id);
  if (detectorid != DetId::HGCalEE && detectorid != DetId::HGCalHSi && detectorid != DetId::HGCalHSc)
  {
    LOG(DEBUG) << "wrong detector " << detectorid << "\n";
    return false;
  }
  if (subdetid != 0)
  {
    LOG(DEBUG) << "wrong subdetector " << subdetid << "\n";
    return false;
  }
  // if (detectorid != DetId::HGCalEE)
  // {
  //   if (layerid < 1 || layerid > 28)
  //   {
  //     LOG(DEBUG) << "layer ee "<< layerid << "\n";
  //     return false;
  //   }
  // }
  // if (detectorid != DetId::HGCalHSi)
  // {
  //   if (layerid < 1 || layerid > 22)
  //   {
  //     LOG(DEBUG) << "layer si "<< layerid << "\n";
  //     return false;
  //   }
  // }
  // if (detectorid != DetId::HGCalHSc)
  // {
  //   if (layerid < 9 || layerid > 22)
  //   {
  //     LOG(DEBUG) << "layer hsc "<< layerid << "\n";
  //     return false;
  //   }
  // }
  return true;
}
void GeoExtractor::validateId(DetId id)
{
  if (!validId(id))
  {
    LOG(DEBUG) << "Invalid cell Id:" << id.rawId()<<"\n";
    LOG(DEBUG) << *getCellptr(id)<<"\n";
    throw std::invalid_argument("Invalid Cell");
  }
}
bool GeoExtractor::isSiliconDet(int ndet)
{
  if (ndet == DetId::HGCalEE || ndet == DetId::HGCalHSi)
  {
    return true;
  }
  else if (ndet == DetId::HGCalHSc)
  {
    return false;
  }
  else
  {
    throw std::invalid_argument("Invalid detector");
    return false;
  }
}
// ------------ method called once each job just before starting event loop  ------------
void GeoExtractor::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void GeoExtractor::endJob() {}

//define this as a plug-in
DEFINE_FWK_MODULE(GeoExtractor);