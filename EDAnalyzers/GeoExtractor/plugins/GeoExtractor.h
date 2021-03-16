#pragma once

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

  //The map, that contains the detector structure
  //Det -> SubDet -> Layer -> Wafer -> Cell
  DetColl detcol;
  //function to the data strution containing the detector/subdetector/wafers/cells
  void instanciateMapForCell(DetId &iterId);

  //vector to save the ids
  std::vector<DetId> v_validHGCalIds;
  //function to filter for cells in the HGCAL
  std::vector<DetId> filterCellIds(const std::vector<DetId> v_allCellIds);

  //Funtions for searching the neighbors in z direction
  void assignZneighbors(std::vector<DetId> &v_validHGCalIds);
  DetId findNextCell(DetId cellId);
  DetId findPreviousCell(DetId cellId);
  std::pair<DetId, float> searchInLayer(DetId cellId, CellHash hash, unsigned int detectorid, unsigned int subdetid, unsigned int layerid);
  DetId getstartcell(unsigned int detectorid, unsigned int subdetid, unsigned int layerid, std::pair<int, int> wafer, std::pair<int, int> cell);

  void fixNeighborsBoundery(std::vector<DetId> &v_validHGCalIds);
  DetId findGapNeighbors(Cell *cellptr);

  //Get the det/subdet/wafer/cell id as a tuple
  CellHash getCellHashKeys(DetId &iterId);

  std::string printcell(unsigned int detectorid, unsigned int subdetid, unsigned int layerid, std::pair<int, int> wafer, std::pair<int, int> cell);

  Cell *getCellptr(DetId &iterId);
  bool validId(DetId id);
  void validateId(DetId id);
  bool isSiliconDet(int ndet);

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