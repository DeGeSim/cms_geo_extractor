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

  //// IO
  // file to write the yaml structured detector information out
  std::ofstream yamlfile;
  edm::Service<TFileService> fs;

  // Debug level
  enum typelog internalDebugLevel;
  // output tree
  TreeOutputInfo::TreeOutput *treeOutput;
  // additionalOutputInfo::additionalOutput *addOutput;

  //Tools
  edm::ESHandle<CaloGeometry> geom;
  hgcal::RecHitTools recHitTools;

  // container for the topologyies
  std::map<int, edm::ESHandle<HGCalTopology>> m_topo;
  std::map<int, edm::ESHandle<HGCalGeometry>> m_geom;

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
  void assignZNeighbors(std::vector<DetId> &v_validHGCalIds);
  DetId findNextCell(DetId cellId);
  DetId findPreviousCell(DetId cellId);
  std::pair<DetId, float> searchInLayer(DetId cellId, CellHash hash, unsigned int detectorid, unsigned int subdetid, int layerid);
  DetId getStartCell(CellHash hash);

  // For the gapfixing:
  void fixGap(std::vector<DetId> &v_validHGCalIds);
  // detector subdet layer -> list [ (xpos, cellptr) ]
  std::map<unsigned int, std::map<unsigned int, std::map<int, std::vector<PosListTup>>>> xdistmap;
  static bool sort_xdistmap_vectors(const PosListTup &first, const PosListTup &second);
  void assignGapNeighbors(Cell *cellptr);
  void setupXLists();
  double cellsDelta(Cell *cp1, Cell *cp2);
  bool rangecond(std::vector<PosListTup>::iterator iter, Cell *cellptr);
  double dist(GlobalPoint &p1, GlobalPoint &p2);
  double xposdiff(std::vector<PosListTup>::iterator iter, Cell *cellptr);
  
  ////Utils
  //Get the det/subdet/wafer/cell id as a tuple
  CellHash getCellHash(const DetId &iterId);

  std::string printCell(unsigned int detectorid, unsigned int subdetid, int layerid, std::pair<int, int> wafer, std::pair<int, int> cell);
  int realLayerFromId (const DetId &iterId);
  Cell *getCellPtr(const DetId &iterId);
  bool validId(DetId id);
  void validateId(DetId id);
  bool isSiliconDet(int ndet);

  //vector with the numbers of the detector part of the hgcal
  std::vector<int> v_HGCalDets;
  //map that saves which cell are rejected in which step
  std::map<int, std::map<std::string, int>> m_rej;

  //Variables for the parameters to be passed
  double maxSearchDelta;
  double maxDeltaHScHSiGap;
  bool filter_cylinder;
  double selectionRadius;
  double selection_x;
  double selection_y;
};

GeoExtractor::GeoExtractor(const edm::ParameterSet &iConfig)
{
  //
  maxSearchDelta = iConfig.getParameter<double>("maxSearchDelta");
  maxDeltaHScHSiGap = iConfig.getParameter<double>("maxDeltaHScHSiGap");
  selectionRadius = iConfig.getParameter<double>("selectionRadius");
  selection_y = iConfig.getParameter<double>("selection_y");
  selection_x = iConfig.getParameter<double>("selection_x");
  filter_cylinder = iConfig.getParameter<bool>("filter_cylinder");

  usesResource("TFileService");
  treeOutput = new TreeOutputInfo::TreeOutput("tree", fs);

  m_topo[DetId::HGCalEE];
  m_topo[DetId::HGCalHSi];
  m_topo[DetId::HGCalHSc];

  v_HGCalDets.push_back(DetId::HGCalEE);
  v_HGCalDets.push_back(DetId::HGCalHSi);
  v_HGCalDets.push_back(DetId::HGCalHSc);

  //Set the loglevel here DEBUG(0) < INFO(1) < WARN(2) < ERROR(3)
  //Some typecasting to get the integer from cmssw to a the desired enum type.
  internalDebugLevel = static_cast<typelog>(iConfig.getParameter<int>("internalDebugLevel"));
  std::cout << "Starting with loglevel " << internalDebugLevel << "\n";
  LOGCFG.level = internalDebugLevel;
  LOGCFG.headers = false;
}
GeoExtractor::~GeoExtractor()
{
  treeOutput->fill();
  // addOutput->fill();
  yamlfile.open("output/geometry.yaml");
  yamlfile.clear();
  detcol.toyaml(yamlfile, 0);
  yamlfile.close();
}