#ifndef TreeOutputInfo_H
#define TreeOutputInfo_H 1

#include <TTree.h>

namespace TreeOutputInfo
{
  class TreeOutput
  {
  public:
    TTree *tree;

    std::vector<unsigned int> globalid;
    std::vector<unsigned int> detectorid;
    std::vector<unsigned int> subdetid;
    std::vector<unsigned int> layerid;
    std::vector<std::pair<int, int> > waferortileid;
    std::vector<std::pair<int, int> > cellid;

    std::vector<float> x;
    std::vector<float> y;
    std::vector<int> celltype;
    std::vector<bool> issilicon;
    std::vector<unsigned int> next;
    std::vector<unsigned int> previous;
    std::vector<unsigned int> nneighbors;
    std::vector<unsigned int> ngapneighbors;
    std::vector<unsigned int> n0;
    std::vector<unsigned int> n1;
    std::vector<unsigned int> n2;
    std::vector<unsigned int> n3;
    std::vector<unsigned int> n4;
    std::vector<unsigned int> n5;
    std::vector<unsigned int> n6;
    std::vector<unsigned int> n7;

    TreeOutput(std::string details, edm::Service<TFileService> fs)
    {
      tree = fs->make<TTree>(details.c_str(), details.c_str());

      // Run info //
      tree->Branch("globalid", &globalid);
      tree->Branch("detectorid", &detectorid);
      tree->Branch("subdetid", &subdetid);
      tree->Branch("layerid", &layerid);
      tree->Branch("waferortileid", &waferortileid);
      tree->Branch("cellid", &cellid);
      tree->Branch("x", &x);
      tree->Branch("y", &y);
      tree->Branch("celltype", &celltype);
      tree->Branch("issilicon", &issilicon);
      tree->Branch("next", &next);
      tree->Branch("previous", &previous);
      tree->Branch("nneighbors", &nneighbors);
      tree->Branch("ngapneighbors", &ngapneighbors);
      tree->Branch("n0", &n0);
      tree->Branch("n1", &n1);
      tree->Branch("n2", &n2);
      tree->Branch("n3", &n3);
      tree->Branch("n4", &n4);
      tree->Branch("n5", &n5);
      tree->Branch("n6", &n6);
      tree->Branch("n7", &n7);
    }

    void fill() { tree->Fill(); }

    void clear()
    {
      globalid.clear();
      cellid.clear();
      detectorid.clear();
      subdetid.clear();
      layerid.clear();
      waferortileid.clear();
      x.clear();
      y.clear();
      celltype.clear();
      issilicon.clear();
      next.clear();
      previous.clear();
      nneighbors.clear();
      ngapneighbors.clear();
      n0.clear();
      n1.clear();
      n2.clear();
      n3.clear();
      n3.clear();
      n4.clear();
      n5.clear();
      n6.clear();
      n7.clear();
    }
  };
}

#endif