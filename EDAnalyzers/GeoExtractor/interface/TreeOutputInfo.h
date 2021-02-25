#ifndef TreeOutputInfo_H
#define TreeOutputInfo_H 1

// #include <TH1F.h>
// #include <TH2F.h>
// #include <TMatrixD.h>
// #include <TROOT.h>
#include <TTree.h>
// #include <TVectorD.h>
// #include <stdlib.h>

// #include <iostream>
// #include <map>
// #include <string>
// #include <type_traits>
// #include <utility>
// #include <vector>

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
    std::vector<std::pair<int, int>> waferortileid;
    std::vector<std::pair<int, int>> cellid;

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
    }
  };
}

#endif