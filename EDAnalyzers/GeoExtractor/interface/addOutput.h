#ifndef additionalOutputInfo_H
#define additionalOutputInfo_H 1

#include <TTree.h>

namespace additionalOutputInfo
{
  class additionalOutput
  {
  public:
    TTree *tree;

    std::vector<unsigned int> originid;
    std::vector<unsigned int> targetid;
    std::vector<float> delta;

    additionalOutput(std::string details, edm::Service<TFileService> fs)
    {
      tree = fs->make<TTree>(details.c_str(), details.c_str());

      // Run info //
      tree->Branch("originid", &originid);
      tree->Branch("targetid", &targetid);
      tree->Branch("subdetid", &delta);

    }

    void fill() { tree->Fill(); }

    void clear()
    {
      originid.clear();
      targetid.clear();
      delta.clear();
    }
  };
}

#endif