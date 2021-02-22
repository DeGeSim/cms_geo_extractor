#include <fstream>
#include <string>

std::string tabs(int n)
{
  return std::string(2*n, ' ');
}

class yamlwo
{
public:
  virtual void printmembers(std::ofstream &outfile, int indentlevel) { return; };
  virtual void printmap(std::ofstream &outfile, int indentlevel) { return; };
  // virtual ~yamlwo();
  void toyaml(std::ofstream &outfile, int indentlevel);
};

void yamlwo::toyaml(std::ofstream &outfile, int indentlevel = 0)
{
  // print the members
  printmembers(outfile, indentlevel);
  // print the subsystem
  printmap(outfile, indentlevel);
}

class Cell : public yamlwo
{
public:
  float x;
  float y;
  DetId globalid;
  DetId next;
  DetId previous;
  DetId north, south, west, east;
  std::vector<DetId> neighbors;
  void printmembers(std::ofstream &outfile, int indentlevel)
  {
    outfile << "x: " << x << "\n";
    // outfile << tabs(indentlevel) << "x: " << x << "\n";
    outfile << tabs(indentlevel) << "y: " << x << "\n";
    outfile << tabs(indentlevel) <<  "next: " << next.rawId() << "\n";
    // outfile << tabs(indentlevel) <<  "previous: " << previous.rawId() << "\n";
    // outfile << tabs(indentlevel) <<  "north: " << north.rawId() << "\n";
    // outfile << tabs(indentlevel) <<  "south: " << south.rawId() << "\n";
    // outfile << tabs(indentlevel) <<  "west: " << west.rawId() << "\n";
    // outfile << tabs(indentlevel) <<  "east: " << east.rawId() << "\n";
  }
  void printmap(std::ofstream &outfile, int indentlevel = 0)
  {
    outfile << tabs(indentlevel) << "neighbors: [";
    for (DetId const &e : neighbors)
    {
      outfile << e.rawId();
      if (e != neighbors[neighbors.size()-1])
      {
        outfile << ", ";
      }
    }
    outfile << "]\n";
  }
};

class Wafer : public yamlwo
{
public:
  float middle_x;
  float middle_y;
  int si_thickness;
  std::map<std::pair< int, int >, Cell> cells;
  void printmembers(std::ofstream &outfile, int indentlevel)
  {
    outfile << "middle_x: " << middle_x << "\n";
    outfile << tabs(indentlevel) << "middle_y: " << middle_y << "\n";
    outfile << tabs(indentlevel) << "si_thickness: " << si_thickness << "\n";
  }
  void printmap(std::ofstream &outfile, int indentlevel = 0)
  {
    for (auto &[key, val] : cells)
    {
      // outfile << tabs(indentlevel) << key.rawId() << ":\n";
      outfile << tabs(indentlevel) << "? !!python/tuple [" << key.first << ", "<< key.second << "]\n";
      outfile << tabs(indentlevel) << ": ";
      val.toyaml(outfile, indentlevel+1);
    }
  }
};

class Layer : public yamlwo
{
public:
  float z;
  std::map<std::pair<int, int>, Wafer> wafers;
  void printmembers(std::ofstream &outfile, int indentlevel)
  {
    outfile << tabs(indentlevel) << "z: " << z << "\n";
  }
  void printmap(std::ofstream &outfile, int indentlevel = 0)
  {
    for (auto &[key, val] : wafers)
    {
      outfile << tabs(indentlevel) << "? !!python/tuple [" << key.first << ", "<< key.second << "]\n";
      // outfile << tabs(indentlevel) << "? !!python/tuple\n";
      // outfile << tabs(indentlevel) << "- " << key.first << "\n";
      // outfile << tabs(indentlevel) << "- " << key.second << "\n";
      outfile << tabs(indentlevel) << ": ";
      val.toyaml(outfile, indentlevel+1);
    }
  }
};

class Subdet : public yamlwo
{
public:
  std::map<int, Layer> layers;
  void printmembers(std::ofstream &outfile, int indentlevel)
  {
    return;
  }
  void printmap(std::ofstream &outfile, int indentlevel)
  {
    for (auto &[key, val] : layers)
    {
      outfile << tabs(indentlevel) << key << ":\n";
      val.toyaml(outfile, indentlevel+1);
    }
  }
};

class Det : public yamlwo
{
public:
  std::map<int, Subdet> subdetectors;
  void printmembers(std::ofstream &outfile, int indentlevel)
  {
    return;
  }
  void printmap(std::ofstream &outfile, int indentlevel)
  {
    for (auto &[key, val] : subdetectors)
    {
      outfile << tabs(indentlevel) << key << ":\n";
      val.toyaml(outfile, indentlevel+1);
    }
  }
};

class DetColl : public yamlwo
{
public:
  std::map<int, Det> detectors;
  void printmembers(std::ofstream &outfile, int indentlevel)
  {
    return;
  }
  void printmap(std::ofstream &outfile, int indentlevel)
  {
    for (auto &[key, val] : detectors)
    {
      outfile << tabs(indentlevel) << key << ":\n";
      val.toyaml(outfile, indentlevel+1);
    }
  }
};
