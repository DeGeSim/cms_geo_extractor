#include <fstream>
#include <string>

class yamlwo
{
public:
  void printmembers(std::ofstream outfile, int indentlevel);

  void toyaml(std::ofstream outfile, int indentlevel = 0)
  {
    outfile << tabs(indentlevel) << "-{\n";
    indentlevel++;
    // print the members
    printmembers(outfile, indentlevel);
    // print the subsystem
    printmap(outfile, indentlevel);
    indentlevel--;
    outfile << tabs(indentlevel) << "}\n";
  }
};

class Cell : public yamlwo
{
public:
  float x;
  float y;
  DetId Id;
  std::vector<DetId> neighbors;
  void printmembers(std::ofstream outfile, int indentlevel)
  {
    outfile << tabs(indentlevel) << "x:" << (std::string)x << "\n";
    outfile << tabs(indentlevel) << "y:" << (std::string)x << "\n";
    outfile << tabs(indentlevel) << "Id:" << (std::string)Id.re << "\n";
  }
  void printmap(std::ofstream outfile, int indentlevel = 0)
  {
    for (auto const &e : neighbors)
    {
      outfile << tabs(indentlevel) << "-" << e.rawId() << "\n";
    }
  }
};
// class Tile: public yamlwo
// {
// public:
//   float middle_x;
//   float middle_y;
//   std::map<DetId, Cell> cells;

// };

class Wafer : public yamlwo
{
public:
  float middle_x;
  float middle_y;
  float si_thickness;
  std::map<DetId, Cell> cells;
  void printmembers(std::ofstream outfile, int indentlevel)
  {
    outfile << tabs(indentlevel) << "middle_x:" << (std::string)middle_x << "\n";
    outfile << tabs(indentlevel) << "middle_y:" << (std::string)middle_y << "\n";
    outfile << tabs(indentlevel) << "si_thickness:" << (std::string)si_thickness << "\n";
  }
  void printmap(std::ofstream outfile, int indentlevel = 0)
  {
    for (auto const &[key, val] : cells)
    {
      outfile << tabs(indentlevel) << "(" << key.first << "," << key.second << ") :\n";
      val.toyaml(outfile, indentlevel);
    }
  }
};

class Layer : public yamlwo
{
public:
  float z;
  std::map<std::pair<int, int>, Wafer> wafers;
  void printmembers(std::ofstream outfile, int indentlevel)
  {
    outfile << tabs(indentlevel) << "z:" << (std::string)z << "\n";
  }
  void printmap(std::ofstream outfile, int indentlevel = 0)
  {
    for (auto const &[key, val] : map)
    {
      outfile << tabs(indentlevel) << (std::string)key << ":\n";
      val.toyaml(outfile, indentlevel);
    }
  }

};

class Subdet : public yamlwo
{
public:
  std::map<int, Layer> layers;
  void printmembers(std::ofstream outfile, int indentlevel)
  {
    return;
  }
};

class Det : public yamlwo
{
public:
  std::map<int, Subdet> subdetectors;
  void printmembers(std::ofstream outfile, int indentlevel)
  {
    return;
  }
};

class DetColl : public yamlwo
{
public:
  std::map<int, Det> detectors;
  void printmembers(std::ofstream outfile, int indentlevel)
  {
    return;
  }
};

std::string tabs(int n)
{
  return std::string(n, '\t');
}
