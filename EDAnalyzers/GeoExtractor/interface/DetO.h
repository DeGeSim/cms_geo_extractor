#ifndef DetO_Hcust
#define DetO_Hcust 1
#include <fstream>
#include <string>
#include <iostream>

std::string tabs(int n)
{
  return std::string(2 * n, ' ');
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
  bool issilicon;
  int type;
  float x;
  float y;
  DetId globalid;
  DetId next;
  DetId previous;
  std::vector<DetId> neighbors;
  void printmembers(std::ofstream &outfile, int indentlevel)
  {
    outfile << "issilicon: " << issilicon << "\n";
    if (issilicon)
    {
      outfile << tabs(indentlevel) << "type: " << type << "\n";
    }
    outfile << tabs(indentlevel) << "x: " << x << "\n";
    outfile << tabs(indentlevel) << "y: " << x << "\n";
    outfile << tabs(indentlevel) << "next: " << next.rawId() << "\n";
  }
  void printmap(std::ofstream &outfile, int indentlevel = 0)
  {
    outfile << tabs(indentlevel) << "neighbors: [";
    for (DetId const &e : neighbors)
    {
      outfile << e.rawId();
      if (e != neighbors[neighbors.size() - 1])
      {
        outfile << ", ";
      }
    }
    outfile << "]\n";
  }
  friend std::ostream &operator<<(std::ostream &os, const Cell &c)
  {
    os << "Cell {";
    os << "issilicon: " << c.issilicon << ", ";
    os << "type: " << c.type << ", ";
    os << "globalid: " << c.globalid.rawId() << ", ";
    os << "x: " << c.x << ", ";
    os << "y: " << c.y << ", ";
    os << "next: " << c.next.rawId() << ", ";
    os << "previous: " << c.previous.rawId() << ", ";
    os << "neighbors: [";
    for (auto &val : c.neighbors)
    {
      os << val.rawId() << ", ";
    }
    os << "]}\n";
    return os;
  }
};

class Wafer : public yamlwo
{
public:
  bool issilicon = true;
  int si_thickness;
  std::map<std::pair<int, int>, Cell> cells;
  void printmembers(std::ofstream &outfile, int indentlevel)
  {
    outfile << "issilicon: " << issilicon << "\n";
    outfile << tabs(indentlevel) << "si_thickness: " << si_thickness << "\n";
  }
  void printmap(std::ofstream &outfile, int indentlevel = 0)
  {
    for (auto &[key, val] : cells)
    {
      // outfile << tabs(indentlevel) << key.rawId() << ":\n";
      outfile << tabs(indentlevel) << "? !!python/tuple [" << key.first << ", " << key.second << "]\n";
      outfile << tabs(indentlevel) << ": ";
      val.toyaml(outfile, indentlevel + 1);
    }
  }
};

class Layer : public yamlwo
{
public:
  float z;
  std::map<std::pair<int, int>, Wafer> wafers;
  std::map<std::pair<int, int>, Cell> tiles;
  void printmembers(std::ofstream &outfile, int indentlevel)
  {
    outfile << tabs(indentlevel) << "z: " << z << "\n";
  }
  void printmap(std::ofstream &outfile, int indentlevel = 0)
  {
    for (auto &[key, val] : wafers)
    {
      outfile << tabs(indentlevel) << "? !!python/tuple [" << key.first << ", " << key.second << "]\n";
      outfile << tabs(indentlevel) << ": ";
      val.toyaml(outfile, indentlevel + 1);
    }
    for (auto &[key, val] : tiles)
    {
      outfile << tabs(indentlevel) << "? !!python/tuple [" << key.first << ", " << key.second << "]\n";
      outfile << tabs(indentlevel) << ": ";
      val.toyaml(outfile, indentlevel + 1);
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
      val.toyaml(outfile, indentlevel + 1);
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
      val.toyaml(outfile, indentlevel + 1);
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
      val.toyaml(outfile, indentlevel + 1);
    }
  }
};

#endif