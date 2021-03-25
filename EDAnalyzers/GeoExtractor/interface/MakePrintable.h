#ifndef MakePrintable_H
#define MakePrintable_H 1
#include <iostream>
#include <vector>

// C++ template to print std::vector container elements
template <typename T>
std::ostream &operator<<(std::ostream &os, const std::vector<T> &v)
{
    os << "[";
    for (int i = 0; i < (int)v.size(); ++i)
    {
        os << v[i];
        if (i != (int)v.size() - 1)
            os << ", ";
    }
    os << "]\n";
    return os;
}

// C++ template to print std::set container elements
template <typename T>
std::ostream &operator<<(std::ostream &os, const std::set<T> &v)
{
    os << "[";
    for (auto it : v)
    {
        os << it;
        if (it != *v.rbegin())
            os << ", ";
    }
    os << "]\n";
    return os;
}
// C++ template to print std::list container elements
template <typename T>
std::ostream &operator<<(std::ostream &os, const std::list<T> &v)
{
    os << "[";
    for (auto it : v)
    {
        os << it;
        if (it != *v.rbegin())
            os << ", ";
    }
    os << "]\n";
    return os;
}
// C++ template to print std::map container elements
template <typename T, typename S>
std::ostream &operator<<(std::ostream &os, const std::map<T, S> &v)
{
    for (auto it : v)
        os << it.first << " : "
           << it.second << "\n";

    return os;
}

template <typename T, typename S>
std::ostream &operator<<(std::ostream &os, const std::pair<T, S> &v)
{
    os << "(";
    os << v.first << ", "
       << v.second << ")";
    return os;
}

typedef std::tuple<unsigned int, unsigned int, unsigned int, std::pair<int, int>, std::pair<int, int> > CellHash;

std::ostream &operator<<(std::ostream &os, const CellHash &hash)
{
    auto [detectorid, subdetid, layerid, waferortileid, cellid] = hash;
    os << "Det " << detectorid;
    os << " Subdet " << subdetid;
    os << " Layer " << layerid;
    os << " Wafer (" << waferortileid.first << "," << waferortileid.second << ")";
    os << " Cell (" << cellid.first << "," << cellid.second << ")";
    return os;
}

typedef std::tuple<double, Cell *> PosListTup;
std::ostream &operator<<(std::ostream &os, const PosListTup &tup)
{
    auto [xpos, cellptr] = tup;
    os << "(x: " << xpos << " id: " << cellptr->globalid.rawId() << ") ";
    return os;
}

bool operator<(const PosListTup &a, const PosListTup &b)
{
    return (std::get<0>(a) < std::get<0>(b));
}

std::ostringstream &operator<<(std::ostringstream &os, const CellHash &hash)
{
    auto [detectorid, subdetid, layerid, waferortileid, cellid] = hash;
    os << "Det " << detectorid;
    os << " Subdet " << subdetid;
    os << " Layer " << layerid;
    os << " Wafer (" << waferortileid.first << "," << waferortileid.second << ")";
    os << " Cell (" << cellid.first << "," << cellid.second << ")";
    return os;
}

#endif