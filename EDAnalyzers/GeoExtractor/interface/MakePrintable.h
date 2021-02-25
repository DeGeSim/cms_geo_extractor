#ifndef MakePrintable_H
#define MakePrintable_H 1
#include <iostream>
#include <vector>
using namespace std;

// C++ template to print vector container elements
template <typename T>
ostream &operator<<(ostream &os, const vector<T> &v)
{
    os << "[";
    for (int i = 0; i < v.size(); ++i)
    {
        os << v[i];
        if (i != v.size() - 1)
            os << ", ";
    }
    os << "]\n";
    return os;
}

// C++ template to print set container elements
template <typename T>
ostream &operator<<(ostream &os, const set<T> &v)
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
// C++ template to print map container elements
template <typename T, typename S>
ostream &operator<<(ostream &os, const map<T, S> &v)
{
    for (auto it : v)
        os << it.first << " : "
           << it.second << "\n";

    return os;
}

template <typename T, typename S>
ostream &operator<<(ostream &os, const pair<T, S> &v)
{
    os << "(";
    os << v.first << ", "
       << v.second << ")";
    return os;
}

typedef std::tuple<unsigned int, unsigned int, unsigned int, std::pair<int, int>, std::pair<int, int> > CellHash;

ostream &operator<<(ostream &os, const CellHash &hash)
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