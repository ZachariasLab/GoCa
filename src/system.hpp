#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <tuple>
#include <map>
#include <limits>
#include <math.h>

#include "const.hpp"
#include "config.hpp"
#include "chain.hpp"

struct CellIndex {
    int x, y, z;

    CellIndex(const int& x, const int& y, const int& z) : x(x), y(y), z(z) {};
    inline bool operator=(CellIndex other) const {return x == other.x && y == other.y && z == other.z;};
    inline bool operator<(CellIndex other) const {return (x < other.x) || (x == other.x && y < other.y) || (x == other.x && y == other.y && z < other.z);}
};

enum FileType {
    PDB, CIF
};

struct AtomEntry {
    std::string chainName;
    std::string residueName;
    std::string atomName;
    unsigned int residueId;
    unsigned int modelId;
    double x;
    double y;
    double z;
    bool altAtom;

    AtomEntry(const std::string& line, const FileType fileType);
};

class System {
public:
    std::set<char> atomTypes;
    std::vector<Chain> chains;
    unsigned int modelCount;

    System(Config& config);
    void translateBeads(const Coordinates& vec);
    const Coordinates getBoxSize(const double& boxPadding);
    const Coordinates getBoxTranslationVector(const double& boxPadding);
    const std::map<std::string, std::vector<std::vector<std::pair<unsigned int, unsigned int>>>> getChainClusters(const double& clusterCutOff);
    std::map<unsigned int, unsigned int> mergeChainsAndIntermolecularContacts(
        const std::map<std::string, std::vector<std::vector<std::pair<unsigned int, unsigned int>>>>& clusters,
        std::map<std::pair<int, int>, std::set<NativePair>>& intermolecularNativeContacts,
        bool varyingInteractionFactors
    );
    const int size();

private:
    Coordinates minBoxCorner = Coordinates(std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max());
    Coordinates maxBoxCorner = Coordinates(std::numeric_limits<double>::min(), std::numeric_limits<double>::min(), std::numeric_limits<double>::min());
};