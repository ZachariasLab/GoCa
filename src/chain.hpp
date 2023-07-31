#pragma once

#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <limits>
#include <math.h>
#include <algorithm>
#include <cstdlib>
#include <Eigen/Geometry>

#include "const.hpp"
#include "config.hpp"
#include "bead.hpp"

struct RmsdResult {
    double rmsd;
    Eigen::Matrix3Xd fittedCoordinates;

    RmsdResult(const double& rmsd, Eigen::Matrix3Xd& fittedCoordinates) : rmsd(rmsd), fittedCoordinates(fittedCoordinates) {};
};

class Chain {
    friend class System;

private:
    std::vector<std::vector<Bead>> mBeads;
    unsigned int mCount = 1;
    Coordinates minBoxCorner = Coordinates(std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max());
    Coordinates maxBoxCorner = Coordinates(std::numeric_limits<double>::min(), std::numeric_limits<double>::min(), std::numeric_limits<double>::min());

    Chain() {};

    const bool boxOverlap(const Coordinates& minCornerOther, const Coordinates& maxCornerOther, const double& boxPadding = 0.0) const;
    const Eigen::Matrix3Xd combineCoordinates(const Chain& other);
    
public:
    const std::vector<std::vector<Bead>>& beads() const {return mBeads;};
    const unsigned int& count() const {return mCount;};

    void add(const unsigned int& modelId, const Bead& bead);
    void translateBeads(const Coordinates& vec);
    void mergeDirectAdd(const Chain& other, const unsigned int& ownModelId, const unsigned int& otherModelId);
    void mergeDirectDivide(const unsigned int& ownModelId, const unsigned int& count);
    void mergeAsModel(const Chain& other);
    void deleteAlternativeModels();
    const std::set<NativePair> getIntramolecularContacts(const Config& config) const;
    const std::set<NativePair> getIntermolecularContacts(const Chain& other, const Config& config, const unsigned int& modelCount) const;

    const bool valid() const;
    const std::string getChainString() const;
    const RmsdResult rmsd(const Chain& other, const int& modelId1, const int& modelId2) const;

    const inline size_t numberOfModels() const { return mBeads.size(); };
    const inline size_t size() const { if (mBeads.size() == 0) { return 0; } return mBeads[0].size(); };
};