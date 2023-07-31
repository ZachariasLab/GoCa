#include "chain.hpp"
#include "helper.hpp"

void Chain::add(const unsigned int& modelId, const Bead& bead) {
    if (mBeads.size() < modelId + 1) {
        mBeads.resize(modelId + 1);
    }
    mBeads[modelId].push_back(bead);
    for (auto atomIt = bead.atoms.begin(); atomIt != bead.atoms.end(); atomIt++){
        if (atomIt->x < minBoxCorner.x) minBoxCorner.x = atomIt->x;
        if (atomIt->x > maxBoxCorner.x) maxBoxCorner.x = atomIt->x;
        if (atomIt->y < minBoxCorner.y) minBoxCorner.y = atomIt->y;
        if (atomIt->y > maxBoxCorner.y) maxBoxCorner.y = atomIt->y;
        if (atomIt->z < minBoxCorner.z) minBoxCorner.z = atomIt->z;
        if (atomIt->z > maxBoxCorner.z) maxBoxCorner.z = atomIt->z;
    }
};

const bool Chain::valid() const {
    if (numberOfModels() == 0) {
        return false;
    } else if (numberOfModels() == 1) {
        return true;
    } else {
        for (unsigned int m = 1; m < numberOfModels(); m++) {
            if (mBeads[0].size() != mBeads[m].size()) {
                return false;
            }
            for (unsigned int b = 0; b < mBeads[m].size(); b++) {
                if (mBeads[0][b].residueName != mBeads[m][b].residueName) {
                    return false;
                }
            }
        }
    }
    return true;
};

const std::string Chain::getChainString() const {
    if (size() == 0) {
        return "";
    }
    std::stringstream ss;
    for (unsigned int i = 0; i < size(); i++) {
        ss << mBeads[0][i].residueName;
    }
    return ss.str();
};

void Chain::translateBeads(const Coordinates& vec) {
    for (auto modelIt = mBeads.begin(); modelIt != mBeads.end(); modelIt++) {
        for (auto chainIt = modelIt->begin(); chainIt != modelIt->end(); chainIt++) {
            chainIt->coordinates.translate(vec);
        }
    }
    minBoxCorner.translate(vec);
    maxBoxCorner.translate(vec);
};

const std::set<NativePair> Chain::getIntramolecularContacts(const Config& config) const {
    std::map<std::pair<int, int>, std::pair<int, double>> intramolecularContactBuilder;
    for (unsigned int modelId = 0; modelId < mBeads.size(); modelId++) {
        for (unsigned int beadId1 = 0; beadId1 < mBeads[modelId].size(); beadId1++) {
            for (unsigned int beadId2 = beadId1 + config.excludedNumber + 1; beadId2 < mBeads[modelId].size(); beadId2++) {
                if (mBeads[modelId][beadId1].isNativeContact(mBeads[modelId][beadId2], config)) {
                    double distance = sqrt(
                        pow2(mBeads[modelId][beadId1].coordinates.x - mBeads[modelId][beadId2].coordinates.x) +
                        pow2(mBeads[modelId][beadId1].coordinates.y - mBeads[modelId][beadId2].coordinates.y) +
                        pow2(mBeads[modelId][beadId1].coordinates.z - mBeads[modelId][beadId2].coordinates.z));
                    std::pair<int, int> key(beadId1, beadId2);
                    if (intramolecularContactBuilder.find(key) == intramolecularContactBuilder.end()) {
                        intramolecularContactBuilder.insert(std::make_pair(key, std::make_pair(1, distance)));
                    } else {
                        std::pair<int, double> dvalue = intramolecularContactBuilder[key];
                        intramolecularContactBuilder[key] = std::make_pair(dvalue.first + 1, dvalue.second + distance);
                    }
                }
            }
        }
    }
    std::set<NativePair> intramolecularContacts;
    for (auto it = intramolecularContactBuilder.begin(); it != intramolecularContactBuilder.end(); it++) {
        intramolecularContacts.emplace(it->first.first, it->first.second, it->second.second / it->second.first);
    }
    return intramolecularContacts;
};

const std::set<NativePair> Chain::getIntermolecularContacts(const Chain& other, const Config& config, const unsigned int& modelCount) const {
    std::set<NativePair> intermolecularContacts;
    std::map<std::pair<int, int>, std::pair<int, double>> intermolecularContactBuilder;
    double boxDistance = std::max(atomicVdwRadii.at('S') * 2 + config.atomicVdwDistance, config.atomicCutOff);
    if (this == &other || !boxOverlap(other.minBoxCorner, other.maxBoxCorner, boxDistance)) {
        /** Return an empty set if other chain is the same or if the other chain is too far away **/
        return intermolecularContacts;
    }
    for (unsigned int modelId = 0; modelId < modelCount; modelId++) {
        int ownModelId = (modelId < mBeads.size()) ? modelId : 0;
        int otherModelId = (modelId < other.mBeads.size()) ? modelId : 0;
        for (unsigned int beadId1 = 0; beadId1 != mBeads[ownModelId].size(); beadId1++) {
            for (unsigned int beadId2 = 0; beadId2 != other.mBeads[otherModelId].size(); beadId2++) {
                if (mBeads[ownModelId][beadId1].isNativeContact(other.mBeads[otherModelId][beadId2], config)) {
                    double distance = sqrt(
                        pow2(mBeads[ownModelId][beadId1].coordinates.x - other.mBeads[otherModelId][beadId2].coordinates.x) +
                        pow2(mBeads[ownModelId][beadId1].coordinates.y - other.mBeads[otherModelId][beadId2].coordinates.y) +
                        pow2(mBeads[ownModelId][beadId1].coordinates.z - other.mBeads[otherModelId][beadId2].coordinates.z));
                    std::pair<int, int> key(beadId1, beadId2);
                    if (intermolecularContactBuilder.find(key) == intermolecularContactBuilder.end()) {
                        intermolecularContactBuilder.insert(std::make_pair(key, std::make_pair(1, distance)));
                    } else {
                        std::pair<int, double> dvalue = intermolecularContactBuilder[key];
                        intermolecularContactBuilder[key] = std::make_pair(dvalue.first + 1, dvalue.second + distance);
                    }
                }
            }
        }
    }
    for (auto it = intermolecularContactBuilder.begin(); it != intermolecularContactBuilder.end(); it++) {
        intermolecularContacts.emplace(it->first.first, it->first.second, it->second.second / it->second.first);
    }
    return intermolecularContacts;
};

const bool Chain::boxOverlap(const Coordinates& min2, const Coordinates& max2, const double& boxPadding) const {
    const Coordinates min1 = minBoxCorner - Coordinates(boxPadding, boxPadding, boxPadding);
    const Coordinates max1 = maxBoxCorner + Coordinates(boxPadding, boxPadding, boxPadding);
    return ( /** x intersection **/
        (min1.x < min2.x && min2.x < max1.x) || (min1.x < max2.x && max2.x < max1.x) ||
        (min2.x < min1.x && min1.x < max2.x) || (min2.x < max1.x && max1.x < max2.x)
    ) && ( /** y intersection **/
        (min1.y < min2.y && min2.y < max1.y) || (min1.y < max2.y && max2.y < max1.y) ||
        (min2.y < min1.y && min1.y < max2.y) || (min2.y < max1.y && max1.y < max2.y)
    ) && ( /** z intersection **/
        (min1.z < min2.z && min2.z < max1.z) || (min1.z < max2.z && max2.z < max1.z) ||
        (min2.z < min1.z && min1.z < max2.z) || (min2.z < max1.z && max1.z < max2.z)
    );
};

void Chain::deleteAlternativeModels() {
    mBeads.erase(std::next(mBeads.begin()), mBeads.end());
};

void Chain::mergeDirectAdd(const Chain& other, const unsigned int& ownModelId, const unsigned int& otherModelId) {
    Eigen::Matrix3Xd fitted = rmsd(other, ownModelId, otherModelId).fittedCoordinates;
    for (unsigned int b = 0; b < size(); b++) {
        mBeads[ownModelId][b].coordinates += Coordinates(fitted(0, b), fitted(1, b), fitted(2, b));
    }
    auto temp = rmsd(*this, ownModelId, ownModelId);
    mCount += other.count();
};

void Chain::mergeDirectDivide(const unsigned int& ownModelId, const unsigned int& count) {
    for (unsigned int b = 0; b < size(); b++) {
        mBeads[ownModelId][b].coordinates /= count;
    }
};

void Chain::mergeAsModel(const Chain& other) {
    // We don't need to fit the models because only internal coordinates are relevant (and it might not be possible to do it accuratly)
    mBeads.insert(mBeads.end(), other.mBeads.begin(), other.mBeads.end());
    mCount += other.count();
};

const RmsdResult Chain::rmsd(const Chain& other, const int& modelId1, const int& modelId2) const {
    // Adapted from from: https://github.com/oleg-alexandrov/projects/blob/master/eigen/Kabsch.cpp
    // Info: http://en.wikipedia.org/wiki/Kabsch_algorithm
    if (size() != other.size()) {
        std::cerr << "Error: Both chains must have the same length for RMSD calculations (" << size() << " and " << other.size() << ")" << std::endl;
        exit(EXIT_FAILURE);
    }
    /** Get bead coordinates **/
    Eigen::Matrix3Xd m1(3, size());
    Eigen::Matrix3Xd m2(3, size());
    for (unsigned int i = 0; i < size(); i++) {
        m1(0, i) = mBeads[modelId1][i].coordinates.x / mCount;
        m1(1, i) = mBeads[modelId1][i].coordinates.y / mCount;
        m1(2, i) = mBeads[modelId1][i].coordinates.z / mCount;
        m2(0, i) = other.mBeads[modelId2][i].coordinates.x / other.mCount;
        m2(1, i) = other.mBeads[modelId2][i].coordinates.y / other.mCount;
        m2(2, i) = other.mBeads[modelId2][i].coordinates.z / other.mCount;
    }
    /** Transform matrix for testing purposes **/
    // Eigen::Matrix3d RTest = Eigen::Quaternion<double>(1, 3, 5, 2).normalized().toRotationMatrix();
    // Eigen::Vector3d TTest(-5, 6, -27);
    // double STest = 5;
    // for (int col = 0; col < m1.cols(); col++) {
    //     m1.col(col) = STest * (RTest * m1.col(col)) + TTest;
    // } 
    /** First find the scale, by finding the ratio of sums of some distances, then bring the datasets to the same scale. **/
    double dist1 = 0, dist2 = 0;
    for (int col = 0; col < m1.cols()-1; col++) {
        dist1  += (m1.col(col+1) - m1.col(col)).norm();
        dist2 += (m2.col(col+1) - m2.col(col)).norm();
    }
    double scale = dist2/dist1;
    m2 /= scale;
    /** Find the centroids then shift to the origin **/
    Eigen::Vector3d ctr1 = Eigen::Vector3d::Zero();
    Eigen::Vector3d ctr2 = Eigen::Vector3d::Zero();
    for (int col = 0; col < m1.cols(); col++) {
        ctr1  += m1.col(col);
        ctr2 += m2.col(col);
    }
    ctr1 /= m1.cols();
    ctr2 /= m2.cols();
    for (int col = 0; col < m1.cols(); col++) {
        m1.col(col)  -= ctr1;
        m2.col(col) -= ctr2;
    }
    /** Singular value decomposition **/
    Eigen::MatrixXd Cov = m1 * m2.transpose();
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(Cov, Eigen::ComputeThinU | Eigen::ComputeThinV);
    /** Find the rotation **/
    Eigen::Matrix3d R = svd.matrixV() * svd.matrixU().transpose();
    double d = (R.determinant() > 0) ? 1.0 : -1.0;
    Eigen::Matrix3d I = Eigen::Matrix3d::Identity(3, 3);
    I(2, 2) = d;
    R = svd.matrixV() * I * svd.matrixU().transpose();
    /* Calculate the RMSD */
    m2 = R.transpose() * m2;
    Eigen::Matrix3Xd diff = m2 - m1;
    double sum = 0;
    for (int col = 0; col < m2.cols(); col++) {
        sum += diff(0, col) * diff(0, col) + diff(1, col) * diff(1, col) + diff(2, col) * diff(2, col);
        m2.col(col) += ctr1;
    }
    double result = sqrt(sum / m2.cols());
    return RmsdResult(result, m2);
};
