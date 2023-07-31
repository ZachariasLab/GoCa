#pragma once

#include <vector>
#include <iostream>
#include <math.h>

#include "const.hpp"
#include "config.hpp"

struct Coordinates {
    double x, y, z;

    Coordinates(){};
    Coordinates(const double& x, const double& y, const double& z) : x(x), y(y), z(z) {};
    Coordinates(const double& x, const double& y, const double& z, const double& factor) : x(x/factor), y(y/factor), z(z/factor) {};
    inline void translate(const Coordinates& vec) {x += vec.x; y += vec.y; z += vec.z;};
    inline Coordinates operator-(const Coordinates& other) const {return Coordinates(x - other.x, y - other.y, z - other.z);};
    inline Coordinates operator+(const Coordinates& other) const {return Coordinates(x + other.x, y + other.y, z + other.z);};
    inline Coordinates operator/(const double& divisor) const {return Coordinates(x / divisor, y / divisor, z / divisor);};
    inline Coordinates min(const Coordinates& other) const {return Coordinates(std::min(x, other.x), std::min(y, other.y), std::min(z, other.z));};
    inline Coordinates max(const Coordinates& other) const {return Coordinates(std::max(x, other.x), std::max(y, other.y), std::max(z, other.z));};
    inline Coordinates& operator+=(const Coordinates& other) {x += other.x; y += other.y; z += other.z; return *this;};
    inline Coordinates& operator-=(const Coordinates& other) {x -= other.x; y -= other.y; z -= other.z; return *this;};
    inline Coordinates& operator/=(const double& divisor) {x /= divisor; y /= divisor; z /= divisor; return *this;};
    inline double norm() const {return sqrt(x * x + y * y + z * z);};

    friend std::ostream& operator<<(std::ostream& os, const Coordinates& coordinates) { 
        return os << "(" << coordinates.x << " " << coordinates.y << " " << coordinates.z << ")";
    }
};

struct Bead {
    unsigned int modelId = -1;
    unsigned int chainId;
    unsigned int residueId;
    char residueName;
    Coordinates coordinates;
    double mass;
    std::vector<Coordinates> atoms;
    std::vector<std::string> atomsTypes;

    Bead(){};
    const bool isNativeContact(const Bead& other, const Config& config) const;
    const ContactType getContactType(const Bead& other) const;
private:
    const bool isAnion(const unsigned int& index) const;
    const bool isCation(const unsigned int& index) const;
    const unsigned int getClostestResidueAtom(const unsigned int& index) const;
    const bool hydrophobicResidue() const;
};

struct NativePair {
    int i, j;
    double distance;
    mutable int id;
    mutable double factor;

    NativePair(const int& i, const int& j, const double& distance, const int& id = 0, const double& factor=1) : i(i), j(j), distance(distance), id(id), factor(factor) {};
    inline bool operator<(const NativePair& other) const { return (i == other.i) ? j < other.j : i < other.i; };
};