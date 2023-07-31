#pragma once
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <sstream>
#include <math.h>
#include <algorithm>

#include "system.hpp"

/** String manipulation methods **/
static inline std::string trim(const std::string& str) {
    size_t endpos = str.find_last_not_of(' ');
    size_t startpos = str.find_first_not_of(' ');
    if (endpos == std::string::npos) {
        return "";
    }
    return str.substr(startpos, endpos - startpos + 1).substr();
};
static inline std::vector<unsigned int> splitToUInt(const std::string& text, const char delim) {
    std::string line;
    std::vector<unsigned int> vec;
    std::stringstream ss(text);
    while(std::getline(ss, line, delim)) {
        vec.push_back(std::stoi(line));
    }
    return vec;
};
static inline std::vector<double> splitToDouble(const std::string& text, const char delim) {
    std::string line;
    std::vector<double> vec;
    std::stringstream ss(text);
    while(std::getline(ss, line, delim)) {
        vec.push_back(std::stod(line));
    }
    return vec;
};
static inline std::vector<std::string> splitToString(const std::string& text) {
    std::stringstream ss(text);
    std::istream_iterator<std::string> begin(ss);
    std::istream_iterator<std::string> end;
    return std::vector<std::string>(begin, end);
};
template<typename ... Args>
static std::string string_format(const std::string& format, Args ... args) {
    int size_s = std::snprintf( nullptr, 0, format.c_str(), args ... ) + 1; // Extra space for '\0'
    if (size_s <= 0){
        throw std::runtime_error( "Error during formatting." );
    }
    auto size = static_cast<size_t>( size_s );
    std::unique_ptr<char[]> buf( new char[ size ] );
    std::snprintf( buf.get(), size, format.c_str(), args ... );
    return std::string( buf.get(), buf.get() + size - 1 ); // We don't want the '\0' inside
};
static bool similarStrings(const std::string& first, const std::string& second,
        const unsigned int windowSize=6, const unsigned int stride = 2, double matchValue=0) {
    unsigned int matches = 0;
    unsigned int count = 0;
    for (unsigned int i = 0; i < first.size() - windowSize; i += stride) {
        if (second.find(first.substr(i, windowSize)) != std::string::npos) {
            matches++;
        }
        count++;
    }
    if (matchValue <= 0) {
        matchValue = std::min(2 + (windowSize / stride) * 2, (unsigned int)first.size() / 5);
    }
    return matches + matchValue >= count;
}

/** Vector helper methods **/
static inline int findStringInVector(const std::vector<std::string>& vector, const std::string& s){
    auto it = std::find(vector.begin(), vector.end(), s);
    if (it == vector.end()) {
        return -1;
    } else {
        return std::distance(vector.begin(), it);
    }
};
template<typename T>
static inline std::string vectorToString(const std::vector<T>& vec, std::string formater, std::string delimiter = " ") {
    std::stringstream ss;
    for (auto it = vec.begin(); it != vec.end(); it++) {
        if (it != vec.begin()) {
            ss << delimiter;
        }
        ss << string_format(formater, *it);

    }
    return ss.str();
};

/** Math helper methods **/
static inline bool abs0or1(const int& value) {
    return value == 1 || value == -1 || value == 0;
};
static inline double pow2(const double& value) {
    return value * value;
};
static inline std::string intToBase26(int input) {
    std::stringstream ss;
    while (input / 26 > 0) {
        ss << static_cast<char>(65 + input % 26);
        input = input / 26 - 1;
    }
    ss << static_cast<char>(65 + input % 26);
    std::string result = ss.str();
    std::reverse(result.begin(), result.end());
    return result;
};

/** Operations on Coordinates objects which are effectively 3d-vectors **/
static inline double getDistance(const Coordinates& co1, const Coordinates& co2) {
    return sqrt(pow(co1.x - co2.x, 2) + pow(co1.y - co2.y, 2) + pow(co1.z - co2.z, 2));
};
static inline double vecNorm(const Coordinates& co) {
    return sqrt(pow(co.x, 2) + pow(co.y, 2) + pow(co.z, 2));
};
static inline double dotProduct(const Coordinates& co1, const Coordinates& co2) {
    return co1.x * co2.x + co1.y * co2.y + co1.z * co2.z;
};
static inline Coordinates crossProduct(const Coordinates& co1, const Coordinates& co2) {
    return Coordinates(co1.y * co2.z - co1.z * co2.y, co1.z * co2.x - co1.x * co2.z, co1.x * co2.y - co1.y * co2.x);
};

/** Angle and dihedral calculation **/
static double getAngle(const Coordinates& co1, const Coordinates& co2, const Coordinates& co3) {
    Coordinates d1 = co1 - co2;
    Coordinates d2 = co3 - co2;
    double d = dotProduct(d1, d2) / vecNorm(d1) / vecNorm(d2);
    d = (d > 1.0) ? 1.0 : (d < -1) ? -1 : d;
    return acos(d) * 180.0 / M_PI; // Returns degrees
};
static double getDihedrals(const Coordinates& co1, const Coordinates& co2, const Coordinates& co3, const Coordinates& co4) {
    Coordinates b2 = co3 - co2;
    Coordinates n1 = crossProduct(co2 - co1, b2);
    Coordinates n2 = crossProduct(b2, co4 - co3);
    Coordinates m1 = crossProduct(n1 / vecNorm(n1), b2 / vecNorm(b2));
    double x = dotProduct(n1 / vecNorm(n1), n2 / vecNorm(n2));
    double y = dotProduct(m1, n2 / vecNorm(n2));
    double result = (-atan2(y, x) * 180.0 / M_PI);
    return result;  // Returns degrees
};

/** Debug message printing methods **/
static void print2dRMSD(std::vector<Chain>& chains) {
    std::cout << "RMSD:" << std::endl;
    for (unsigned int i = 0; i < chains.size(); i++) {
        for (unsigned int j = 0; j < chains.size(); j++) {
            double rmsd  = chains[i].rmsd(chains[j], 0, 0).rmsd;
            std::cout << std::setprecision(2) << std::setw(8) << rmsd << " ";
        }
        std::cout << std::setprecision(0) << std::setw(0) << std::endl;
    }
};
static void printClusteringInfo(std::map<std::string, std::vector<std::vector<std::pair<unsigned int, unsigned int>>>>& clusters) {
    std::cout << "Info: Chain clustering results:" << std::endl;
    for (auto clIt = clusters.begin(); clIt != clusters.end(); clIt++) {
        for (auto chIt = clIt->second.begin(); chIt != clIt->second.end(); chIt++) {
            if (chIt != clIt->second.begin()) {
                std::cout << ", ";
            }
            std::cout << chIt->size();
        }
        std::cout << " of size " << clIt->first.size() << std::endl;
    }
    std::cout << std::endl;
};

/** Other methods **/
template<typename T>
static inline std::pair<T, T> makeOrderedPair(const T& first, const T& second) {
    if (first < second) {
        return std::make_pair(first, second);
    }
    return std::make_pair(second, first);
};
