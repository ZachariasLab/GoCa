#pragma once
#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <regex>
#include <set>

class Bead;

enum OutputType {
    GROMACS
};

struct Config {
    // General
    std::string inputPdb = "input.pdb";
    OutputType outputType = GROMACS;
    std::string coordinateFilename = "output.gro";
    std::string topologyFilename = "output.top";
    std::string moleculeName = "Mol";
    // Force field
    double epsilon = 1.0;
    double KBondLength = 20000.0;
    double KBondAngle = 40.0;
    double KDihedral1 = 1.0;
    double KDihedral2 = 0.5;
    double LJIntermolecular = 1.0;
    double LJRadius = 0.4;
    bool nonNativeAttraction = false;
    // Model details
    double atomicVdwDistance = 0.05;
    double beadVdwDistance = 0.0;
    bool includeBeadCutOff = false;
    double atomicCutOff = 0.0;
    double beadCutOff = 0.0;
    int excludedNumber = 3;
    double boxPadding = 3.0;
    double uniformMass = true;
    double angleDihedralCutOff = 155.0;
    bool varyingInteractionFactors = false;
    bool excludeHAtomsForCutoff = false;
    // Multi model
    std::vector<unsigned int> modelIds = {1};
    bool multiModel = false;
    std::string tabelDirectory = "tables";
    double minAngleDifference = 8.0;
    double minDihedralDifference = 8.0;
    unsigned int tabulatedPoints = 100;
    // Other
    bool saveConfig = false;
    bool saveNonbondedInfo = true;
    std::string splitTopologyFiles = "";
    double clusteringCutOff = 0;
    bool logChains = false;
    bool deleteAlternativeModels = false;

    Config(const std::string configFilename);
    std::string printString() const;
    bool isConfigured(const std::string& parameterName) const;

    private:
        void unknownParameterWarning(std::string group, std::string key, std::string value) const;
        std::string outputTypeString = "gromacs";
        bool isValid() const;
        std::set<std::string> configuredParameters;

        inline bool stringToBool(std::string& value) const {
            return (value == "true") || (value == "yes") || (value == "1") || (value == "True");
        };
        inline std::string boolToString(bool value) const {
            return ((value) ? "yes" : "no");
        };
        inline std::string valueToStringOrEmpty(bool condition, double value) const {
            return condition ? std::to_string(value) : "";
        }
};