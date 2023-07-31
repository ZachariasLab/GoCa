#include "config.hpp"
#include "helper.hpp"

Config::Config(const std::string configFilenameInput) {
    std::string configFilename = configFilenameInput;
    if (configFilenameInput.size() == 0) {
        configFilename = "GoCa.ini";
    }
    std::ifstream infile(configFilename);
    if (!infile.good()) {
        if (configFilenameInput.size() == 0) {
            return;
        } else {
            std::cerr << "Error: Configuration file \"" << configFilename << "\" does not exist" << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    std::string group;
    for(std::string line; getline(infile, line);) {
        std::string input = line.substr(0, std::min(line.find('#'), line.find(';')));
        if (input.length() == 0) {
            continue;
        }
        if (input.find("[") == 0) {
            group = input.substr(1, input.find("]") - 1);
        } else {
            int index = input.find("=");
            std::string key = trim(input.substr(0, index));
            std::string value = trim(input.substr(index + 1, input.length() - index - 1));
            configuredParameters.insert(key);
            if (group == "General") {
                if (key == "input") {
                    inputPdb = value;
                } else if (key == "type") {
                    outputTypeString = value;
                    if (value == "gromacs") {
                        outputType = GROMACS;
                    } else {
                        std::cerr << "Error: Output type \"" << value << "\" currently not supported" << std::endl;
                        exit(EXIT_FAILURE);
                    }
                }
                else if (key == "coordinates") coordinateFilename = value;
                else if (key == "topology") topologyFilename = value;
                else if (key == "name") moleculeName = value;
                else unknownParameterWarning(group, key, value);
            } else if (group == "Force field") {
                if (key == "epsilon") epsilon = std::stod(value);
                else if (key == "bond-length") KBondLength = std::stod(value);
                else if (key == "bond-angle") KBondAngle = std::stod(value);
                else if (key == "dihedral-1") KDihedral1 = std::stod(value);
                else if (key == "dihedral-2") KDihedral2 = std::stod(value);
                else if (key == "intermolecular") LJIntermolecular = std::stod(value);
                else if (key == "radius") LJRadius = std::stod(value);
                else if (key == "non-native-attraction") nonNativeAttraction = stringToBool(value);
                else unknownParameterWarning(group, key, value);
            } else if (group == "Model details") {
                if (key == "atomic-vdw-distance") atomicVdwDistance = std::stod(value);
                else if (key == "bead-vdw-distance") beadVdwDistance = std::stod(value);
                else if (key == "include-bead-cutoff") includeBeadCutOff = stringToBool(value);
                else if (key == "atomic-cutoff") atomicCutOff = std::stod(value);
                else if (key == "bead-cutoff") beadCutOff = std::stod(value);
                else if (key == "exclude-h-for-cutoff") excludeHAtomsForCutoff = stringToBool(value);
                else if (key == "excluded-number") excludedNumber = std::stoi(value);
                else if (key == "angle-dihedral-cutoff") angleDihedralCutOff = std::stod(value);
                else if (key == "box-padding") boxPadding = std::stod(value);
                else if (key == "uniform-mass") uniformMass = stringToBool(value);
                else if (key == "varying-interactions") varyingInteractionFactors = stringToBool(value);
                else unknownParameterWarning(group, key, value);
            } else if (group == "Multi model") {
                if (key == "model-indices") {
                    modelIds = splitToUInt(value, ',');
                    multiModel = modelIds.size() > 1;
                }
                else if (key == "table-directory") tabelDirectory = value;
                else if (key == "min-angle-dif") minAngleDifference = std::stod(value);
                else if (key == "min-dihedral-dif") minDihedralDifference = std::stoi(value);
                else if (key == "tabulation-points") tabulatedPoints = std::stoi(value);
                else unknownParameterWarning(group, key, value);
            } else if (group == "Other") {
                if (key == "save-config") saveConfig = stringToBool(value);
                else if (key == "save-nonbonded-info") saveNonbondedInfo = stringToBool(value);
                else if (key == "split-topology") splitTopologyFiles = value;
                else if (key == "cluster-cutoff") clusteringCutOff = std::stod(value);
                else if (key == "delete-other-models") deleteAlternativeModels = stringToBool(value);
                else if (key == "log-chains") logChains = stringToBool(value);
                else unknownParameterWarning(group, key, value);
            }
            else unknownParameterWarning(group, key, value);
        }
    }
    if (!isValid()) {
        exit(EXIT_FAILURE);
    }
};

bool Config::isValid() const {
    /** Warnings **/
    if (!includeBeadCutOff && beadCutOff > 0) {
        std::cout << "Warning: Parameter \"bead-cutoff\" is ignored when \"include-bead-cutoff\" is false" << std::endl;
    }
    /** Errors **/
    if (moleculeName.find(' ') != std::string::npos) {
        std::cerr << "Error: Moleculename must not contain any whitespace characters." << std::endl;
        return false;
    }
    return true;
};

void Config::unknownParameterWarning(std::string group, std::string key, std::string value) const {
    std::cout << "Warning: Unkown parameter \"" << key << "\" with value \"" << value << "\" in category \"" << group << "\"" << std::endl; 
};

bool Config::isConfigured(const std::string& parameterName) const {
    return configuredParameters.find(parameterName) != configuredParameters.end();
};

std::string Config::printString() const {
    std::stringstream ss;
    ss << ";   [General]\n";
    ss << ";   input                  = " << inputPdb << "\n";
    ss << ";   topology               = " << topologyFilename << "\n";
    ss << ";   coordinates            = " << coordinateFilename << "\n";
    ss << ";   type                   = " << outputTypeString << "\n";
    ss << ";   name                   = " << moleculeName << "\n";
    ss << ";\n;   [Force field]\n";
    ss << ";   epsilon                = " << epsilon << "\n";
    ss << ";   bond-length            = " << KBondLength << "\n";
    ss << ";   bond-angle             = " << KBondAngle << "\n";
    ss << ";   dihedral-1             = " << KDihedral1 << "\n";
    ss << ";   dihedral-2             = " << KDihedral2 << "\n";
    ss << ";   intermolecular         = " << LJIntermolecular << "\n";
    ss << ";   radius                 = " << LJRadius << "\n";
    ss << ";   non-native-attraction  = " << boolToString(nonNativeAttraction) << "\n";
    ss << ";\n;   [Model details]\n";
    ss << ";   box-padding            = " << boxPadding << "\n";
    ss << ";   uniform-mass           = " << boolToString(uniformMass) << "\n";
    ss << ";   atomic-vdw-distance    = " << atomicVdwDistance << "\n";
    ss << ";   bead-vdw-distance      = " << valueToStringOrEmpty(includeBeadCutOff, beadVdwDistance) << "\n";
    ss << ";   include-bead-cutoff   = " << boolToString(includeBeadCutOff) << "\n";
    ss << ";   atomic-cutoff         = " << valueToStringOrEmpty(atomicCutOff > 0, atomicCutOff) << "\n";
    ss << ";   bead-cutoff           = " << valueToStringOrEmpty(beadCutOff > 0, beadCutOff) << "\n";
    ss << ";   exclude-h-for-cutoff   = " << boolToString(excludeHAtomsForCutoff) << "\n";
    ss << ";   excluded-number        = " << excludedNumber << "\n";
    ss << ";   angle-dihedral-cutoff = " << angleDihedralCutOff << "\n";
    ss << ";   varying-interactions   = " << boolToString(varyingInteractionFactors) << "\n";
    ss << ";\n;   [Multi model]\n";
    ss << ";   model-indices          = " << vectorToString(modelIds, "%d", ", ") << "\n";
    ss << ";   table-directory        = " << tabelDirectory << "\n";
    ss << ";   min-angle-dif          = " << minAngleDifference << "\n";
    ss << ";   min-dihedral-dif       = " << minDihedralDifference << "\n";
    ss << ";   tabulation-points      = " << tabulatedPoints << "\n";
    ss << ";\n;   [Other]\n";
    ss << ";   save-config            = " << boolToString(saveConfig) << "\n";
    ss << ";   save-nonbonded-info    = " << boolToString(saveNonbondedInfo) << "\n";
    ss << ";   split-topology         = " << splitTopologyFiles << "\n";
    ss << ";   cluster-cutoff        = " << valueToStringOrEmpty(clusteringCutOff > 0, clusteringCutOff) << "\n";
    ss << ";   delete-other-models    = " << boolToString(deleteAlternativeModels) << "\n";
    ss << ";   log-chains             = " << boolToString(logChains) << "\n";
    return ss.str();
};