#include <iostream>
#include <fstream>
#include <regex>
#include <algorithm>

#include "const.hpp"
#include "chain.hpp"
#include "config.hpp"
#include "helper.hpp"
#include "files.hpp"
#include "system.hpp"

int main(int argc, char *argv[]) {

    /** Read config file **/
    const std::string configFilename = (argc > 1) ? argv[1] : "";
    Config config(configFilename);

    std::cout << "1. Generating coarse grained model" << std::endl;
    
    /** Get bead system from protein structure **/
    System system(config);

    /** Get periodic box dimensions and translate beads to box center **/
    Coordinates boxSize = system.getBoxSize(config.boxPadding);
    system.translateBeads(system.getBoxTranslationVector(config.boxPadding));

    /** Create a copy of the chains before chain merging **/
    std::vector<Chain> originalChains = system.chains;
    unsigned int systemSize = system.size();

    std::cout << "2. Calculating native contacts" << std::endl;

    /** Calculate intermolecular native contacts **/
    std::map<std::pair<int, int>, std::set<NativePair>> intermolecularNativeContacts;
    std::vector<std::pair<unsigned int, unsigned int>> pairIndices;
    for (unsigned int mi = 0; mi < system.chains.size(); mi++) {
        for (unsigned int mj = mi + 1; mj < system.chains.size(); mj++) {
            pairIndices.push_back(std::pair<unsigned int, unsigned int>(mi, mj));
        }
    }
    #pragma omp parallel for
    for (int i = 0; i < pairIndices.size(); i++) {
        std::set<NativePair> chainPairNativePairs = system.chains[pairIndices[i].first].getIntermolecularContacts(system.chains[pairIndices[i].second], config, system.modelCount);
        if (chainPairNativePairs.size() > 0) {
            #pragma omp critical
            intermolecularNativeContacts.insert(std::make_pair(std::pair<int, int>(pairIndices[i].first, pairIndices[i].second), chainPairNativePairs));
        }
    }

    /* Perform chain merging if clustering cutoff is provided and we have more than one chain **/
    std::map<unsigned int, unsigned int> chainIndicesMap;
    if (config.clusteringCutOff > 0 && system.chains.size() > 1) {
        const std::map<std::string, std::vector<std::vector<std::pair<unsigned int, unsigned int>>>> clusters = system.getChainClusters(config.clusteringCutOff);
        /** Print clustering results **/
        if (clusters.size() == system.chains.size()) {
            std::cout << "Warning: Chain merging is enabled but no chains with an identical sequence were found" << std::endl;
        } else if (clusters.size() == 1) {
            std::cout << "Info: Got 1 merged chain from originally " << system.chains.size() << " chains" << std::endl;
        } else {
            std::stringstream sizeInfo;
            for (auto clit = clusters.begin(); clit != clusters.end(); clit++) {
                if (clit == clusters.begin()) {
                    sizeInfo << "Sizes: ";
                } else {
                    sizeInfo << " ";
                }
                sizeInfo << clit->first.size();
            }
            std::cout << "Info: Got " << clusters.size() << " merged chains (" << sizeInfo.str();
            std::cout << " residues) from originally " << system.chains.size() << " chains" << std::endl;
        }
        /** Merge clustered chains **/
        chainIndicesMap = system.mergeChainsAndIntermolecularContacts(clusters, intermolecularNativeContacts, config.varyingInteractionFactors);
        /** Warn if chains are similar but not identical **/
        unsigned int similarChains = 0;
        for (auto it1 = clusters.begin(); it1 != clusters.end(); it1++) {
            for (auto it2 = std::next(it1); it2 != clusters.end(); it2++) {
                if (similarStrings(it1->first, it2->first)) {
                    similarChains++;
                }
            }
        }
        if (similarChains == 1) {
            std::cout << "Warning: There is one pair of chains which are very similar but not identical." << std::endl;
        } else if (similarChains > 1) {
            std::cout << "Warning: There are " << similarChains << " pairs of chains which are very similar but not identical." << std::endl;
        }
    }

    if (config.logChains) {
        for (unsigned int i = 0; i < system.chains.size(); i++) {
            std::cout << "Chain " << (i+1) << ": " << system.chains[i].getChainString() << std::endl;
        }
    }

    if (config.deleteAlternativeModels) {
        for (auto cIt = system.chains.begin(); cIt != system.chains.end(); cIt++) {
            cIt->deleteAlternativeModels();
        }
    }

    if (config.varyingInteractionFactors && intermolecularNativeContacts.size() > 1) {
        std::cout << "Intermolecular interaction modification factors between different chains:" << std::endl;
        for (auto contactIt = intermolecularNativeContacts.begin(); contactIt != intermolecularNativeContacts.end(); contactIt++) {
            if (contactIt->first.first != contactIt->first.second) {
                std::cout << "  Please enter a modification factor for interactions between the chains with ";
                std::cout << system.chains[contactIt->first.first].size() << " and " << system.chains[contactIt->first.second].size() << " residues";
                std::cout << " (" << contactIt->second.size() << " contacts): ";
                std::string input;
                std::getline(std::cin, input);
                double factor;
                try {
                    factor = std::stod(input);
                } catch (const std::invalid_argument& e) {
                    std::cout << "  Warning: Invalid input. Use no modification factor." << std::endl;
                    factor = 1;
                }
                for (auto npIt = contactIt->second.begin(); npIt != contactIt->second.end(); npIt++) {
                    npIt->factor = factor;
                }
            }
        }
    }

    std::cout << "3. Writing topology file to " << config.topologyFilename << std::endl;

    /** Write topology file **/
    std::map<std::pair<int, int>, std::string> specialAtomTypes;
    std::ofstream topologyFile(config.topologyFilename);
    switch (config.outputType) {
    case GROMACS:
        /** General information **/
        topologyFile << "; Structure based coarse grained model for use with GROMACS\n";
        auto t = std::time(nullptr);
        auto tm = *std::localtime(&t); // 2022-05-31 10:36:52.170316
        topologyFile << "; Date and time: " << std::put_time(&tm, "%F %T") << "\n";
        topologyFile << "; Reference: https://manual.gromacs.org/documentation/current/reference-manual/topologies/topology-file-formats.html\n\n";
        /** Save configuration **/
        if (config.saveConfig) {
            topologyFile << "; Configuration file:\n";
            topologyFile << config.printString() << "\n";
        }
        /** Set defaults **/
        topologyFile << "[ defaults ]\n";
        topologyFile << "; Use combination rule 1 (https://manual.gromacs.org/documentation/current/reference-manual/topologies/parameter-files.html#nbpar)\n";
        topologyFile << "; nbfunc comb-rule gen-pairs fudgeLJ fudgeQQ\n";
        topologyFile << "  1      1         no        1       1\n\n";
        /** Write normal atomtypes **/
        topologyFile << "[ atomtypes ]\n;   name       mass   charge    ptype  c6             c12\n";
        double V = 0.0;
        if (config.nonNativeAttraction) {
            V = config.epsilon * pow(config.LJRadius, 6);
        }
        double W = config.epsilon * pow(config.LJRadius, 12);
        if (config.uniformMass) {
            topologyFile << string_format("%8s %10.4f %10.6f  %-5s  %.6E   %.6E\n", "CA", 1.0, 0.0, "A", V , W);
        } else {
            for (auto atomTypeIt = system.atomTypes.begin(); atomTypeIt != system.atomTypes.end(); atomTypeIt++) {
                std::string aminoAcid = aminoAcidNameTable1to3.at(*atomTypeIt);
                topologyFile << string_format("%8s %10.4f %10.6f  %-5s  %.6E   %.6E\n", aminoAcid.c_str(), massTable.at(*atomTypeIt), 0.0, "A", V , W);
            }
        }
        double factor = pow(2.0, 1.0/6.0);
        unsigned int count = 0;
        if (intermolecularNativeContacts.size() > 0) {
            /** Write special atom types **/
            for (auto chainPairIt = intermolecularNativeContacts.begin(); chainPairIt != intermolecularNativeContacts.end(); chainPairIt++) {
                for (auto contactIt = chainPairIt->second.begin(); contactIt != chainPairIt->second.end(); contactIt++) {
                    std::pair<int, int> keys[2] = {std::make_pair(contactIt->i , chainPairIt->first.first), std::make_pair(contactIt->j, chainPairIt->first.second)};
                    for (unsigned int ki = 0; ki < 2; ki++) {
                        if (specialAtomTypes.find(keys[ki]) == specialAtomTypes.end()) {
                            std::string name = string_format("CA%d", count);
                            // If there are more than 999 atom types the format "CAxxx" does not work for the coordinate file and we use other characters
                            if (count > 999) {
                                name = intToBase26(count - 1000);
                            }
                            specialAtomTypes.insert(std::make_pair(keys[ki], name));
                            std::string aminoAcidName = aminoAcidNameTable1to3.at(system.chains[chainPairIt->first.first].beads()[0][contactIt->i].residueName);
                            double mass = (config.uniformMass) ? 1.0 : massTable.at(system.chains[chainPairIt->first.first].beads()[0][contactIt->i].residueName);
                            topologyFile << string_format("%8s %10.4f %10.6f  %-5s  %.6E   %.6E\n", name.c_str(), mass, 0.0, "A", V , W);
                            count++;
                        }
                    }
                }
            }
            /** Write intermolecular native contacts **/
            topologyFile << "\n[ nonbond_params ]\n";        
            if (config.saveNonbondedInfo) {
                topologyFile << "; r0  : Native contact equilibrium distance\n";
                topologyFile << "; id  : Interaction surface number\n";
                topologyFile << "; f   : Interaction modification factor\n";
                topologyFile << "; type: Possible atomistic interaction type\n";
                topologyFile << ";     ai       aj type                V                W                r0 id  f    type\n";
            } else {
                topologyFile << ";     ai       aj type                V                W\n";
            }
            unsigned int idCount = 1;
            std::map<unsigned int, std::array<unsigned int, 5>> nonbondedInfo = {{1, {0, 0, 0, 0, 0}}};
            for (auto chainPairIt = intermolecularNativeContacts.begin(); chainPairIt != intermolecularNativeContacts.end(); chainPairIt++) {
                int addToCount = 1;
                for (auto contactIt = chainPairIt->second.begin(); contactIt != chainPairIt->second.end(); contactIt++) {
                    double V = config.LJIntermolecular * contactIt->factor * pow(contactIt->distance / factor, 6);
                    double W = config.LJIntermolecular * contactIt->factor * pow(contactIt->distance / factor, 12);                
                    topologyFile << string_format("%8s %8s %4d  %.9E  %.9E",
                        specialAtomTypes.at(std::make_pair(contactIt->i, chainPairIt->first.first)).c_str(),
                        specialAtomTypes.at(std::make_pair(contactIt->j, chainPairIt->first.second)).c_str(), 1, V, W);
                    if (config.saveNonbondedInfo) {
                        Bead bead1 = system.chains[chainPairIt->first.first].beads()[0][contactIt->i];
                        Bead bead2 = system.chains[chainPairIt->first.second].beads()[0][contactIt->j];
                        topologyFile << string_format(" ; %.9E %-3d", contactIt->distance, idCount + contactIt->id);
                        if (contactIt->factor != 1.0) {
                            topologyFile << string_format(" %01.2f", contactIt->factor);
                        } else {
                            topologyFile << "    ";
                        }
                        ContactType contactType = bead1.getContactType(bead2);
                        if (nonbondedInfo.find(idCount + contactIt->id) == nonbondedInfo.end()) {
                            nonbondedInfo.insert(std::pair<unsigned int, std::array<unsigned int, 5>>(idCount + contactIt->id, {0, 0, 0, 0, 0}));
                        }
                        nonbondedInfo[idCount + contactIt->id][contactType - 1] += 1;
                        topologyFile << " " << contactTypeNames.at(contactType);
                    }  
                    addToCount = std::max(std::max(1, contactIt->id), addToCount);
                    topologyFile << "\n";
                }
                idCount += addToCount;
            }
            if (config.saveNonbondedInfo) {
                topologyFile << "; Summary of potential bond types\n";
                topologyFile << ";    id Total Other H-bond  Salt Disulf Hydrophobic\n";
                for (auto infoIt = nonbondedInfo.begin(); infoIt != nonbondedInfo.end(); infoIt++) {
                    unsigned int total = infoIt->second.at(0) + infoIt->second.at(1) + infoIt->second.at(2) + infoIt->second.at(3) + infoIt->second.at(4);
                    topologyFile << string_format("; %5d %5d %5d %6d %5d %6d %5d\n",
                        infoIt->first, total, infoIt->second.at(0), infoIt->second.at(1), infoIt->second.at(2), infoIt->second.at(3), infoIt->second.at(4));
                }
            }
        }
        /** Prepare molecule splitting **/
        std::vector<std::pair<std::string, unsigned int>> molecules;
        std::string topologySplitFiles = "";
        if (config.splitTopologyFiles.length() > 0) {
            size_t ind = config.splitTopologyFiles.find_last_of('.');
            if (ind != std::string::npos) {
                topologySplitFiles = config.splitTopologyFiles.replace(ind, 1, "%d.");
            } else {
                topologySplitFiles = config.splitTopologyFiles + "%d";
            }
        }
        std::ofstream topologySubFile;
        std::ofstream* topologySubFilePtr = &topologyFile;
        unsigned int ignoredDihedralsCount = 0;
        /** Iterate over all chains **/
        for (unsigned int c = 0; c < system.chains.size(); c++) {
            if (chainIndicesMap.find(c) != chainIndicesMap.end()) {
                continue;
            }
            if (topologySplitFiles.length() > 0) {
                std::string topologySubFilename = string_format(topologySplitFiles, c);
                topologySubFile = std::ofstream(topologySubFilename);
                topologySubFilePtr = &topologySubFile;
                topologyFile << string_format("\n#include \"%s\"\n", topologySubFilename.c_str());
            }
            /** Write molecule type **/
            (*topologySubFilePtr) << "\n[ moleculetype ]\n";
            (*topologySubFilePtr) << "; name               nrexcl\n";
            std::string moleculeName = (system.chains.size() > 1) ? string_format("%s%d", config.moleculeName.c_str(), c) : config.moleculeName;
            molecules.push_back(std::make_pair(moleculeName, system.chains[c].count()));
            (*topologySubFilePtr) << string_format("%s        3\n", moleculeName.c_str());
            /** Write beads **/
            (*topologySubFilePtr) << "\n[ atoms ]\n";
            (*topologySubFilePtr) << ";    nr     type   resnr residue atom    cgnr   charge\n";
            for (unsigned int i = 0; i < system.chains[c].size(); i++) {
                std::string beadType = (specialAtomTypes.find(std::make_pair(i, c)) != specialAtomTypes.end()) ?
                    specialAtomTypes.at(std::make_pair(i, c)) :
                    ((config.uniformMass) ? "CA" : aminoAcidNameTable1to3.at(system.chains[c].beads()[0][i].residueName));
                (*topologySubFilePtr) << string_format("%7d %8s %7d %7s %4s %7d      0.0\n",
                    i+1, beadType.c_str(), i+1, aminoAcidNameTable1to3.at(system.chains[c].beads()[0][i].residueName).c_str(), "CA", i+1);
            }
            /** Write bonds **/
            (*topologySubFilePtr) << "\n[ bonds ]\n";
            (*topologySubFilePtr) << ";    ai      aj func           r0(nm)              Kb\n";
            for (unsigned int i = 0; i < system.chains[c].size() - 1; i++) {
                double distance = getDistance(system.chains[c].beads()[0][i].coordinates, system.chains[c].beads()[0][i+1].coordinates);
                (*topologySubFilePtr) << string_format("%7d %7d %4d  %.9E %.9E\n", i+1, i+2, 1, distance, config.KBondLength * config.epsilon);
            }
            /** Write angles **/
            (*topologySubFilePtr) << "\n[ angles ]\n";
            (*topologySubFilePtr) << ";    ai      aj      ak func         th0(deg)              Ka\n";
            count = 0;
            for (unsigned int i = 0; i < system.chains[c].size() - 2; i++) {
                std::vector<double> angles = {getAngle(system.chains[c].beads()[0][i].coordinates, system.chains[c].beads()[0][i+1].coordinates, system.chains[c].beads()[0][i+2].coordinates)};
                for (unsigned int modelId = 1; modelId < system.chains[c].numberOfModels(); modelId++) {
                    double angle = getAngle(system.chains[c].beads()[modelId][i].coordinates, system.chains[c].beads()[modelId][i+1].coordinates, system.chains[c].beads()[modelId][i+2].coordinates);
                    for (auto angleIt = angles.begin(); angleIt != angles.end(); angleIt++) {
                        if (abs(angle - *angleIt) > config.minAngleDifference) {
                            angles.push_back(angle);
                            break;
                        }
                    }
                }
                if (angles.size() == 1) {
                    (*topologySubFilePtr) << string_format("%7d %7d %7d %4d  %.9E %.9E\n", i+1, i+2, i+3, 1, angles[0], config.KBondAngle * config.epsilon);
                } else {
                    tabulateAngles(angles, config.tabelDirectory, string_format("table_a%d.xvg", count), config.tabulatedPoints);
                    (*topologySubFilePtr) << string_format("%7d %7d %7d %4d  %15d %.9E ; %s\n", i+1, i+2, i+3, 8, count, config.KBondAngle * config.epsilon, vectorToString(angles, "%.9E").c_str());
                    count++;
                }
            }
            /** Write dihedrals **/
            (*topologySubFilePtr) << "\n[ dihedrals ]\n";
            (*topologySubFilePtr) << ";    ai      aj      ak      al func        phi0(deg)               Kd  mult\n";
            count = 0;
            for (unsigned int i = 0; i < system.chains[c].size() - 3; i++) {
                if (config.angleDihedralCutOff > 0) {
                    double angle1 = getAngle(system.chains[c].beads()[0][i].coordinates, system.chains[c].beads()[0][i+1].coordinates, system.chains[c].beads()[0][i+2].coordinates);
                    double angle2 = getAngle(system.chains[c].beads()[0][i+1].coordinates, system.chains[c].beads()[0][i+2].coordinates, system.chains[c].beads()[0][i+3].coordinates);
                    if (angle1 > config.angleDihedralCutOff || angle2 > config.angleDihedralCutOff) {
                        ignoredDihedralsCount++;
                        continue;
                    }
                }
                std::vector<double> dihedrals = {180 + getDihedrals(system.chains[c].beads()[0][i].coordinates, system.chains[c].beads()[0][i+1].coordinates, system.chains[c].beads()[0][i+2].coordinates, system.chains[c].beads()[0][i+3].coordinates)};
                for (unsigned int modelId = 1; modelId < system.chains[c].numberOfModels(); modelId++) {
                    double dihedral = 180 + getDihedrals(system.chains[c].beads()[modelId][i].coordinates, system.chains[c].beads()[modelId][i+1].coordinates, system.chains[c].beads()[modelId][i+2].coordinates, system.chains[c].beads()[modelId][i+3].coordinates);
                    for (auto dihedralIt = dihedrals.begin(); dihedralIt != dihedrals.end(); dihedralIt++) {
                        if (abs(dihedral - *dihedralIt) > config.minDihedralDifference && abs(360 - dihedral + *dihedralIt) > config.minDihedralDifference) {
                            dihedrals.push_back(dihedral);
                            break;
                        }
                    }
                }
                if (dihedrals.size() == 1) {
                    (*topologySubFilePtr) << string_format("%7d %7d %7d %7d %4d  %.9E  %.9E  %d\n", i+1, i+2, i+3, i+4, 1, dihedrals[0], config.KDihedral1 * config.epsilon, 1);
                    (*topologySubFilePtr) << string_format("%7d %7d %7d %7d %4d  %.9E  %.9E  %d\n", i+1, i+2, i+3, i+4, 1, dihedrals[0] * 3, config.KDihedral2 * config.epsilon, 3);
                } else {
                    tabulateDihedrals(dihedrals, config.tabelDirectory, string_format("table_d%d.xvg", count), 1, config.tabulatedPoints);
                    (*topologySubFilePtr) << string_format("%7d %7d %7d %7d %4d  %15d  %.9E ; %s\n", i+1, i+2, i+3, i+4, 8, count, config.KDihedral1 * config.epsilon, vectorToString(dihedrals, "%.9E").c_str());
                    count++;
                }
            }
            /** Write pairs and exclusions for intramolecular native contacts **/
            (*topologySubFilePtr) << "\n[ pairs ]\n";
            (*topologySubFilePtr) << ";    ai      aj type                V                W                r0\n";
            std::stringstream exclusionStream;
            std::set<NativePair> intramolecularContacts = system.chains[c].getIntramolecularContacts(config);
            for (auto nativeContactIt = intramolecularContacts.begin(); nativeContactIt != intramolecularContacts.end(); nativeContactIt++) {
                double V = config.epsilon * pow(nativeContactIt->distance / factor, 6);
                double W = config.epsilon * pow(nativeContactIt->distance / factor, 12);
                (*topologySubFilePtr) << string_format("%7d %7d %4d  %.9E  %.9E ; %.9E\n", nativeContactIt->i + 1, nativeContactIt->j + 1, 1, V, W, nativeContactIt->distance);
                exclusionStream << string_format("%7d %7d\n", nativeContactIt->i + 1, nativeContactIt->j + 1);    
            }
            (*topologySubFilePtr) << "\n[ exclusions ]\n";
            (*topologySubFilePtr) << ";    ai      aj\n";
            (*topologySubFilePtr) << exclusionStream.str();
        }
        if (ignoredDihedralsCount == 1) {
            std::cout << "Warning: 1 dihedral was ignored because adjoining bond angles were close to 180°." << std::endl;
        } else if (ignoredDihedralsCount > 1) {
            std::cout << "Warning: " << ignoredDihedralsCount << " dihedrals were ignored because adjoining bond angles were close to 180°." << std::endl;
        }
        /** Write general system information **/
        topologyFile << "\n[ system ]\n";
        topologyFile << "; name\n";
        topologyFile << string_format("%s\n", config.moleculeName.c_str());
        topologyFile << "\n[ molecules ]\n";
        topologyFile << "; name               molecule-number\n";
        for (auto moleculeIt = molecules.begin(); moleculeIt != molecules.end(); moleculeIt++) {
            topologyFile << string_format("%s        %d\n", moleculeIt->first.c_str(), moleculeIt->second);
        }
        break;
    }
    topologyFile.close();

    std::cout << "4. Writing coordinate file to " << config.coordinateFilename << std::endl;

    /** Write coordinate file **/
    std::ofstream coordinateFile(config.coordinateFilename);
    switch (config.outputType) {
        case GROMACS:
            coordinateFile << string_format("Coordinate file for a structure based model\n%lu\n", systemSize);
            for (unsigned int c = 0; c < originalChains.size(); c++) {
                for (unsigned int i = 0; i < originalChains[c].size(); i++) {
                    Bead bead = originalChains[c].beads()[0][i];
                    int index = (i+1) % 100000;
                    /** Determine atomtype **/
                    auto chainIndicesIt = chainIndicesMap.find(c);
                    unsigned int newChainIndex = (chainIndicesIt == chainIndicesMap.end()) ? c : chainIndicesIt->second;
                    std::pair<int, int> atomTypesKey = std::make_pair(i, newChainIndex);
                    auto specialAtomTypesIt = specialAtomTypes.find(atomTypesKey);
                    std::string atomType = (specialAtomTypesIt == specialAtomTypes.end()) ? "CA" : specialAtomTypesIt->second;
                    /** Write to file **/
                    coordinateFile << string_format("%5i%-5s%5s%5i%8.3f%8.3f%8.3f\n", 
                        index, aminoAcidNameTable1to3.at(bead.residueName).c_str(), "CA",
                        index, bead.coordinates.x, bead.coordinates.y, bead.coordinates.z);
                }
            }
            coordinateFile << string_format("%10.5f %10.5f %10.5f\n", boxSize.x, boxSize.y, boxSize.z);
            break;
    }
    coordinateFile.close();
}