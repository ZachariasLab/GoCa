#include "system.hpp"
#include "helper.hpp"
#include "chain.hpp"

System::System(Config& config) {
    std::ifstream infile(config.inputPdb);

    /** Check if input structure file exists **/
    if (!infile.good()) {
        std::cerr << "Error: The file \"" << config.inputPdb << "\" does not exist" << std::endl;
        exit(EXIT_FAILURE);
    }

    /** Determine input file type **/
    FileType fileType;
    if (config.inputPdb.find(".pdb") != std::string::npos) {
        fileType = PDB;
    } else if (config.inputPdb.find(".cif") != std::string::npos) {
        fileType = CIF;
    } else {
        std::cerr << "Error: Filetype of \"" << config.inputPdb << "\" not supported." << std::endl;
        exit(EXIT_FAILURE);
    }
    
    /** Prepare variables for system generation **/
    std::cout << "Info: Reading protein file: " << config.inputPdb << std::endl;
    unsigned int modelId = 0;
    unsigned int chainId = 0;
    unsigned int residueId = 1;
    std::vector<std::string> pdbChainNames;
    int lastResiduePdbId = -1;
    std::string lastChainPdbName;
    Bead bead;
    bool firstAtom = true;
    bool printedOccupancyWarning = false;
    unsigned int atomCount = 0;
    unsigned int beadCount = 0;
    std::vector<unsigned int> modelIndices;
    bool foundHatom = false;

    /** Iterate over lines in file **/
    for(std::string line; getline(infile, line);) {
        if (line.substr(0, 4) == "ATOM") {
            if (modelIndices.size() == 0 && modelId == 0) {
                modelIndices.push_back(modelId);
            }
            if (std::find(config.modelIds.begin(), config.modelIds.end(), modelId) == config.modelIds.end() && modelId != 0) {
                /** Skip atom because it is from an excluded model **/
                continue;
            }
            AtomEntry atomEntry(line, fileType);
            if (atomEntry.altAtom) {
                /** If there are multiple coordinate sets for the same atom only the first one is used and others are ignored **/
                if (!printedOccupancyWarning) {
                    std::cout << "Info: Found occupancy values unequal to one. Only first atoms are used." << std::endl;
                    printedOccupancyWarning = true;
                }
                continue;
            }
            if (firstAtom) {
                /** If it is the first atom we set some variable which help to recognize residue and chain changes **/
                lastResiduePdbId = atomEntry.residueId;
                lastChainPdbName = atomEntry.chainName;
                pdbChainNames.push_back(atomEntry.chainName);
                chains.push_back(Chain());
                firstAtom = false;
            }
            atomCount++;
            if (atomEntry.residueId != lastResiduePdbId) {
                /** If the atom has a new residueId the previous bead is stored in the chains list **/
                if (bead.modelId == -1) {
                    std::cerr << "Error: No Cα atom in residue " << lastResiduePdbId << "." << std::endl;
                    exit(EXIT_FAILURE);
                }
                chains[chainId].add(bead.modelId, bead);
                bead = Bead();
                beadCount++;
                residueId++;
                lastResiduePdbId = atomEntry.residueId;
            }
            if (atomEntry.chainName != lastChainPdbName) {
                /** If the atom has a new chain name we start a new entry in the chains list **/
                int index = findStringInVector(pdbChainNames, atomEntry.chainName);
                if (index < 0) {
                    pdbChainNames.push_back(atomEntry.chainName);
                    chains.push_back(Chain());
                    chainId = pdbChainNames.size() - 1;
                } else {
                    chainId = index;
                }
                residueId = 1;
                lastChainPdbName = atomEntry.chainName;
            }
            if (atomEntry.atomName == "CA") {
                /** If the atom is the C-Alpha atom of the residue we add all important information to the new Bead **/
                bead.coordinates = Coordinates(atomEntry.x, atomEntry.y, atomEntry.z);
                bead.atoms.emplace(bead.atoms.begin(), atomEntry.x, atomEntry.y, atomEntry.z);
                bead.atomsTypes.insert(bead.atomsTypes.begin(), atomEntry.atomName);
                auto aminoAcidIt = aminoAcidNameTable3to1.find(atomEntry.residueName);
                if (aminoAcidIt == aminoAcidNameTable3to1.end()) {
                    std::cerr << "Error: Unknown amino acid \"" << atomEntry.residueName << "\" in input structure." << std::endl;
                    exit(EXIT_FAILURE);
                }
                bead.residueName = aminoAcidIt->second;
                if (beadVdwRadii.find(bead.residueName) == beadVdwRadii.end()) {
                    std::cerr << "Error: Unknown residue type: " << bead.residueName << std::endl;
                    exit(EXIT_FAILURE);
                }
                bead.mass = (config.uniformMass) ? 1.0 : massTable.at(bead.residueName);
                atomTypes.insert(bead.residueName);
                bead.modelId = modelIndices.size() - 1;
                bead.chainId = chainId;
                bead.residueId = residueId;
                /** Adjust box size **/
                if (atomEntry.x < minBoxCorner.x) minBoxCorner.x = atomEntry.x;
                if (atomEntry.x > maxBoxCorner.x) maxBoxCorner.x = atomEntry.x;
                if (atomEntry.y < minBoxCorner.y) minBoxCorner.y = atomEntry.y;
                if (atomEntry.y > maxBoxCorner.y) maxBoxCorner.y = atomEntry.y;
                if (atomEntry.z < minBoxCorner.z) minBoxCorner.z = atomEntry.z;
                if (atomEntry.z > maxBoxCorner.z) maxBoxCorner.z = atomEntry.z;
            } else {
                /** Add atom to the list of non-C-Alpha atom coordinates **/
                bead.atoms.emplace_back(atomEntry.x, atomEntry.y, atomEntry.z);
                bead.atomsTypes.push_back(atomEntry.atomName);
                if (atomicVdwRadii.find(atomEntry.atomName.at(0)) == atomicVdwRadii.end()) {
                    std::cerr << "Error: Unknown atom element: " << atomEntry.atomName.at(0) << std::endl;
                    exit(EXIT_FAILURE);
                }
                if (!foundHatom && atomEntry.atomName.at(0) == 'H') {
                    foundHatom = true;
                }
            }
        } else if (line.substr(0, 5) == "MODEL") {
            modelId++;
            if (std::find(config.modelIds.begin(), config.modelIds.end(), modelId) != config.modelIds.end()) {
                modelIndices.push_back(modelId);
                residueId = 0;
            }
        }
    }
    /** Check model validity **/
    for (auto modelIdIt = config.modelIds.begin(); modelIdIt != config.modelIds.end(); modelIdIt++) {
        if (std::find(modelIndices.begin(), modelIndices.end(), *modelIdIt) == modelIndices.end() && !(modelIndices.size() == 1 && modelIndices[0] == 0)) {
            std::cerr << "Error: Model id " << *modelIdIt << " does not exist in structure file" << std::endl;
            exit(EXIT_FAILURE);
        }
    }
    /** Finish list of beads **/
    chains[chainId].add(bead.modelId, bead);
    beadCount++;
    modelCount = modelIndices.size();

    /** Check for unexpected chain size **/
    for (auto it = chains.begin(); it != chains.end(); it++) {
        if (it->size() < 5) {
            std::cout << "Warning: Found small chain with less than five residues. If this is not intended please check your structure file." << std::endl;
            break;
        }
    }

    /** Print warning if structure does not contain H atoms **/
    if (!foundHatom && (!config.isConfigured("atomic-vdw-distance") || !config.isConfigured("atomic-cutoff"))) {
        std::cout << "Warning: Structure does not contain H atoms. It is recommended to adjust default parameters for native contact cutoff." << std::endl;
    }    

    /** Print model info **/
    if (modelIndices.size() == 1) {
        std::cout << "Info: Got a total of " << beadCount << " Cα-atoms from " << atomCount << " atoms in " << chains.size() << " chain(s)" << std::endl;
    } else {
        std::cout << "Info: Got a total of " << beadCount << " Cα-atoms in " << chains.size() << " chain(s) and " << modelIndices.size() << " model(s)" << std::endl;
    }
    /** Check if all chains are valid **/
    for (auto chainIt = chains.begin(); chainIt != chains.end(); chainIt++) {
        if (!chainIt->valid()) {
            std::cerr << "Error: All models per chain must have the same amino acid sequence" << std::endl;
            exit(EXIT_FAILURE);
        }
    }
};

void System::translateBeads(const Coordinates& vec) {
    for (auto chainIt = chains.begin(); chainIt != chains.end(); chainIt++) {
        chainIt->translateBeads(vec);
    }
    minBoxCorner.translate(vec);
    maxBoxCorner.translate(vec);
};

const Coordinates System::getBoxSize(const double& boxPadding) {
    return Coordinates(
        maxBoxCorner.x - minBoxCorner.x + 2 * boxPadding,
        maxBoxCorner.y - minBoxCorner.y + 2 * boxPadding,
        maxBoxCorner.z - minBoxCorner.z + 2 * boxPadding);
};

const Coordinates System::getBoxTranslationVector(const double& boxPadding) {
    return Coordinates(
        boxPadding - minBoxCorner.x,
        boxPadding - minBoxCorner.y,
        boxPadding - minBoxCorner.z);
};

const int System::size() {
    int sum = 0;
    for (auto chainIt = chains.begin(); chainIt != chains.end(); chainIt++) {
        sum += chainIt->size();
    }
    return sum;
};

const std::map<std::string, std::vector<std::vector<std::pair<unsigned int, unsigned int>>>> System::getChainClusters(const double& clusterCutOff) {
    std::map<std::string, std::vector<std::vector<std::pair<unsigned int, unsigned int>>>> clusters;
    for (unsigned int c = 0; c < chains.size(); c++) {
        std::string chainString = chains[c].getChainString();
        if (chainString.size() <= 1) {
            std::cerr << "Error: Chain without residues. Please check your structure file." << std::endl;
            exit(EXIT_FAILURE);
        }
        if (clusters.find(chainString) == clusters.end()) {
            for (unsigned int m = 0; m < chains[c].numberOfModels(); m++) {
                clusters.insert(std::make_pair(chainString, std::vector<std::vector<std::pair<unsigned int, unsigned int>>>(1, std::vector<std::pair<unsigned int, unsigned int>>(1, std::make_pair(c, m)))));
            }
        } else {
            for (unsigned int m = 0; m < chains[c].numberOfModels(); m++){
                std::vector<double> rmsdSum(clusters[chainString].size(), 0.0);
                unsigned int tooBigDifferenceCount = 0;
                for (unsigned int k = 0; k < clusters[chainString].size(); k++) {
                    for (auto jIt = clusters[chainString][k].begin(); jIt != clusters[chainString][k].end(); jIt++) {
                        double rmsd = chains[c].rmsd(chains[jIt->first], m, jIt->second).rmsd;
                        if (rmsd < clusterCutOff) {
                            rmsdSum[k] += rmsd;
                        } else {
                            rmsdSum[k] = -1.0;
                            tooBigDifferenceCount++;
                            break;
                        }
                    }
                }
                if (tooBigDifferenceCount == clusters[chainString].size()) {
                    clusters[chainString].push_back(std::vector<std::pair<unsigned int, unsigned int>>(1, std::make_pair(c, m)));
                } else {
                    int minIndex = 0;
                    double minAverageRmsd = std::numeric_limits<double>::max();
                    for (unsigned int k = 0; k < clusters[chainString].size(); k++) {
                        if (rmsdSum[k] >= 0) {
                            double averageRmsd = rmsdSum[k] / clusters[chainString][k].size();
                            if (averageRmsd < minAverageRmsd) {
                                minAverageRmsd = averageRmsd;
                                minIndex = k;
                            }
                        }
                    }
                    clusters[chainString][minIndex].push_back(std::make_pair(c, m));
                }
            }
        }
    }
    return clusters;
};

std::map<unsigned int, unsigned int> System::mergeChainsAndIntermolecularContacts(
        const std::map<std::string, std::vector<std::vector<std::pair<unsigned int, unsigned int>>>>& clusters,
        std::map<std::pair<int, int>, std::set<NativePair>>& intermolecularNativeContacts,
        bool varyingInteractionFactors
    ) {
    bool wroteInteractionTitle = false;
    std::map<unsigned int, unsigned int> chainIndicesMap;
    for (auto clIt = clusters.begin(); clIt != clusters.end(); clIt++) {
        /** Merge chains **/
        int chainId = clIt->second.begin()->begin()->first;
        std::set<int> mergedChainIds = {chainId};
        for (auto mIt = clIt->second.begin(); mIt != clIt->second.end(); mIt++) {
            for (auto cIt = std::next(mIt->begin()); cIt != mIt->end(); cIt++) {
                chains[mIt->begin()->first].mergeDirectAdd(chains[cIt->first], mIt->begin()->second, cIt->second);
                mergedChainIds.insert(cIt->first);
                chainIndicesMap.insert(std::make_pair(cIt->first, chainId));
            }
            chains[mIt->begin()->first].mergeDirectDivide(mIt->begin()->second, mIt->size());
            if (mIt != clIt->second.begin()) {
                mergedChainIds.insert(mIt->begin()->first);
                chains[chainId].mergeAsModel(chains[mIt->begin()->first]);
                chainIndicesMap.insert(std::make_pair(mIt->begin()->first, chainId));
            }
        }
        /** Merge native contacts **/
        std::map<std::tuple<int, int, int, int>, std::tuple<int, double, double>> intermolecularContactBuilder;
        std::vector<std::set<std::pair<int, int>>> intermolecularContactAreasBuilder;
        for (auto contactListIt = intermolecularNativeContacts.begin(); contactListIt != intermolecularNativeContacts.end(); contactListIt++) {
            bool foundFirst = std::find(mergedChainIds.begin(), mergedChainIds.end(), contactListIt->first.first) != mergedChainIds.end();
            bool foundSecond = std::find(mergedChainIds.begin(), mergedChainIds.end(), contactListIt->first.second) != mergedChainIds.end();
            /** Get contact areas **/
            if (varyingInteractionFactors && mergedChainIds.size() > 1 && foundFirst && foundSecond) {
                if (intermolecularContactAreasBuilder.size() == 0) {
                    intermolecularContactAreasBuilder.resize(1);
                    for (auto contactIt = contactListIt->second.begin(); contactIt != contactListIt->second.end(); contactIt++) {
                        intermolecularContactAreasBuilder[0].insert(makeOrderedPair(contactIt->i, contactIt->j));
                    }
                } else {
                    bool found = false;
                    for (auto groupIt = intermolecularContactAreasBuilder.begin(); groupIt != intermolecularContactAreasBuilder.end(); groupIt++){
                        for (auto contactIt = contactListIt->second.begin(); contactIt != contactListIt->second.end(); contactIt++) {
                            if (groupIt->find(makeOrderedPair(contactIt->i, contactIt->j)) != groupIt->end()) {
                                found = true;
                                break;
                            }
                        }
                        if (found) {
                            for (auto contactIt = contactListIt->second.begin(); contactIt != contactListIt->second.end(); contactIt++) {
                                groupIt->insert(makeOrderedPair(contactIt->i, contactIt->j));
                            }
                            break;
                        }
                    }
                    if (!found) {
                        int size = intermolecularContactAreasBuilder.size();
                        intermolecularContactAreasBuilder.resize(size + 1);
                        for (auto contactIt = contactListIt->second.begin(); contactIt != contactListIt->second.end(); contactIt++) {
                            intermolecularContactAreasBuilder[size].insert(makeOrderedPair(contactIt->i, contactIt->j));
                        }
                    }
                }
            }
            /** Merge contacts **/
            for (auto contactIt = contactListIt->second.begin(); contactIt != contactListIt->second.end(); contactIt++) {
                std::tuple<int, int, int, int> key;
                if (foundFirst && foundSecond) {
                    if (contactIt->i < contactIt->j) {
                        key = std::make_tuple(chainId, chainId, contactIt->i, contactIt->j);
                    } else {
                        key = std::make_tuple(chainId, chainId, contactIt->j, contactIt->i);
                    }
                } else if (foundFirst) {
                    if (chainId < contactListIt->first.second) {
                        key = std::make_tuple(chainId, contactListIt->first.second, contactIt->i, contactIt->j);
                    } else {
                        key = std::make_tuple(contactListIt->first.second, chainId, contactIt->i, contactIt->j);
                    }
                } else if (foundSecond) {
                    if (contactListIt->first.first < chainId) {
                        key = std::make_tuple(contactListIt->first.first, chainId, contactIt->i, contactIt->j);
                    } else {
                        key = std::make_tuple(chainId, contactListIt->first.first, contactIt->i, contactIt->j);
                    }
                } else {
                    key = std::make_tuple(contactListIt->first.first, contactListIt->first.second, contactIt->i, contactIt->j);
                }
                auto findIt = intermolecularContactBuilder.find(key);
                if (findIt == intermolecularContactBuilder.end()) {
                    intermolecularContactBuilder.insert(std::make_pair(key, std::make_tuple(1, contactIt->distance, contactIt->factor)));
                } else {
                    std::tuple<int, double, double> dvalue = findIt->second;
                    findIt->second = std::make_tuple(std::get<0>(dvalue) + 1, std::get<1>(dvalue) + contactIt->distance, std::get<2>(dvalue) + contactIt->factor);
                }
            }
        }
        /** Process areas and ask user for modification factors **/
        std::map<std::pair<int, int>, std::tuple<double, int>> contactAreaInfo;        
        if (varyingInteractionFactors && mergedChainIds.size() > 1) {
            if (!wroteInteractionTitle) {
                wroteInteractionTitle = true;
                std::cout << "Intermolecular interaction modification factors between merged chains:" << std::endl;
            }
            if (intermolecularContactAreasBuilder.size() == 1) {
                std::cout << "  Found 1 contact area";
                std::cout << " for the chain with " << clIt->first.size() << " residues with the following number of contacts: ";
            } else {
                std::cout << "  Found " << intermolecularContactAreasBuilder.size() << " contact areas";
                std::cout << " for the chain with " << clIt->first.size() << " residues each with the following number of contacts: ";
            }
            for (auto areaIt = intermolecularContactAreasBuilder.begin(); areaIt != intermolecularContactAreasBuilder.end(); areaIt++) {
                std::cout << areaIt->size() << " ";
            }
            std::cout << std::endl << "  Please enter space seperated modification factors for each interaction group: ";
            std::string input;
            std::getline(std::cin, input);
            std::vector<double> factors;
            try {
                factors = splitToDouble(input, ' ');
            } catch (const std::invalid_argument& e) {
                // Ignore error. We use default values and print an warning message (see below)
            }
            if (factors.size() == intermolecularContactAreasBuilder.size()) {
                for (auto areaBuilderId = 0; areaBuilderId < intermolecularContactAreasBuilder.size(); areaBuilderId++) {
                    for (auto contactIt = intermolecularContactAreasBuilder[areaBuilderId].begin(); contactIt != intermolecularContactAreasBuilder[areaBuilderId].end(); contactIt++) {
                        contactAreaInfo.insert(std::make_pair(std::make_pair(contactIt->first, contactIt->second), std::make_tuple(factors[areaBuilderId], areaBuilderId)));
                    }
                }
            } else {
                std::cout << "  Warning: Invalid input. No modification factors are used." << std::endl;
            }
        }
        /** Rewrite intermolecular contacts **/
        intermolecularNativeContacts.clear();
        for (auto contactIt = intermolecularContactBuilder.begin(); contactIt != intermolecularContactBuilder.end(); contactIt++) {
            std::pair<int, int> key(std::get<0>(contactIt->first), std::get<1>(contactIt->first));
            double factor = std::get<2>(contactIt->second) / std::get<0>(contactIt->second);
            unsigned int id = 0;
            if (key.first == key.second) {
                auto factorIt = contactAreaInfo.find(makeOrderedPair(std::get<2>(contactIt->first), std::get<3>(contactIt->first)));
                if (factorIt != contactAreaInfo.end()) {
                    factor = std::get<0>(factorIt->second);
                    id = std::get<1>(factorIt->second);
                }
            }
            NativePair contact(std::get<2>(contactIt->first), std::get<3>(contactIt->first),  std::get<1>(contactIt->second) / std::get<0>(contactIt->second), id, factor);
            if (intermolecularNativeContacts.find(key) == intermolecularNativeContacts.end()) {
                intermolecularNativeContacts.insert(std::make_pair(key, std::set<NativePair>{contact}));
            } else {
                intermolecularNativeContacts[key].insert(contact);
            }
        }
    }
    return chainIndicesMap;
};

AtomEntry::AtomEntry(const std::string& line, const FileType fileType) {
    switch (fileType) {
        case PDB: {
            chainName = line.substr(21, 1);
            residueName = line.substr(17, 3);
            atomName = trim(line.substr(12, 4));
            residueId = std::stoi(line.substr(22, 4));
            modelId = -1;
            x = std::stod(line.substr(30, 8)) / 10.0; // Convert from Å to nm
            y = std::stod(line.substr(38, 8)) / 10.0;
            z = std::stod(line.substr(46, 8)) / 10.0;
            std::string altId = line.substr(16, 1);
            altAtom = !(altId == " " || altId == "A");
            break;
        } case CIF: {
            std::vector<std::string> entries = splitToString(line);
            chainName = entries[6];
            residueName = entries[5];
            atomName = entries[3];
            residueId = std::stoi(entries[8]);
            modelId = std::stoi(entries[20]);
            x = std::stod(entries[10]) / 10.0; // Convert from Å to nm
            y = std::stod(entries[11]) / 10.0;
            z = std::stod(entries[12]) / 10.0;
            altAtom = !(entries[4] == "." || entries[4] == "A");
            break;
        }
    }
};