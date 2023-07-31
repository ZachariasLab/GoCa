#include "bead.hpp"
#include "helper.hpp"

const bool Bead::isNativeContact(const Bead& other, const Config& config) const {
    double caDistance = 0;
    if (config.includeBeadCutOff) {
        double CaSqDistance = pow2(coordinates.x - other.coordinates.x) + pow2(coordinates.y - other.coordinates.y) + pow2(coordinates.z - other.coordinates.z);
        if (config.beadCutOff > 0 && CaSqDistance > pow2(config.beadCutOff)) {
            return false;
        }
        if (CaSqDistance > pow2(beadVdwRadii.at(residueName) + beadVdwRadii.at(other.residueName) + config.beadVdwDistance)) {
            return false;
        }
    }
    for (unsigned int i1 = 0; i1 < atoms.size(); i1++) {
        for (unsigned int i2 = 0; i2 < other.atoms.size(); i2++) {
            if ((this != &other || &atoms[i1] != &other.atoms[i2]) &&
                (!config.excludeHAtomsForCutoff || (atomsTypes[i1].at(0) != 'H' && other.atomsTypes[i2].at(0) != 'H'))) {
                double sqDistance = pow2(atoms[i1].x - other.atoms[i2].x) + pow2(atoms[i1].y - other.atoms[i2].y) + pow2(atoms[i1].z - other.atoms[i2].z);
                if (config.atomicCutOff > 0) {
                    if (sqDistance < pow2(config.atomicCutOff)) {
                        return true;
                    }
                } else {
                    double vdwCutOffSquared = pow2(atomicVdwRadii.at(atomsTypes[i1].at(0)) + atomicVdwRadii.at(other.atomsTypes[i2].at(0)) + config.atomicVdwDistance);
                    if (sqDistance < vdwCutOffSquared) {
                        return true;
                    }
                }
            }
        }
    }
    return false;
};

const bool Bead::isAnion(const unsigned int& index) const {
    return (residueName == 'D' /*ASP*/ && (atomsTypes[index] == "OD1" || atomsTypes[index] == "OD2")) ||
        (residueName == 'E' /*GLU*/ && (atomsTypes[index] == "OE1" || atomsTypes[index] == "OE2"));
};
const bool Bead::isCation(const unsigned int& index) const {
    return (residueName == 'K' /*LYS*/ && atomsTypes[index] == "NZ") ||
        (residueName == 'R' /*ARG*/ && (atomsTypes[index] == "NH1" || atomsTypes[index] == "NH2"));
};

const unsigned int Bead::getClostestResidueAtom(const unsigned int& index) const {
    double closestSqDistance = std::numeric_limits<double>::max();
    unsigned int closestIndex = 0;
    for (unsigned int l = 0; l < atoms.size(); l++) {
        if (l != index) {
            double sqDistance = pow2(atoms[index].x - atoms[l].x) + pow2(atoms[index].y - atoms[l].y) + pow2(atoms[index].z - atoms[l].z);
            if (sqDistance < closestSqDistance) {
                closestSqDistance = sqDistance;
                closestIndex = l;
            }
        }
    }
    return closestIndex;
};

const bool Bead::hydrophobicResidue() const {
    // ALA CYS PHE GLY ILE LEU MET PRO VAL TRP
    return residueName == 'A' || residueName == 'C' || residueName == 'F' || residueName == 'G' || residueName == 'I' ||
        residueName == 'L' || residueName == 'M' || residueName == 'P' || residueName == 'V' || residueName == 'W';
}

const ContactType Bead::getContactType(const Bead& other) const {
    /** Check if disulfide bridge **/
    if (residueName == 'C' && other.residueName == 'C') {
        return DISULF;
    }
    /** Get indices and distance of closest and second closest atom **/
    for (unsigned int i1 = 0; i1 < atoms.size(); i1++) {
        for (unsigned int i2 = 0; i2 < other.atoms.size(); i2++) {
            double distanceSquared = pow2(atoms[i1].x - other.atoms[i2].x) + pow2(atoms[i1].y - other.atoms[i2].y) + pow2(atoms[i1].z - other.atoms[i2].z);
            /** Check if salt bridge **/
            if (distanceSquared < 0.4*0.4 && ((isAnion(i1) && other.isCation(i2)) || (isCation(i1) && other.isAnion(i2)))) {
                return SALTB;
            }
            /** Check if hydrogen bond **/
            if (atomsTypes[i1].at(0) == 'H' && other.atomsTypes[i2].at(0) != 'H' && other.atomsTypes[i2].at(0) != 'C' && other.atomsTypes[i2].at(0) != 'S') {
                unsigned int donorIndex = getClostestResidueAtom(i1);
                double acceptorDonorSqDistance = pow2(atoms[donorIndex].x - other.atoms[i2].x) + pow2(atoms[donorIndex].y - other.atoms[i2].y) + pow2(atoms[donorIndex].z - other.atoms[i2].z);                
                if (atomsTypes[donorIndex].at(0) != 'H' && atomsTypes[donorIndex].at(0) != 'C' && atomsTypes[donorIndex].at(0) != 'S' && acceptorDonorSqDistance < 0.35*0.35) {
                    double angle = 180 - getAngle(other.atoms[i2], atoms[i1], atoms[donorIndex]);
                    if (angle < 70) {
                        return HBOND;
                    }
                }
            }
            if (other.atomsTypes[i2].at(0) == 'H' && atomsTypes[i1].at(0) != 'H' && atomsTypes[i1].at(0) != 'C' && atomsTypes[i1].at(0) != 'S') {
                unsigned int donorIndex = other.getClostestResidueAtom(i2);
                double acceptorDonorSqDistance = pow2(other.atoms[donorIndex].x - atoms[i1].x) + pow2(other.atoms[donorIndex].y - atoms[i1].y) + pow2(other.atoms[donorIndex].z - atoms[i1].z);  
                if (other.atomsTypes[donorIndex].at(0) != 'H' && other.atomsTypes[donorIndex].at(0) != 'C' && other.atomsTypes[donorIndex].at(0) != 'S' && acceptorDonorSqDistance < 0.35*0.35) {
                    double angle = 180 - getAngle(atoms[i1], other.atoms[i2], other.atoms[donorIndex]);
                    if (angle < 70) {
                        return HBOND;
                    }
                }
            }
            /** Check if hydrophobic contact **/
            if ((atomsTypes[i1].at(0) == 'C' || atomsTypes[i1].at(0) == 'S') &&
                (other.atomsTypes[i2].at(0) == 'C' || other.atomsTypes[i2].at(0) == 'S') &&
                hydrophobicResidue() && other.hydrophobicResidue() && 
                distanceSquared < pow2(atomicVdwRadii.at(atomsTypes[i1].at(0)) + atomicVdwRadii.at(other.atomsTypes[i2].at(0)) + 0.05)) {
                return HYDROPHOBIC;
            }
        }
    }    
    return OTHER;
};
