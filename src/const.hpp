#pragma once

#include <map>
#include <string>

static std::map<std::string, char> const aminoAcidNameTable3to1 = { 
    {"ALA", 'A'}, {"ARG", 'R'}, {"ASN", 'N'}, {"ASP", 'D'}, {"CYS", 'C'},
    {"GLU", 'E'}, {"GLN", 'Q'}, {"GLY", 'G'}, {"HIS", 'H'}, {"ILE", 'I'},
    {"LEU", 'L'}, {"LYS", 'K'}, {"MET", 'M'}, {"PHE", 'F'}, {"PRO", 'P'},
    {"SER", 'S'}, {"THR", 'T'}, {"TRP", 'W'}, {"TYR", 'Y'}, {"VAL", 'V'},
};

static std::map<char, std::string> const aminoAcidNameTable1to3 = { 
    {'A', "ALA"}, {'R', "ARG"}, {'N', "ASN"}, {'D', "ASP"}, {'C', "CYS"},
    {'E', "GLU"}, {'Q', "GLN"}, {'G', "GLY"}, {'H', "HIS"}, {'I', "ILE"},
    {'L', "LEU"}, {'K', "LYS"}, {'M', "MET"}, {'F', "PHE"}, {'P', "PRO"},
    {'S', "SER"}, {'T', "THR"}, {'W', "TRP"}, {'Y', "TYR"}, {'V', "VAL"},
};

static std::map<char, double> const massTable = {
    // https://www.thermofisher.com/de/de/home/references/ambion-tech-support/rna-tools-and-calculators/proteins-and-amino-acids.html
    {'A',  89.1}, {'R', 174.2}, {'N', 132.1}, {'D', 133.1}, {'C', 121.2},
    {'E', 147.1}, {'Q', 146.2}, {'G', 75.10}, {'H', 155.2}, {'I', 131.2},
    {'L', 131.2}, {'K', 146.2}, {'M', 149.2}, {'F', 165.2}, {'P', 115.1},
    {'S', 105.1}, {'T', 119.1}, {'W', 204.2}, {'Y', 181.2}, {'V', 117.1},
};

static std::map<char, double> const atomicVdwRadii = {
    // From https://github.com/getcontacts/getcontacts/blob/b0777f7148a327387133d9632f8ed1e34bfaa282/contact_calc/atom.py#L111
    // or https://en.wikipedia.org/wiki/Van_der_Waals_radius
    {'H', 0.12}, {'C', 0.17}, {'N', 0.155}, {'O', 0.152}, {'S', 0.18}
};

static std::map<char, double> const beadVdwRadii = {
    // https://doi.org/10.1088/1478-3975/12/4/046002
    // Mateusz Chwastyk, Adolfo Poma Bernaola and Marek Cieplak 
    // Statistical radii associated with amino acids to determine the contact map: fixing the structure of a type I cohesin domain in the Clostridium thermocellum cellulosome
    {'G', 3.15}, {'A', 3.35}, {'V', 4.00}, {'L', 4.60}, {'I', 4.50},
    {'M', 4.50}, {'P', 3.70}, {'H', 4.00}, {'F', 4.60}, {'Y', 4.50},
    {'W', 4.70}, {'C', 3.70}, {'S', 3.30}, {'T', 3.60}, {'N', 3.65},
    {'Q', 3.90}, {'D', 3.50}, {'E', 3.65}, {'K', 3.65}, {'R', 3.95}
};


enum ContactType {
    NONE = 0,
    OTHER = 1,
    HBOND = 2,
    SALTB = 3,
    DISULF = 4,
    HYDROPHOBIC = 5,
};

static std::map<ContactType, std::string> const contactTypeNames = {
    {NONE, "No contact"}, {OTHER, ""}, {HBOND, "H-bond"}, {SALTB, "Salt-bridge"}, {DISULF, "Disulfide bond"}, {HYDROPHOBIC, "Hydrophobic"}
};

static std::vector<char> nonDonorAcceptor = {'H', 'C', 'S'};