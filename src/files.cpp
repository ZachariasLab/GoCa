#include "files.hpp"
#include "helper.hpp"

/** Directory helper methods **/
bool dirExists(const std::string path) {
    struct stat info;
    if (stat(path.c_str(), &info ) != 0)
        return false;
    else if(info.st_mode & S_IFDIR)
        return true;
    else
        return false;
};

int createDir(std::string s, mode_t mode) {
    size_t pos = 0;
    std::string dir;
    int mdret;
    if(s[s.size()-1] != '/'){
        s += '/';
    }
    while((pos=s.find_first_of('/',pos)) != std::string::npos){
        dir = s.substr(0, pos++);
        if (dir.size() == 0) continue;
        if ((mdret = mkdir(dir.c_str(),mode)) && errno != EEXIST) {
            return mdret;
        }
    }
    return mdret;
};

/** Tabulation of angles and dihedrals **/
#define STRANGE_ANGLE_FACTOR 0.000152295524691358
void tabulateAngles(const std::vector<double>& angles, const std::string& directoryName, const std::string& filename, const unsigned int N) {
    if (!dirExists(directoryName)) {
        createDir(directoryName, 0777);
    }
    std::ofstream file(directoryName + "/" + filename);
    file << "@ title \"Tabulated angle interaction\"\n";
    for (double x = 0; x <= 180; x += 180.0 / (N-1)) {
        double fx = std::numeric_limits<double>::max();
        double dfx = 0;
        for (auto it = angles.begin(); it != angles.end(); it++) {
            double value = pow(x - *it, 2) * STRANGE_ANGLE_FACTOR;
            if (value < fx) {
                fx = value;
                dfx = -2 * (x - *it) * STRANGE_ANGLE_FACTOR;
            }
        }
        file << string_format("%16.9E %16.9E %16.9E\n", x, fx, dfx);
    }
};

void tabulateDihedrals(const std::vector<double>& angles, const std::string& directoryName, const std::string& filename, const unsigned int multiplicity, const unsigned int N) {
    /** Create tabulation directory if we have multiple models or chain merging enabled **/
    if (!dirExists(directoryName)) {
        createDir(directoryName, 0777);
    }
    std::ofstream file(directoryName + "/" + filename);
    file << "@ title \"Tabulated dihedral interaction\"\n";
    for (double x = -180; x <= 180; x += 360.0 / (N-1)) {
        double fx = std::numeric_limits<double>::max();
        double dfx = 0;
        for (auto it = angles.begin(); it != angles.end(); it++) {
            double value = 1 + cos(multiplicity * (x - *it) * M_PI/180.0);
            if (value < fx) {
                fx = value;
                // Two minus signs cancel each other for the derivative
                dfx = multiplicity * sin(multiplicity * (x - *it) * M_PI/180.0) * M_PI/180.0;
            }
        }
        file << string_format("%16.9E %16.9E %16.9E\n", x, fx, dfx);
    }
};
