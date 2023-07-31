#pragma once
#include <string>
#include <vector>
#include <fstream>
#include <sys/stat.h>
#include <limits>
#include <math.h>

bool dirExists(const std::string path);
int createDir(std::string s, mode_t mode);
void tabulateAngles(const std::vector<double>& angles, const std::string& directoryName, const std::string& filename, const unsigned int N);
void tabulateDihedrals(const std::vector<double>& angles, const std::string& directoryName, const std::string& filename, const unsigned int multiplicity, const unsigned int N);
