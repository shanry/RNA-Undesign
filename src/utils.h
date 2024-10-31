
#ifndef UTILS_H
#define UTILS_H

#include <set>
#include <vector>
#include <map>
#include <unordered_map>
#include <utility>
#include <string>
#include <random>
#include <cmath>
// #include <ctime>

/* Old compatibility names for C types.  */
typedef unsigned long int ulong;

bool compareByFirstStringLength(const std::pair<std::string, std::string> &a, const std::pair<std::string, std::string> &b); 
std::vector<int> findAllOccurrences(const std::string& mainString, const std::string& subString); 
std::vector<std::string> readLinesFromFile(const std::string& filePath); 
bool isMFE(std::vector<std::string>& subopts, std::string& target);
bool isUMFE(std::vector<std::string>& subopts, std::string& target);
std::vector<std::string> split_string(const std::string& input, char delimiter); 
std::vector<int> ref2pairs(std::string& ref);
std::vector<std::tuple<int, int>> idx2pair(std::set<int>& positions, std::string& ref);
std::set<std::pair<int, int>> ref2pairset(std::string& ref);
ulong count_enum(std::vector<std::tuple<int, int>>& pairs_diff);
bool check_compatible(std::string seq, std::string ss);
template <typename T>
std::set<T> setIntersection(const std::set<T>& set1, const std::set<T>& set2); 
std::string getSubstrings(std::set<int>& indices, std::string& str);
std::string exec_command(const char* cmd); 
std::vector<std::string> subopt(std::string& seq, std::string constr);
std::string tg_init(std::string& y, int seed = 0);
int countOccurrences(const std::string& str, char target); 
std::string compose_args4plot(std::string id, std::string y, std::vector<std::pair<int, int>>& pairs_outside, std::vector<std::pair<int, int>>& pairs_inside);
std::string compose_pairsplot(std::string id, std::string y, std::vector<std::pair<int, int>>& pairs_ds, std::vector<std::pair<int, int>>& pairs_ud);
std::string compose_pairstr(std::vector<std::pair<int, int>>& pairs_inside, std::vector<std::pair<int, int>>& pairs_outside);
std::string fl2str(float x, int d=4);
std::string replaceFileExtension(const std::string& filePath, const std::string& newExtension);
std::map<std::string, std::set<std::string>> readMotif(const char* file);
std::vector<std::vector<int>> PowerSet3(int start, int set_size);
std::vector<std::vector<int>> PowerSet(int start, int set_size, int subset_size);
std::vector<std::vector<std::pair<int, int>>> pairSubSet(std::vector<std::pair<int, int>> pairSet);
std::vector<std::vector<std::string>> read_csv(const char* file);
std::unordered_map<std::string, std::string> loadlib_eterna(std::string csv);
std::string pairs2string(std::vector<std::pair<int, int>> pairs);
std::string genHelix(int len);
inline int countDotBrackets(const std::string& text) {
    // Count occurrences of '(', '.', and ')'
    int count = std::count(text.begin(), text.end(), '(') +
                std::count(text.begin(), text.end(), '.') +
                std::count(text.begin(), text.end(), ')');
    return count;
}
inline std::string getRandomIntString(int min = std::pow(10, 8), int max = std::pow(10, 9) - 1) {
    static std::random_device rd;       // Seed generator (only initialized once)
    static std::mt19937 gen(rd());      // Mersenne Twister engine
    std::uniform_int_distribution<> dist(min, max);
    return std::to_string(dist(gen));   // Convert the integer to a string
}
#endif