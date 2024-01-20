// #ifndef UTILS
// #define UTILS

#include <iostream>
#include <fstream>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <sstream>
#include <array>
#include <vector>
#include <string>
#include <cstring>
#include <utility>
#include <set>
#include <map>
#include <ctime>
#include <stack>
#include <cassert>
#include <cmath>
#include <algorithm>


bool compareByFirstStringLength(const std::pair<std::string, std::string> &a, const std::pair<std::string, std::string> &b) {
    return a.first.length() < b.first.length();
}

std::vector<int> findAllOccurrences(const std::string& mainString, const std::string& subString) {
    std::vector<int> pos_vec;
    size_t pos = mainString.find(subString); // Find the first occurrence
    while (pos != std::string::npos) {
        // std::cout << "Substring found at index: " << pos << std::endl;
        pos_vec.push_back(pos);
        pos = mainString.find(subString, pos + 1); // Find the next occurrence
    }
    return pos_vec;
}

std::vector<std::string> readLinesFromFile(const std::string& filePath) {
    std::vector<std::string> lines;
    std::string line;
    
    std::ifstream file(filePath);

    if (file.is_open()) {
        while (std::getline(file, line)) {
            lines.push_back(line);
        }
        file.close();
    } else {
        std::cerr << "Error: Unable to open the file " << filePath << std::endl;
    }

    return lines;
}

bool isMFE(std::vector<std::string>& subopts, std::string& target){
    auto it = std::find(subopts.begin(), subopts.end(), target);
    if (it != subopts.end())
        return true;
    else
        return false;
}

bool isUMFE(std::vector<std::string>& subopts, std::string& target){
    auto it = std::find(subopts.begin(), subopts.end(), target);
    if (it != subopts.end() and subopts.size() == 1)
        return true;
    else
        return false;
}

std::vector<std::string> split_string(const std::string& input, char delimiter) {
    std::vector<std::string> tokens;
    std::istringstream tokenStream(input);
    std::string token;
    
    while (std::getline(tokenStream, token, delimiter)) {
        tokens.push_back(token);
    }
    
    return tokens;
}

std::vector<int> ref2pairs(std::string& ref){
    std::vector<int> pairs(ref.length(), -1);
    std::stack<int> brackets;
    for(int i = 0; i < ref.length(); i++){
        if(ref[i] == '.')
            pairs[i] = i;
        else if (ref[i] == '(')
            brackets.push(i);
        else{
            int j = brackets.top();
            pairs[j] = i;
            pairs[i] = j;
            brackets.pop();
        }
    }
    return pairs;
}

std::vector<std::tuple<int, int>> idx2pair(std::set<int>& positions, std::string& ref){
    std::vector<int> pairs_all = ref2pairs(ref);
    std::vector<std::tuple<int, int>> pairs_diff;
    for(auto& idx: positions){
        if(pairs_all[idx]==idx)
            pairs_diff.push_back(std::make_tuple(idx, idx));
        else if (pairs_all[idx]>idx)
            pairs_diff.push_back(std::make_tuple(idx, pairs_all[idx]));
    }
    return pairs_diff;
}

std::set<std::pair<int, int>> ref2pairset(std::string& ref){
    std::set<std::pair<int, int>> pairset;
    std::stack<int> brackets;
    for(int i = 0; i < ref.length(); i++){
        if (ref[i] == '(')
            brackets.push(i);
        else if (ref[i] == ')'){
            int j = brackets.top();
            pairset.insert(std::make_pair(j, i));
            brackets.pop();
        }
    }
    return pairset;
}

ulong count_enum(std::vector<std::tuple<int, int>>& pairs_diff){
    ulong count = 1;
    for(auto& pair: pairs_diff){
        if (std::get<0>(pair) == std::get<1>(pair))
            count *= 4;
        else
            count *= 6;
    }
    return count;
}

bool check_compatible(std::string seq, std::string ss){
    std::vector<int> pairs_list = ref2pairs(ss);
    for(int i=0; i<seq.length(); i++){
        int j = pairs_list[i];
        if(i!=j){
            char nuc_ij[3];
            nuc_ij[0] = seq[i];
            nuc_ij[1] = seq[j];
            nuc_ij[2] = '\0';
            if( strcmp(nuc_ij, "CG")&&strcmp(nuc_ij, "GC")&&strcmp(nuc_ij, "AU")&&strcmp(nuc_ij, "UA")&&strcmp(nuc_ij, "GU")&&strcmp(nuc_ij, "UG") ){
                return false;
            }
        }
    }
    return true;
}

template <typename T>
std::set<T> setIntersection(const std::set<T>& set1, const std::set<T>& set2) {
    std::set<T> intersection;
    std::set_intersection(set1.begin(), set1.end(), set2.begin(), set2.end(), std::inserter(intersection, intersection.begin()));
    return intersection;
}

std::string getSubstrings(std::set<int>& indices, std::string& str){
    std::string substr;
    for (int index : indices) {
        substr += str[index];
    }
    return substr;
}

std::string exec_command(const char* cmd) {
    std::array<char, 1280> buffer;
    std::string result;
    std::shared_ptr<FILE> pipe(popen(cmd, "r"), pclose);
    
    if (!pipe) {
        throw std::runtime_error("popen() failed!");
    }
    
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        result += buffer.data();
    }
    
    return result;
}

std::vector<std::string> subopt(std::string& seq, std::string constr){
    // std::string cmd_str = "echo -ne" + "\"" + seq + "\\n" + constr + "\"" + "|" + "/nfs/stak/users/zhoutian/acl/rna_tool/v251/bin/RNAsubopt -C --enforceConstraint -e 0";
    std::replace(constr.begin(), constr.end(), '.', 'x');
    std::replace(constr.begin(), constr.end(), '?', '.');
    std::string cmd_str = "echo -ne \"" + seq + "\\n" + constr + "\" | /nfs/stak/users/zhoutian/acl/rna_tool/v251/bin/RNAsubopt -C --enforceConstraint -e 0";
    const char* cmd_cstr = cmd_str.c_str();
    // printf("%s\n", cmd_cstr);
    std::string output = exec_command(cmd_cstr);
    std::vector<std::string> array = split_string(output, '\n');
    std::vector<std::string> structures;
    for(int i = 1; i < array.size(); i++){
        structures.push_back(split_string(array[i], ' ')[0]);
    }
    return structures;
}

std::string tg_init(std::string& y){
    std::string x(y.length(), 'A');
    std::vector<int> pair_map = ref2pairs(y);
    assert(pair_map.size() == y.length());
    // Seed the random number generator with the current time
    std::srand(static_cast<unsigned int>(std::time(nullptr)));
    for(int i = 0; i < pair_map.size(); i ++){
        int j = pair_map[i];
        if(i < j){
            if(rand()%2){
                x[i] = 'C';
                x[j] = 'G';
            }else{
                x[i] = 'G';
                x[j] = 'C';
            }
        }
    }
    return x;
}

int countOccurrences(const std::string& str, char target) {
    int count = 0;
    for (char c : str) {
        if (c == target) {
            count++;
        }
    }
    return count;
}

std::string compose_args4plot(std::string id, std::string y, std::vector<std::pair<int, int>>& pairs_outside, std::vector<std::pair<int, int>>& pairs_inside){
    std::string prestring;
    // color the outside span
    if(pairs_outside[0].first >= 0)
        prestring += std::to_string(pairs_outside[0].first+1) + " " + std::to_string(pairs_outside[0].second+1) + " GREEN Fomark ";
    else
        prestring += std::to_string(1) + " " + std::to_string(y.length()) + " GREEN Fomark ";
    // uncolor inside spans
    for(int i = 1; i < pairs_outside.size(); i++){
        prestring += std::to_string(pairs_outside[i].first+1) + " " + std::to_string(pairs_outside[i].second+1) + " WHITE Fomark ";
    }
    // color outside pairs
    for(int i = 0; i < pairs_outside.size(); i++){
        std::pair<int, int> pair = pairs_outside[i];
        if(i == 0){
            if (pair.first >=0)
                prestring += std::to_string(pair.first+1) + " " + std::to_string(pair.second+1) + " 0.70 0.5 colorpair ";
        }else{
            prestring += std::to_string(pair.first+1) + " " + std::to_string(pair.second+1) + " 0.667 0.5 colorpair ";
        }
    }
    // color inside pairs
    for(std::pair<int, int> pair: pairs_inside)
        prestring += std::to_string(pair.first+1) + " " + std::to_string(pair.second+1) + " 0.1667 1.0 colorpair ";
        
    return "\"" + id + "," + y + "," + prestring + "\"";
}

std::string compose_pairsplot(std::string id, std::string y, std::vector<std::pair<int, int>>& pairs_ds, std::vector<std::pair<int, int>>& pairs_ud){
    std::string prestring;
    // color designable pairs
    for(int i = 0; i < pairs_ds.size(); i++){
        std::pair<int, int> pair = pairs_ds[i];
        prestring += std::to_string(pair.first+1) + " " + std::to_string(pair.second+1) + " 0.0 0.5 colorpair ";
    }
    // color undesignable pairs
    for(std::pair<int, int> pair: pairs_ud)
        prestring += std::to_string(pair.first+1) + " " + std::to_string(pair.second+1) + " 0.6667 0.5 colorpair ";

    return "\"" + id + "," + y + "," + prestring + "\"";
}

std::string compose_pairstr(std::vector<std::pair<int, int>>& pairs_inside, std::vector<std::pair<int, int>>& pairs_outside){
    std::string pairstring;
    // color designable pairs
    pairstring += " inpairs:";
    for(auto pair: pairs_inside){
        pairstring += std::to_string(pair.first) + "," + std::to_string(pair.second) + ";";
    }
    // color undesignable pairs
    pairstring += " outpairs:";
    for(auto pair: pairs_outside)
        pairstring += std::to_string(pair.first) + "," + std::to_string(pair.second) + ";";

    return pairstring;
}

std::string fl2str(float x, int d=4){
    // Set the precision to `d` decimals and convert float to string
    std::string str = std::to_string(x);
    size_t found = str.find(".");
    if (found != std::string::npos) {
        size_t decimals = str.size() - found - 1;
        if (decimals > d) {
            decimals = d; // Limit to 4 decimals
        }
        str = str.substr(0, found + 1 + decimals);
    }
    return str;
}

// Function to replace file extension
std::string replaceFileExtension(const std::string& filePath, const std::string& newExtension) {
    size_t lastDotPos = filePath.find_last_of('.');
    
    if (lastDotPos != std::string::npos) {
        std::string stem = filePath.substr(0, lastDotPos);
        return stem + "." + newExtension;
    }

    // If there's no dot in the original file path, just append the new extension
    return filePath + "." + newExtension;
}

// get pairset
// get pairset
std::map<std::string, std::set<std::string>> readMotif(const char* file){
    std::map<std::string, std::set<std::string>> id2motif;
    std::ifstream inputFile(file); // Replace "example.txt" with the name of your file
    if (!inputFile.is_open()) {
        std::cerr << "Failed to open the file." << std::endl;
        return id2motif;
    }

    std::string line;

    while (std::getline(inputFile, line)){
        std::vector<std::string> splits = split_string(line, ':');
        std::string key = splits[0];
        std::string val = splits[1];
        // Check if the key already exists
        auto it = id2motif.find(key);
        if (it == id2motif.end()) {
            // Key doesn't exist, create a new entry
            id2motif[key] = {val};
        } else {
            // Key exists, add the pair to the existing entry
            id2motif[key].insert(val);
        }
    }

    return id2motif;
}

std::vector<std::vector<int>> PowerSet3(int start, int set_size){ 
    std::vector<std::vector<int>> powset_all;
    // Set_size of power set of a set with set_size 
    // n is (2^n-1) 
    unsigned int pow_set_size = pow(2, set_size); 
    int counter, j; 
  
    // Run from counter 000..0 to 111..1 
    for (counter = 0; counter < pow_set_size; counter++) { 
        std::vector<int> powset;
        for (j = 0; j < set_size; j++) { 
            // Check if jth bit in the counter is set 
            // If set then print jth element from set 
            if (counter & (1 << j)) 
                powset.push_back(j+start);
        } 
        if(powset.size() > 2 && powset.size() < 5)
            powset_all.push_back(powset);
    } 
    return powset_all;
}

// #endif