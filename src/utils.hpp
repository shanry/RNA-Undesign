#ifndef UTILS
#define UTILS

#include <iostream>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <sstream>
#include <array>
#include <vector>
#include <string>
#include <utility>
#include <set>
#include <algorithm>

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

#endif