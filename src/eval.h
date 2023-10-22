#ifndef EVAL_H
#define EVAL_H

#include <vector>
#include <set>
#include <string>

// Function to find critical values
std::vector<std::vector<int>> find_critical(std::string ref1, std::string ref2, bool is_verbose);
std::vector<std::vector<int>> find_critical_plus(std::string ref1, std::string ref2, std::set<int>& critical_positions, bool is_verbose);

// Function for diff evaluation
long diff_eval(std::string seq, std::vector<std::vector<int>>& cr_loops, bool is_verbose, int dangle_model);
long max_diff(int n, std::vector<std::vector<int>>& cr_loops, bool is_verbose, int dangle_model);
long linear_eval(std::string& seq, std::string& ref, bool& is_verbose, int& dangle_model);


#endif // EVAL_H
