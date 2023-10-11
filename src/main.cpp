#include <iostream>
#include <vector>
#include <tuple>
#include <set>
#include <stack>
#include "eval.h"

# define BIGN 100000000

int* ref2pairs(std::string& ref){
    int* pairs = new int[ref.length()];
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
    int* pairs_all = ref2pairs(ref);
    std::vector<std::tuple<int, int>> pairs_diff;
    for(auto& idx: positions){
        if(pairs_all[idx]==idx)
            pairs_diff.push_back(std::make_tuple(idx, idx));
        else if (pairs_all[idx]>idx)
            pairs_diff.push_back(std::make_tuple(idx, pairs_all[idx]));
    }
    return pairs_diff;
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


void test(std::string& seq, std::string& ref1, std::string& ref2, bool is_verbose, int dangle_model) {
    std::cout << "seq: " << seq << std::endl;
    std::cout << "  y: " << ref1 << std::endl;
    std::cout << " y': " << ref2 << std::endl;
    // std::cout << "alg: " << alg << std::endl;
    std::set<int> critical_positions;
    // std::vector<std::vector<int>> cr_loops = find_critical(ref1, ref2, is_verbose);
    std::vector<std::vector<int>> cr_loops = find_critical_plus(ref1, ref2, critical_positions, is_verbose);
    // long delta_energy = diff_eval(seq, cr_loops, is_verbose, dangle_model);
    long delta_energy = max_diff(seq.size(), cr_loops, is_verbose, dangle_model);
    if (is_verbose){
        printf("delta  : %.2f kcal/mol\n", delta_energy/-100.0);
        printf("critical positions: ");
        for (int x: critical_positions)
            printf("%d, ", x);
        printf("\n");
    }
    std::vector<std::tuple<int, int>> pairs_diff = idx2pair(critical_positions, ref1);
    if (is_verbose)
        for(auto& pair: pairs_diff)
            std::cout<<std::get<0>(pair)<<"\t"<<std::get<1>(pair)<<std::endl;
    ulong n_enum = count_enum(pairs_diff);
    std::cout<<"enumeration count: "<<n_enum<<std::endl;
    if(n_enum < BIGN){
        std::cout<<"alg 1"<<std::endl;
    }
}

int main(int argc, char* argv[]) {
    char* alg = (argc < 2) ? nullptr : argv[1];
    std::cout << "alg: " << ((alg == nullptr) ? "None" : alg) << std::endl;
    bool verbose = true;
    int dangle = 2;
    while (true) {
        std::string seq;
        std::string ref1;
        std::string ref2;

        std::string line;

        // Read input line by line until EOF (end of file) is reached
        while (std::getline(std::cin, seq)) {
            // Process the line as needed
            // std::cout<<"got seq"<<std::endl;
            getline(std::cin, ref1);
            getline(std::cin, ref2);
            test(seq, ref1, ref2, verbose, dangle);
        }
        break;
    }
    return 0;
}