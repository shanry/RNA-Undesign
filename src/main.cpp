#include <iostream>
#include <vector>
#include <tuple>
#include <set>
#include <stack>
#include <ctime>

#include <string.h> 
#include "eval.h"

# define BIGN 100000000

char nuc_all[] = "ACGU";
char nuc_pair_all[][3] = {"GC", "CG", "AU", "UA", "GU", "UG"};

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

std::string enumerate(std::vector<std::tuple<int, int>>& pairs_diff, ulong order, std::string& seq){
    std::string seq_new = seq;
    for(auto& pair: pairs_diff){
        int first = std::get<0>(pair);
        int second = std::get<1>(pair);
        if (first == second){
            int idx_nuc = order%4;
            seq_new[second] = nuc_all[idx_nuc];
            order = order>>2;
        }else{
            int idx_nuc_pair = order%6;
            seq_new[first] = nuc_pair_all[idx_nuc_pair][0];
            seq_new[second] = nuc_pair_all[idx_nuc_pair][1];
            order = order/6;
        }
    }
    return seq_new;
}

bool check_compatible(std::string seq, std::string ss){
    int* pairs_list = ref2pairs(ss);
    for(int i=0; i<seq.length(); i++){
        int j = pairs_list[i];
        if(i!=j){
            char nuc_ij[3];
            nuc_ij[0] = seq[i];
            nuc_ij[1] = seq[j];
            nuc_ij[2] = '\0';
            if( strcmp(nuc_ij, "CG")&&strcmp(nuc_ij, "GC")&&strcmp(nuc_ij, "AU")&&strcmp(nuc_ij, "UA")&&strcmp(nuc_ij, "GU")&&strcmp(nuc_ij, "UG") ){
                printf("i=%d, j=%d, %s, %d\n", i, j, nuc_ij, strcmp(nuc_ij, "GC"));
                return false;
            }
        }
    }
    return true;
}

std::vector<std::string> alg_1(std::string& y, std::string& y_prime, std::vector<std::vector<int>>& cr_loops, std::vector<std::tuple<int, int>>& pairs_diff, std::string& seq, bool is_verbose, int dangle_model){
    time_t startTime = time(nullptr);
    // std::cout << "Current time in seconds: " << startTime << std::endl;
    ulong n_enum = count_enum(pairs_diff);
    std::vector<std::string> X;
    std::set<std::string> seqs;
    for(ulong i=0; i<n_enum; i++){
        std::string seq_i = enumerate(pairs_diff, i, seq);
        // printf("%8d, %s\n", i, seq_i.c_str());
        seqs.insert(seq_i);
        if ((i+1)%10000 == 0){
            time_t pauseTime = time(nullptr);
            printf("%8d, %d seconds\n", i+1, pauseTime-startTime);
        }
        if(check_compatible(seq_i, y_prime)){
            // double e_diff = diff_eval(seq_i, cr_loops, is_verbose, dangle_model)/-100.0;
            double energy_ref1 = linear_eval(seq_i, y, is_verbose, dangle_model)/-100.0;
            double energy_ref2 = linear_eval(seq_i, y_prime, is_verbose, dangle_model)/-100.0;
            double e_diff = energy_ref1 - energy_ref2;
            // printf("e1: %.2f\n", energy_ref1);
            // printf("e2: %.2f\n", energy_ref2);
            if(e_diff <= 0.)
                X.push_back(seq_i);
        }else{
            X.push_back(seq_i);
        }
    }
    printf("seqs size: %d\n", seqs.size());
    return X;
}



void test(std::string& seq, std::string& ref1, std::string& ref2, bool is_verbose, int dangle_model) {
    std::cout << "seq: " << seq << std::endl;
    std::cout << "  y: " << ref1 << std::endl;
    std::cout << " y': " << ref2 << std::endl;
    // std::cout << "alg: " << alg << std::endl;
    std::set<int> critical_positions;
    // std::vector<std::vector<int>> cr_loops = find_critical(ref1, ref2, is_verbose);
    std::vector<std::vector<int>> cr_loops = find_critical_plus(ref1, ref2, critical_positions, is_verbose);
    long delta_energy = diff_eval(seq, cr_loops, is_verbose, dangle_model);
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
        auto X = alg_1(ref1, ref2, cr_loops, pairs_diff, seq, is_verbose, dangle_model);
        printf("X size: %d\n", X.size());
    }
}

int main(int argc, char* argv[]) {
    char* alg = (argc < 2) ? nullptr : argv[1];
    std::cout << "alg: " << ((alg == nullptr) ? "None" : alg) << std::endl;
    bool verbose = false;
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