#include <iostream>
#include <vector>
#include <tuple>
#include <set>
#include <map>
#include <stack>
#include <ctime>
#include <chrono>
#include <string.h> 

#include "eval.h"
#include "fold.h"
// #include "LinearFold.h"

using namespace std;

# define BIGN 100000000

char nuc_all[] = "ACGU";
char nuc_pair_all[][3] = {"GC", "CG", "AU", "UA", "GU", "UG"};
std::string cs_fold(std::string seq, std::string& constr, int beamsize, bool sharpturn, bool verbose, int dangle);
ulong count_enum(std::vector<std::tuple<int, int>>& pairs_diff);
std::string enumerate(std::vector<std::tuple<int, int>>& pairs_diff, ulong order, std::string& seq);
// Define a basic TreeNode structure
struct TreeNode {
    int first; // root: -1
    int second; // 
    // Bracket br;
    std::vector<TreeNode*> children;

    TreeNode(int first_val, int second_val){
        first = first_val;
        second = second_val;
    }

    void printTree(){
        printf("first: %d, second: %d\n", first, second);
        for(int i = 0; i < children.size(); i++){
            printf("child[%d]: first: %d, second: %d\n", i, children[i]->first, children[i]->second);
        }
        printf("\n");
        for(int j = 0; j < children.size(); j++){
            children[j]->printTree();
        }       
    }

    void printTreeEnum(std::string& seq, std::string& y){
        printf("first: %d, second: %d\n", first, second);
        std::vector<std::tuple<int, int>> pairs_diff;
        int n_unpair = second - first - 1;
        int n_pair = 0;
        if (first >= 0){
            n_pair += 1;
            pairs_diff.push_back(std::make_tuple(first, second));
        }
        for(int i = 0; i < children.size(); i++){
            n_unpair -= children[i]->second - children[i]->first + 1;
            pairs_diff.push_back(std::make_tuple(children[i]->first, children[i]->second));
            n_pair += 1;
        }
        ulong count = 1;
        for(int p=0; p<n_pair; p++)
            count *= 6;
        for(int u=0; u<n_unpair; u++)
            count *= 4;
        printf("(%d, %d), pairs: %d, unpairs: %d, enum: %d\n", first, second, n_pair, n_unpair, count);
        std::string subref(seq.length(), '?');
        if (first >= 0){
            int len = second-first+1;
            // std::string subseq = seq.substr(first, len);
            subref[first] = '(';
            subref[second] = ')';
            for(int i = 0; i < children.size(); i++){
                subref[children[i]->first] = '(';
                subref[children[i]->second] = ')';
                for (int k = children[i]->first+1; k < children[i]->second; k++)
                    subref[k] = y[k];
            }
            // printf("constraint: %s\n", subref.c_str());
            std::cout<<subref.substr(first, len)<<std::endl;
            for(int l = 0; l < len; l++){
                if (subref[first+l]=='?')
                    pairs_diff.push_back(std::make_tuple(first+l, first+l));
            }
        }else{
            int len = second-first-1;
            assert( len == seq.length() );
            for(int i = 0; i < children.size(); i++){
                subref[children[i]->first] = '(';
                subref[children[i]->second] = ')';
                for (int k = children[i]->first+1; k < children[i]->second; k++)
                    subref[k] = y[k];
            }
            // printf("constraint: %s\n", subref.c_str());
            std::cout<<subref<<std::endl;
            for(int l = 0; l < len; l++){
                if (subref[l]=='?')
                    pairs_diff.push_back(std::make_tuple(l, l));
            }
        }
        ulong count_v2 =  count_enum(pairs_diff);
        printf("enum v2: %d\n", count_v2);
        assert(count_v2 == count);
        if (first >= 0){
            std::string constr_i = subref.substr(first, second-first+1);
            std::string gold_i = constr_i;
            std::replace(gold_i.begin(), gold_i.end(), '?', '.');
            int count_mfe = 0;
            std::cout<<"gold: "<<gold_i<<std::endl;
            for(ulong i=0; i<count_v2; i++){
                std::string seq_i = enumerate(pairs_diff, i, seq);
                std::string subseq_i = seq_i.substr(first, second-first+1);
                std::string mfe_str = cs_fold(subseq_i, constr_i, 0, false, false, 2);
                if(mfe_str==gold_i)
                    count_mfe += 1;
            }
            printf("count_mfe: %d\n", count_mfe);
        }
        for(int j = 0; j < children.size(); j++){
            children[j]->printTreeEnum(seq, y);
        }       
    }

    void printTree(std::string ref){
        printf("first: %d, second: %d\n", first, second);
        if(first >= 0)
            printf("%s\n", ref.substr(first, second-first+1).c_str());
        else
            printf("%s\n", ref.c_str());
        for(int i = 0; i < children.size(); i++){
            printf("child[%d]: first: %d, second: %d\n", i, children[i]->first, children[i]->second);
        }
        printf("\n");
        for(int j = 0; j < children.size(); j++){
            children[j]->printTree(ref);
        }       
    }
};

// Function to parse a string of nested pairs into a tree
TreeNode* parseStringToTree(const std::string& ref) {
    std::stack<TreeNode*> nodeStack;
    TreeNode* root = new TreeNode(-1, ref.length());
    nodeStack.push(root);
    for (int i = 0; i < ref.length(); i++) {
        char c = ref[i];
        if (c == '(') {
            TreeNode* node = new TreeNode(i, -2);
            nodeStack.top()->children.push_back(node);
            nodeStack.push(node);
        } else if (c == ')') {
            nodeStack.top()->second = i;
            // printf("pop: (%d, %d), %d\n", nodeStack.top()->first, nodeStack.top()->second, nodeStack.top()->children.size());
            nodeStack.pop();
        }
    }
    // printf("root: (%d, %d)\n", nodeStack.top()->first, nodeStack.top()->second);
    return root;
}


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
                // printf("i=%d, j=%d, %s, %d\n", i, j, nuc_ij, strcmp(nuc_ij, "GC"));
                return false;
            }
        }
    }
    return true;
}

std::vector<std::string> alg_1(std::string& y, std::string& y_prime, std::vector<std::vector<int>>& cr_loops, std::vector<std::tuple<int, int>>& pairs_diff, std::string& seq, bool is_verbose, int dangle_model){
    time_t startTime = time(nullptr);
    auto start = std::chrono::high_resolution_clock::now();
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
            auto pause = std::chrono::high_resolution_clock::now();
            printf("%8d, %d seconds\n", i+1, std::chrono::duration_cast<std::chrono::milliseconds>(pause - start));
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
    auto stop = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double, std::milli> fp_ms = stop - start;
    printf("finished: %f seconds\n", fp_ms/1000.f);
    // printf("finished: %f seconds\n", std::chrono::duration_cast<std::chrono::milliseconds>(stop - start)/1000.f);
    printf("seqs size: %d\n", seqs.size());
    return X;
}

std::string cs_fold(std::string seq, std::string& constr, int beamsize, bool sharpturn, bool verbose, int dangle){
    bool consflag = true;
    std::set<char> consSet {'?', '.', '(', ')'};
    if (seq.length() != constr.length()){
        printf("The lengths don't match between sequence and constraints: %s, %s\n", seq.c_str(), constr.c_str());
        return "";
    }
    int n = seq.length();
    std::vector<int> cons(n);
    std::stack<int> leftBrackets;
    consflag = true;
    for (int i=0; i < n; i++){
        char coni = constr[i];
        if (consSet.count(coni) == 0){
            printf("Unrecognized constraint character, should be ? . ( or )\n");
            consflag = false;
            break;
        }
        switch(coni){
            case '.':
                cons[i] = -2;
                break;
            case '?':
                cons[i] = -1;
                break;
            case '(':
                leftBrackets.push(i);
                break;
            case ')':
                int leftIndex = leftBrackets.top();
                leftBrackets.pop();
                cons[leftIndex] = i;
                cons[i] = leftIndex;
                break;
        }
    }
    if (consflag) {
        // printf("%s\n", constr.c_str());
        
        // lhuang: moved inside loop, fixing an obscure but crucial bug in initialization
        BeamCKYParser parser(beamsize, !sharpturn, verbose, true, false, 5.0, "", false, dangle);

        BeamCKYParser::DecoderResult result = parser.parse(seq, &cons);

        #ifdef lv
                double printscore = (result.score / -100.0);
        #else
                double printscore = result.score;
        #endif
                // printf("%s (%.2f)\n", result.structure.c_str(), printscore);
        return result.structure;
    }
    return "";
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
    if ( alg != nullptr && strcmp(alg, "csfold") == 0 ){
        std::cout << alg << std::endl;
        int beamsize = 0;
        bool sharpturn = false;
        std::string seq;
        std::string constr;

        while (std::getline(std::cin, seq))
        {
            std::getline(std::cin, constr);
            // printf("%s", constr.c_str());
            cs_fold(seq, constr, beamsize, sharpturn, verbose, dangle);
        }
        return 0;
    }else if (alg != nullptr && strcmp(alg, "inout") == 0){
        std::cout << alg << std::endl;
        int beamsize = 0;
        bool sharpturn = false;
        std::string seq;
        std::string ref;
        // std::string constr;
        while (std::getline(std::cin, seq))
        {
            std::getline(std::cin, ref);
            TreeNode* root = parseStringToTree(ref);
            root->printTree(ref);
        }
        return 0;
    }else{
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
    }
    return 0;
}