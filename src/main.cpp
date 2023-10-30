#include <iostream>
#include <fstream>
#include <vector>
#include <tuple>
#include <set>
#include <map>
#include <stack>
#include <ctime>
#include <chrono>
#include <string.h> 
#include <cmath>
#include <utility>
#include <algorithm>
#include <omp.h>

#include "eval.cpp"

using namespace std;

#define MAX_ENUM 10000000000
#define MAX_CONSTRAINT 100000

char nuc_all[] = "ACGU";
char nuc_pair_all[][3] = {"GC", "CG", "AU", "UA", "GU", "UG"};
std::vector<std::string> TriHP{"", "CAACG", "GUUAC"}; // middle position can't be C or G
std::vector<std::string> TetraHP{"",
                                 "CAACGG",
                                 "CCAAGG",
                                 "CCACGG",
                                 "CCCAGG",
                                 "CCGAGG",
                                 "CCGCGG",
                                 "CCUAGG",
                                 "CCUCGG",
                                 "CUAAGG",
                                 "CUACGG",
                                 "CUCAGG",
                                 "CUCCGG",
                                 "CUGCGG",
                                 "CUUAGG",
                                 "CUUCGG",
                                 "CUUUGG"}; // middle positions can't be AG, AU, CG, CU, GG, GU, UG
// middle positions can be AC(3), AA(2), CA(2), CC(1), GA(1), GC(2), UA(2), UC(2), UU(1)

std::vector<std::string> cs_fold(std::string seq, std::string& constr, int beamsize, bool sharpturn, bool verbose, int dangle);
std::vector<std::string> fold(std::string seq, int beamsize, bool sharpturn, bool verbose, int dangle, float energy_delta);
ulong count_enum(std::vector<std::tuple<int, int>>& pairs_diff);
std::string enumerate(std::vector<std::tuple<int, int>>& pairs_diff, ulong order, std::string& seq);

// Define a basic constraint structure
struct Constraint {
    std::set<int>* indices;
    std::vector<std::string>* seqs;
    std::string structure;

    // Constructor that accepts two arguments to initialize the members
    Constraint(std::set<int>* index_set, std::vector<std::string>* seq_list){
        indices = index_set;
        seqs = seq_list;
    }

    void setStructure(std::string& ref){
        structure = ref;
    }

    // ~Constraint(){
    //     delete indices;
    //     delete seqs;
    // }
};

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
                std::string mfe_str = cs_fold(subseq_i, constr_i, 0, false, false, 2)[0];
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

    void printTree(std::string& ref, std::string& seq){
        printf("first: %d, second: %d\n", first, second);
        if(first >= 0){
            printf("%s\n", seq.substr(first, second-first+1).c_str());
            printf("%s\n", ref.substr(first, second-first+1).c_str());
        }
        else{
            printf("%s\n", seq.c_str());
            printf("%s\n", ref.c_str());
        }
        for(int i = 0; i < children.size(); i++){
            printf("child[%d]: first: %d, second: %d\n", i, children[i]->first, children[i]->second);
        }
        printf("\n");
        for(int j = 0; j < children.size(); j++){
            children[j]->printTree(ref, seq);
        }       
    }

    void printTree(std::string& ref, std::string& seq, std::vector<std::pair<std::string, std::string>>& subrefs){
        printf("first: %d, second: %d\n", first, second);
        if(first >= 0){
            printf("%s\n", seq.substr(first, second-first+1).c_str());
            printf("%s\n", ref.substr(first, second-first+1).c_str());
            if(children.size()&&(ref[first+1]!='('||ref[second-1]!=')'))
                subrefs.push_back({ref.substr(first, second-first+1), seq.substr(first, second-first+1)});
        }
        else{
            printf("%s\n", seq.c_str());
            printf("%s\n", ref.c_str());
            // if(children.size())
            //     subrefs.push_back({ref, seq});
        }
        for(int i = 0; i < children.size(); i++){
            printf("child[%d]: first: %d, second: %d\n", i, children[i]->first, children[i]->second);
        }
        printf("\n");
        for(int j = 0; j < children.size(); j++){
            children[j]->printTree(ref, seq, subrefs);
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

bool compareByFirstStringLength(const std::pair<std::string, std::string> &a, const std::pair<std::string, std::string> &b) {
    return a.first.length() < b.first.length();
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

vector<int> ref2pairs(std::string& ref){
    vector<int> pairs(ref.length(), -1);
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
    vector<int> pairs_all = ref2pairs(ref);
    std::vector<std::tuple<int, int>> pairs_diff;
    for(auto& idx: positions){
        if(pairs_all[idx]==idx)
            pairs_diff.push_back(std::make_tuple(idx, idx));
        else if (pairs_all[idx]>idx)
            pairs_diff.push_back(std::make_tuple(idx, pairs_all[idx]));
    }
    return pairs_diff;
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

void print(Constraint& cs){
    for(int idx: *cs.indices)
        std::cout<<idx<<"\t";
    std::cout<<std::endl;
    for(std::string str: *cs.seqs){
        std::string substr = getSubstrings(*cs.indices, str);
        for(char c: substr)
            std::cout<<c<<"\t";
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
}

void intersect(Constraint* cs1, Constraint* cs2){
    std::set<int> idx_inter = setIntersection(*cs1->indices, *cs2->indices);
    if(idx_inter.empty())
        return;
    std::set<std::string> substr_cs_1;
    std::set<std::string> substr_cs_2;
    for(auto str: *cs1->seqs){
        substr_cs_1.insert(getSubstrings(idx_inter, str));
    }
    for(auto str: *cs2->seqs){
        substr_cs_2.insert(getSubstrings(idx_inter, str));
    }
    std::set<std::string> substr_inter = setIntersection(substr_cs_1, substr_cs_2);
    std::vector<std::string> seqs_new_1;
    std::vector<std::string> seqs_new_2;
    for (int i = 0; i < cs1->seqs->size(); i++){
        if (substr_inter.count(getSubstrings(idx_inter, (*cs1->seqs)[i])))
            seqs_new_1.push_back( (*cs1->seqs)[i]);
    }
    cs1->seqs = &seqs_new_1;

    for (int j = 0; j < cs2->seqs->size(); j++){
        if (substr_inter.count(getSubstrings(idx_inter, (*cs2->seqs)[j])))
            seqs_new_2.push_back( (*cs2->seqs)[j]);
    }
    cs2->seqs = &seqs_new_2;
}

void intersect(Constraint& cs1, Constraint& cs2){
    std::set<int> idx_inter = setIntersection(*cs1.indices, *cs2.indices);
    if(idx_inter.empty())
        return;
    std::set<std::string> substr_cs_1;
    std::set<std::string> substr_cs_2;
    for(auto str: *cs1.seqs){
        substr_cs_1.insert(getSubstrings(idx_inter, str));
    }
    for(auto str: *cs2.seqs){
        substr_cs_2.insert(getSubstrings(idx_inter, str));
    }
    std::set<std::string> substr_inter = setIntersection(substr_cs_1, substr_cs_2);
    std::vector<std::string> seqs_new_1;
    std::vector<std::string> seqs_new_2;
    for (int i = 0; i < cs1.seqs->size(); i++){
        if (substr_inter.count(getSubstrings(idx_inter, (*cs1.seqs)[i])))
            seqs_new_1.push_back( (*cs1.seqs)[i]);
    }
    *cs1.seqs = seqs_new_1;

    for (int j = 0; j < cs2.seqs->size(); j++){
        if (substr_inter.count(getSubstrings(idx_inter, (*cs2.seqs)[j])))
            seqs_new_2.push_back( (*cs2.seqs)[j]);
    }
    *cs2.seqs = seqs_new_2;
}

std::set<int> loops2positions(std::vector<std::vector<int>>& cr_loops, int length){
    std::set<int> critical_positions;
    for(vector<int> cr_loop: cr_loops){
        for(int i = 2; i < cr_loop.size(); i++)
            if (cr_loop[i]>=0 && cr_loop[i] < length)
                critical_positions.insert(cr_loop[i]);
    }
    return critical_positions;
}

std::vector<std::pair<std::string, std::set<int>>> genSeqCritical(std::string seq,  std::set<int> critical_positions, std::vector<std::pair<int, int>> sp_y){
    std::vector<std::pair<std::string, std::set<int>>> seqCritical;
    std::vector<std::string> first_sp;
    std::vector<std::string> second_sp;
    int L1 = sp_y[0].second - sp_y[0].first - 1;
    int L2 = sp_y[1].second - sp_y[1].first - 1;
    assert (L1 == 3 || L1 == 4);
    assert (L2 == 3 || L2 == 4);
    if ( L1 == 3)
        first_sp = TriHP;
    else{
        first_sp = TetraHP;
    }
    if (L2 == 3)
        second_sp = TriHP;
    else{
        second_sp = TetraHP;
    }
    for(auto sp1: first_sp)
        for(auto sp2: second_sp){
            std::string seq_new = seq;
            std::set<int> critical_positions_new = critical_positions;
            if(sp1 == ""){
                if(L1 == 3){
                    seq_new[sp_y[0].first+2] = 'C';
                    critical_positions_new.erase(sp_y[0].first+2);
                }else{
                    seq_new[sp_y[0].first+2] = 'C';
                    seq_new[sp_y[0].first+3] = 'G';
                    critical_positions_new.erase(sp_y[0].first+2);
                    critical_positions_new.erase(sp_y[0].first+3);
                }
            }else{
                if(L1 == 3){
                    seq_new.replace(sp_y[0].first, sp1.length(), sp1);
                    for(int i = 0; i < 5; i++)
                        critical_positions_new.erase(sp_y[0].first+i);
                }else{
                    seq_new.replace(sp_y[0].first, sp1.length(), sp1);
                    for(int i = 0; i < 6; i++)
                        critical_positions_new.erase(sp_y[0].first+i);
                }
            }

            if(sp2 == ""){
                if(L2 == 3){
                    seq_new[sp_y[1].first+2] = 'C';
                    critical_positions_new.erase(sp_y[1].first+2);
                }else{
                    seq_new[sp_y[1].first+2] = 'C';
                    seq_new[sp_y[1].first+3] = 'G';
                    critical_positions_new.erase(sp_y[1].first+2);
                    critical_positions_new.erase(sp_y[1].first+3);
                }
            }else{
                if(L2 == 3){
                    seq_new.replace(sp_y[1].first, sp2.length(), sp2);
                    for(int i = 0; i < 5; i++)
                        critical_positions_new.erase(sp_y[1].first+i);
                }else{
                    seq_new.replace(sp_y[1].first, sp2.length(), sp2);
                    for(int i = 0; i < 6; i++)
                        critical_positions_new.erase(sp_y[1].first+i);
                }
            }

            seqCritical.push_back({seq_new, critical_positions_new});
        }
    return seqCritical;
}

std::vector<std::pair<int, int>> loops2specialhp(std::vector<std::vector<int>>& cr_loops, int length){
    std::vector<std::pair<int, int>> sp_y;
    std::set<int> critical_positions_1; // y
    std::set<int> critical_positions_0; // y'
    for(vector<int> cr_loop: cr_loops){
        if (cr_loop[0] == 0){  // y'
            for(int i = 2; i < cr_loop.size(); i++){
                if (cr_loop[i]>=0 && cr_loop[i] < length)
                    critical_positions_0.insert(cr_loop[i]);
            }
        }
    }
    for(vector<int> cr_loop: cr_loops){
        if (cr_loop[0] == 1 && cr_loop[1] == 0){ // sp in crloop of y
            int lenHP = cr_loop[3] - cr_loop[2] - 1;
            #ifdef SPECIAL_HP_3
            if(lenHP==3){
                set<int> indices_hp;
                for(int i = 4; i < cr_loop.size(); i++){
                    indices_hp.insert(cr_loop[i]);
                }
                std::set<int> intersection = setIntersection(indices_hp,critical_positions_0);
                if (intersection.empty())
                {
                    sp_y.push_back({cr_loop[2], cr_loop[3]});
                }
            }
            #endif
            #ifdef SPECIAL_HP_4
            if(lenHP==4){
                set<int> indices_hp;
                for(int i = 4; i < cr_loop.size(); i++){
                    indices_hp.insert(cr_loop[i]);
                }
                std::set<int> intersection = setIntersection(indices_hp,critical_positions_0);
                if (intersection.empty())
                {
                    sp_y.push_back({cr_loop[2], cr_loop[3]});
                }
            }
            #endif
        }
    }
    return sp_y;
}

std::vector<std::string> alg_1(std::string& y, std::string& y_prime, std::vector<std::vector<int>>& cr_loops, std::vector<std::tuple<int, int>>& pairs_diff, std::string& seq, bool is_verbose, int dangle_model){
    std::cout<<"inside alg1"<<std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    ulong nEnum = count_enum(pairs_diff);
    std::cout<<"count enum: "<<nEnum<<std::endl;
    std::vector<std::pair<int, std::string>> idX;
    std::vector<std::string> X;
    int flag = 0;
    #pragma omp parallel for
    for(ulong i=0; i < nEnum; i++){
        if(flag)
            continue;
        std::string seq_i = enumerate(pairs_diff, i, seq);
        if ((i+1)%1000000 == 0){
            auto pause = std::chrono::high_resolution_clock::now();
            printf("%8d, %d, %.2f seconds\n", (i+1)/10000, X.size(), std::chrono::duration<double, std::milli>(pause - start)/1000.f);
        }
        if(check_compatible(seq_i, y_prime)){
            long e_diff = -diff_eval(seq_i, cr_loops, is_verbose, dangle_model); // not divided by 100
            if(e_diff < 0){
                #pragma omp critical
                idX.push_back({i, seq_i});
                // X.push_back(seq_i);
            }
        }else{
            #pragma omp critical
            idX.push_back({i, seq_i});
            // X.push_back(seq_i);
        }
        #pragma omp critical
        {
            if (idX.size() > MAX_CONSTRAINT)
                flag = 1;
        }
    }
    auto stop = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double, std::milli> fp_ms = stop - start;
    printf("finished: %.4f seconds\n", fp_ms/1000.f);
    std::sort(idX.begin(), idX.end());
    for(auto p: idX)
        X.push_back(p.second);
    return X;
}

std::vector<std::string> alg_1_v2(std::string& y, std::string& y_prime, std::string& seq, bool is_verbose, int dangle_model){
    std::cout<<"inside alg1_v2"<<std::endl;
    std::set<int> critical_positions;
    std::vector<std::vector<int>> cr_loops = find_critical_plus(y, y_prime, critical_positions, is_verbose);
    std::set<int> critical_positions_v2 = loops2positions(cr_loops, y.length());
    assert(critical_positions == critical_positions_v2);
    
    std::vector<std::pair<int, int>> sh_y = loops2specialhp(cr_loops, y.length());
    if (sh_y.size() == 2){
        std::vector<std::string> X;
        std::vector<std::pair<int, std::string>> idX;
        auto seqCriticals = genSeqCritical(seq, critical_positions, sh_y);
        auto start = std::chrono::high_resolution_clock::now();
        for(int j = 0; j < seqCriticals.size(); j++){
            std::cout<<"SH:"<<j<<"\t"<<seqCriticals[j].first<<"\t"<<seqCriticals[j].second.size()<<std::endl;
            std::cout<<"SH:"<<j<<"\t"<<y<<std::endl;
            std::vector<std::string> X;
            std::vector<std::tuple<int, int>> pairs_diff = idx2pair(seqCriticals[j].second, y);
            if (is_verbose)
                for(auto& pair: pairs_diff)
                    std::cout<<std::get<0>(pair)<<"\t"<<std::get<1>(pair)<<std::endl;
            auto start = std::chrono::high_resolution_clock::now();
            ulong nEnum = count_enum(pairs_diff);
            std::cout<<"count enum: "<<nEnum<<std::endl;
            int flag = 0;
            #pragma omp parallel for
            for(ulong i=0; i < nEnum; i++){
                if(flag)
                    continue;
                std::string seq_i = enumerate(pairs_diff, i, seqCriticals[j].first);
                if ((i+1)%1000000 == 0){
                    auto pause = std::chrono::high_resolution_clock::now();
                    printf("%8d, %d, %.2f seconds\n", (i+1)/10000, X.size(), std::chrono::duration<double, std::milli>(pause - start)/1000.f);
                }
                if(check_compatible(seq_i, y_prime)){
                    long e_diff = -diff_eval(seq_i, cr_loops, is_verbose, dangle_model); // not divided by 100
                    if(e_diff < 0){
                        #pragma omp critical
                        idX.push_back({i, seq_i});
                        // X.push_back(seq_i);
                    }
                }else{
                    #pragma omp critical
                    idX.push_back({i, seq_i});
                    // X.push_back(seq_i);
                }
                #pragma omp critical
                {
                    if (idX.size() > MAX_CONSTRAINT)
                        flag = 1;
                }
            }
        }
        auto stop = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<double, std::milli> fp_ms = stop - start;
        printf("finished: %.4f seconds\n", fp_ms/1000.f);
        std::sort(idX.begin(), idX.end());
        for(auto p: idX)
        X.push_back(p.second);
        return X;
    }
    std::vector<std::string> X;
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<std::tuple<int, int>> pairs_diff = idx2pair(critical_positions, y);
    ulong nEnum = count_enum(pairs_diff);
    std::cout<<"count enum: "<<nEnum<<std::endl;
    std::vector<std::pair<int, std::string>> idX;

    // if(nEnum <= 0 || nEnum >= MAX_ENUM){
    //     printf("too many enumeration needed!");
    //     return X;
    // }

    int flag = 0;
    #pragma omp parallel for
    for(ulong i=0; i < nEnum; i++){
        if(flag)
            continue;
        std::string seq_i = enumerate(pairs_diff, i, seq);
        if ((i+1)%1000000 == 0){
            auto pause = std::chrono::high_resolution_clock::now();
            printf("%8d, %d, %.2f seconds\n", (i+1)/10000, X.size(), std::chrono::duration<double, std::milli>(pause - start)/1000.f);
        }
        if(check_compatible(seq_i, y_prime)){
            long e_diff = -diff_eval(seq_i, cr_loops, is_verbose, dangle_model); // not divided by 100
            if(e_diff < 0){
                #pragma omp critical
                idX.push_back({i, seq_i});
                // X.push_back(seq_i);
            }
        }else{
            #pragma omp critical
            idX.push_back({i, seq_i});
            // X.push_back(seq_i);
        }
        #pragma omp critical
        {
            if (idX.size() > MAX_CONSTRAINT)
                flag = 1;
        }
    }
    auto stop = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double, std::milli> fp_ms = stop - start;
    printf("finished: %.4f seconds\n", fp_ms/1000.f);
    std::sort(idX.begin(), idX.end());
    for(auto p: idX)
        X.push_back(p.second);
    return X;
}

std::string alg_2(std::string& ref1, std::set<std::string>& refs_checked, std::vector<Constraint>& cs_vec, bool verbose, int dangle_model){ // ref1, ref2, X, is_verbose, dangle_model
    int count_cs = cs_vec.size();
    std::vector<std::pair<ulong, std::pair<std::string, std::string>>> y_primes;
    std::cout<<"inside alg2"<<std::endl;
    std::vector<std::string> X;
    for(auto cs: cs_vec){
        X.insert(X.end(), cs.seqs->begin(), cs.seqs->end());
        if(X.size() > 500)
            break;
    }
    std::cout<<"X.size: "<<X.size()<<std::endl;
    for(auto x: X){
        assert (check_compatible(x, ref1));
        std::vector<std::string> subopts = fold(x, 0, false, verbose, dangle_model, 0.0);
        if (isUMFE(subopts, ref1)){
            printf("UMFE found!");
            std::cout<<x<<std::endl;
            return "UMFE";
        }
        for(std::string ref: subopts){
            if(ref != ref1 && !refs_checked.count(ref)){
                std::set<int> critical_positions;
                find_critical_plus(ref1, ref, critical_positions, verbose);
                std::vector<std::tuple<int, int>> pairs_diff = idx2pair(critical_positions, ref1);
                ulong n_enum = count_enum(pairs_diff);
                std::pair<ulong, std::pair<std::string, std::string>> enum_seq(n_enum, std::make_pair(ref, x));
                if (n_enum > 0 && n_enum < MAX_ENUM)
                    y_primes.push_back(enum_seq);
                else{
                    std::cout<<"seq: "<<x<<std::endl;
                    std::cout<<"ref: "<<ref<<std::endl;
                    std::cout<<"enum:"<<n_enum<<std::endl;
                }
            }
        }
    }
    std::sort(y_primes.begin(), y_primes.end());
    std::vector<std::pair<ulong, std::pair<std::string, std::string>>> y_primes_dedup;
    for(auto y_prime: y_primes){
        if (y_primes_dedup.empty() || y_prime.second.first != y_primes_dedup.back().second.first){
            std::cout<<y_prime.first<<std::endl;
            std::cout<<y_prime.second.second<<std::endl;
            std::cout<<y_prime.second.first<<std::endl;
            y_primes_dedup.push_back(y_prime);
        }
    }
    for (auto y_prime: y_primes_dedup){
        std::set<int> critical_positions;
        std::vector<std::vector<int>> cr_loops = find_critical_plus(ref1, y_prime.second.first, critical_positions, verbose);
        std::vector<std::tuple<int, int>> pairs_diff = idx2pair(critical_positions, ref1);
        // std::vector<std::string> X_new = alg_1(ref1, y_prime.second.first, cr_loops, pairs_diff, y_prime.second.second, verbose, dangle_model);
        std::vector<std::string> X_new = alg_1_v2(ref1, y_prime.second.first, y_prime.second.second, verbose, dangle_model);
        if (X_new.size() == 0){
            std::cout<<"y :"<<ref1<<std::endl;
            std::cout<<"y':"<<y_prime.second.first<<std::endl;
            std::cout<<"undesignable!"<<std::endl;
            return "undesignable";
        }else if (X_new.size() > MAX_CONSTRAINT)
            std::cout<<"too many constraints: "<<X_new.size()<<"\t"<<"out of "<<y_prime.first<<std::endl;
        else{
            Constraint cs_new = Constraint(&critical_positions, &X_new);
            cs_new.setStructure(y_prime.second.first);
            std::cout<<"number of constraints: "<<X_new.size()<<std::endl;
            for(Constraint& cs_old: cs_vec){
                // std::cout<<"before intersection: "<<cs_old.seqs->size()<<"\t"<<cs_new.seqs->size()<<std::endl;
                intersect(cs_old, cs_new);
                // std::cout<<"after  intersection: "<<cs_old.seqs->size()<<"\t"<<cs_new.seqs->size()<<std::endl;
                if (cs_old.seqs->empty() || cs_new.seqs->empty()){
                    for(int i = 0; i<cs_vec.size(); i++)
                        std::cout<<cs_vec[i].seqs->size()<<"\t";
                    std::cout<<std::endl;
                    for(int i = 0; i<cs_vec.size(); i++)
                        std::cout<<cs_vec[i].structure<<std::endl;
                    std::cout<<cs_new.structure<<std::endl;
                    std::cout<<"y_prime count: "<<cs_vec.size()+1<<std::endl;
                    std::cout<<"undesignable!"<<std::endl;
                    return "undesignable";
                }
            }
            std::set<int>* idx_new = new std::set<int>(*cs_new.indices);
            std::vector<std::string>* x_new_copy = new std::vector<std::string>(*cs_new.seqs);
            Constraint* cs_new_copy = new Constraint(idx_new, x_new_copy);
            cs_new_copy->setStructure(cs_new.structure);
            cs_vec.push_back(*cs_new_copy);
            for(int i = 0; i<cs_vec.size(); i++)
                std::cout<<cs_vec[i].seqs->size()<<"\t";
            std::cout<<cs_new.seqs->size()<<std::endl;
        }
        refs_checked.insert(y_prime.second.first);
    }
    if (count_cs == cs_vec.size()){
        std::cout<<"no more new y_prime"<<std::endl;
        for(auto cs: cs_vec){
            print(cs);
            std::cout<<std::endl;
        }
        return "no more new y_prime";
    }
    if (cs_vec.size() < 100)
        return alg_2(ref1, refs_checked, cs_vec, verbose, dangle_model);
    else{
        std::cout<<"no conclusion!"<<std::endl;
        return "no conclusion";
    }
}

std::string alg_2_cs(std::string& ref1, std::set<std::string>& refs_checked, std::vector<Constraint>& cs_vec, bool verbose, int dangle_model){ // ref1, ref2, X, is_verbose, dangle_model
    int count_cs = cs_vec.size();
    std::vector<std::pair<ulong, std::pair<std::string, std::string>>> y_primes;
    std::cout<<"inside alg2cs"<<std::endl;
    std::vector<std::string> X;
    for(auto cs: cs_vec){
        X.insert(X.end(), cs.seqs->begin(), cs.seqs->end());
        if(X.size() > 500)
            break;
    }
    if (X.size() > 500)
        X.resize(500);
    std::cout<<"X.size: "<<X.size()<<std::endl;
    std::string constr(ref1.length(), '?');
    constr[0] = '(';
    constr[ref1.length()-1] = ')';
    for(auto x: X){
        assert (check_compatible(x, ref1));
        std::vector<std::string> subopts_raw = cs_fold(x, constr, 0, false, verbose, dangle_model);
        std::vector<std::string> subopts;
        for(std::string subopt: subopts_raw){
            if(subopt[ref1.length()-1]==')' && subopt[0]=='(')
                subopts.push_back(subopt);
        }
        if (isUMFE(subopts, ref1)){
            printf("UMFE found!");
            std::cout<<x<<std::endl;
            return "UMFE";
        }
        for(std::string ref: subopts){
            if(ref != ref1 && !refs_checked.count(ref)){
                std::set<int> critical_positions;
                find_critical_plus(ref1, ref, critical_positions, verbose);
                std::vector<std::tuple<int, int>> pairs_diff = idx2pair(critical_positions, ref1);
                ulong n_enum = count_enum(pairs_diff);
                std::pair<ulong, std::pair<std::string, std::string>> enum_seq(n_enum, std::make_pair(ref, x));
                if (n_enum > 0 && n_enum < MAX_ENUM)
                    y_primes.push_back(enum_seq);
            }
        }
    }
    std::sort(y_primes.begin(), y_primes.end());
    std::vector<std::pair<ulong, std::pair<std::string, std::string>>> y_primes_dedup;
    for(auto y_prime: y_primes){
        if (y_primes_dedup.empty() || y_prime.second.first != y_primes_dedup.back().second.first){
            std::cout<<y_prime.first<<std::endl;
            std::cout<<y_prime.second.second<<std::endl;
            std::cout<<y_prime.second.first<<std::endl;
            y_primes_dedup.push_back(y_prime);
        }
    }
    for (auto y_prime: y_primes_dedup){
        std::set<int> critical_positions;
        std::vector<std::vector<int>> cr_loops = find_critical_plus(ref1, y_prime.second.first, critical_positions, verbose);
        std::vector<std::tuple<int, int>> pairs_diff = idx2pair(critical_positions, ref1);
        // std::vector<std::string> X_new = alg_1(ref1, y_prime.second.first, cr_loops, pairs_diff, y_prime.second.second, verbose, dangle_model);
        std::vector<std::string> X_new = alg_1_v2(ref1, y_prime.second.first, y_prime.second.second, verbose, dangle_model);
        if (X_new.size() == 0){
            std::cout<<"y :"<<ref1<<std::endl;
            std::cout<<"y':"<<y_prime.second.first<<std::endl;
            std::cout<<"undesignable!"<<std::endl;
            return "undesignable";
        }else if (X_new.size() > MAX_CONSTRAINT){
            std::cout<<"too many constraints: "<<X_new.size()<<"\t"<<"out of "<<y_prime.first<<std::endl;
            refs_checked.insert(y_prime.second.first);
            continue;
        }else{
            Constraint cs_new = Constraint(&critical_positions, &X_new);
            cs_new.setStructure(y_prime.second.first);
            std::cout<<"number of constraints: "<<X_new.size()<<std::endl;
            for(Constraint& cs_old: cs_vec){
                intersect(cs_old, cs_new);
                if (cs_old.seqs->empty() || cs_new.seqs->empty()){
                    for(int i = 0; i<cs_vec.size(); i++)
                        std::cout<<cs_vec[i].seqs->size()<<"\t";
                    std::cout<<cs_new.seqs->size()<<std::endl;
                    for(int i = 0; i<cs_vec.size(); i++)
                        std::cout<<cs_vec[i].structure<<std::endl;
                    std::cout<<cs_new.structure<<std::endl;
                    std::cout<<"y_prime count: "<<cs_vec.size()+1<<std::endl;
                    std::cout<<"undesignable!"<<std::endl;
                    return "undesignable";
                }
            }
            std::set<int>* idx_new = new std::set<int>(*cs_new.indices);
            std::vector<std::string>* x_new_copy = new std::vector<std::string>(*cs_new.seqs);
            Constraint* cs_new_copy = new Constraint(idx_new, x_new_copy);
            cs_new_copy->setStructure(cs_new.structure);
            cs_vec.push_back(*cs_new_copy);
            for(int i = 0; i<cs_vec.size(); i++)
                std::cout<<cs_vec[i].seqs->size()<<"\t";
            std::cout<<std::endl;
        }
        refs_checked.insert(y_prime.second.first);
    }
    if (count_cs == cs_vec.size()){
        std::cout<<"no more new y_prime"<<std::endl;
        for(auto cs: cs_vec){
            std::cout<<cs.seqs->size()<<std::endl;
            std::cout<<std::endl;
        }
        return "no more new y_prime";
    }
    if (cs_vec.size() < 100)
        return alg_2_cs(ref1, refs_checked, cs_vec, verbose, dangle_model);
    else{
        std::cout<<"no conclusion!"<<std::endl;
        return "no conclusion";
    }
}

std::string alg_2_helper(std::string& ref1, std::string& ref2, std::string& seq, bool verbose, int dangle_model){
    std::cout << "seq: " << seq << std::endl;
    std::cout << "  y: " << ref1 << std::endl;
    std::cout << " y': " << ref2 << std::endl;
    std::set<int> critical_positions;
    std::vector<std::vector<int>> cr_loops = find_critical_plus(ref1, ref2, critical_positions, verbose);
    long delta_energy = diff_eval(seq, cr_loops, verbose, dangle_model);
    std::vector<std::tuple<int, int>> pairs_diff = idx2pair(critical_positions, ref1);
    ulong n_enum = count_enum(pairs_diff);
    std::cout<<"enumeration count: "<<n_enum<<std::endl;
    if(n_enum > 0 && n_enum < MAX_ENUM){
        std::cout<<"alg 1"<<std::endl;
        // auto X = alg_1(ref1, ref2, cr_loops, pairs_diff, seq, verbose, dangle_model);
        auto X = alg_1_v2(ref1, ref2, seq, verbose, dangle_model);
        printf("X size: %d\n", X.size());
        if (X.size() == 0){
            std::cout<<"undesignable!"<<std::endl;
            return "undesignable";
        }
        // X.push_back("AAAACGGGAACGCUCAACCCGAAGGCCAAAAAGGCCCCGCGACAAUCGACGGCGGGGCGGGGACGAGGAGCGCCAAAAGGACGCCCCGGGGCCAGCACGAGCAAAAGCCGGCCGCCACCGAAAACGAGGAGCGAAGGACCCCCCACGCGAGCGCCGAGCGAGGGAGGGCAAAAGCCCCCCCGCCCGGCGGCAAAGCCAGCACGGACGGCCCG");
        // X.push_back("AAAAGGCCAACGCCCAACCCCAACGCGAAAAACGCGGGCGGUAACACCAGCCGGGGCCCGGCAGGACGACGGGGCAAGCCACCGCGCCCGGCCACCAGCAGCAAAAGCGCGGGCGCACGCAAAAGCAGGACGCAAGGACCGGCCAGGGGGGCGGCGAGGCAGCCAGCGCAAAAGCGCGGCGCCGCCGCCGGAAACCGAGCAGGCACGCGGCC");
        // X.push_back("AACACGGGAACGCUCAACCCGAAGGCCAAAAAGGCCCCGCGACAAUCGACGGCGGGGCGGGGACGAGGAGCGCCAAAAGGACGCCCCGGGGCCAGCACGAGCAAAAGCCGGCCGCCAGCGAAAACGAUGAGCGAAGGACCCCCCACGCGAGCGCCGAGCGAUGGAGGGCAAAAGCCCCUACGCCCGGCGGCAAAGCCAGCACGGACGGCCCG");
        // X.push_back("AAAACUGGAACGCUCAACCCGAAGGCCAUAAAGGCCCCGCGACAAUCGACGGCGGGGCGGGGAGUAGGAGCGCCAAAAGGACGCCCACGGGCCAGCACGAGCAAAAGCCGGCCGCCAUCGAAAACGAAGAGCGAAGGACCCCCCACGCGAGCGCCGAGCGAUGGAGGGCAAAAGCCCCUACGCCCAGCGGCAAAGCCAGCAUGGACGGCCAG");
        // X.push_back("AGAACGGGAACGCUCAACCCGAAGGCCAAAAAGGCCCCCCGACAAUCGAGGGCGGGGCGGGGACGAGGAGCGCCAAAAGGACGCCCCGGGGCCAGCACGAGGAAAUCCCGGCCGCUAGCGAAAACGACGAGCGAAGGACCCCCCACGCGAGCGCCGAGCGAUGCAGGGCAAAAGCCCGUACGCCCAGCCCCAAAGGGAGCAUGGACGGCCCG");
        // std::vector<std::string> X_samfeo = readLinesFromFile("design.txt");
        // X.insert(X.end(), X_samfeo.begin(), X_samfeo.end());
        std::set<std::string> refs_checked;
        refs_checked.insert(ref2);
        std::vector<Constraint> cs_vec;
        Constraint cs_ref2 = Constraint(&critical_positions, &X);
        cs_ref2.setStructure(ref2);
        cs_vec.push_back(cs_ref2);
        return alg_2(ref1, refs_checked, cs_vec, verbose, dangle_model);
    }
    std::cout<<"intial y_prime too bad!"<<std::endl;
    return "intial y_prime too bad";
}

std::string alg_2_cs_helper(std::string& ref1, std::string& ref2, std::string& seq, bool verbose, int dangle_model){
    std::cout << "seq: " << seq << std::endl;
    std::cout << "  y: " << ref1 << std::endl;
    std::cout << " y': " << ref2 << std::endl;
    std::set<int> critical_positions;
    std::vector<std::vector<int>> cr_loops = find_critical_plus(ref1, ref2, critical_positions, verbose);
    long delta_energy = diff_eval(seq, cr_loops, verbose, dangle_model);
    std::vector<std::tuple<int, int>> pairs_diff = idx2pair(critical_positions, ref1);
    ulong n_enum = count_enum(pairs_diff);
    std::cout<<"enumeration count: "<<n_enum<<std::endl;
    if(n_enum > 0 && n_enum < MAX_ENUM){
        std::cout<<"alg 1"<<std::endl;
        // auto X = alg_1(ref1, ref2, cr_loops, pairs_diff, seq, verbose, dangle_model);
        auto X = alg_1_v2(ref1, ref2, seq, verbose, dangle_model);
        printf("X size: %d\n", X.size());
        if (X.size() == 0){
            std::cout<<"undesignable!"<<std::endl;
            return "undesignable";
        }
        // X.push_back("AAAACGGGAACGCUCAACCCGAAGGCCAAAAAGGCCCCGCGACAAUCGACGGCGGGGCGGGGACGAGGAGCGCCAAAAGGACGCCCCGGGGCCAGCACGAGCAAAAGCCGGCCGCCACCGAAAACGAGGAGCGAAGGACCCCCCACGCGAGCGCCGAGCGAGGGAGGGCAAAAGCCCCCCCGCCCGGCGGCAAAGCCAGCACGGACGGCCCG");
        // X.push_back("AAAAGGCCAACGCCCAACCCCAACGCGAAAAACGCGGGCGGUAACACCAGCCGGGGCCCGGCAGGACGACGGGGCAAGCCACCGCGCCCGGCCACCAGCAGCAAAAGCGCGGGCGCACGCAAAAGCAGGACGCAAGGACCGGCCAGGGGGGCGGCGAGGCAGCCAGCGCAAAAGCGCGGCGCCGCCGCCGGAAACCGAGCAGGCACGCGGCC");
        // X.push_back("AACACGGGAACGCUCAACCCGAAGGCCAAAAAGGCCCCGCGACAAUCGACGGCGGGGCGGGGACGAGGAGCGCCAAAAGGACGCCCCGGGGCCAGCACGAGCAAAAGCCGGCCGCCAGCGAAAACGAUGAGCGAAGGACCCCCCACGCGAGCGCCGAGCGAUGGAGGGCAAAAGCCCCUACGCCCGGCGGCAAAGCCAGCACGGACGGCCCG");
        // X.push_back("AAAACUGGAACGCUCAACCCGAAGGCCAUAAAGGCCCCGCGACAAUCGACGGCGGGGCGGGGAGUAGGAGCGCCAAAAGGACGCCCACGGGCCAGCACGAGCAAAAGCCGGCCGCCAUCGAAAACGAAGAGCGAAGGACCCCCCACGCGAGCGCCGAGCGAUGGAGGGCAAAAGCCCCUACGCCCAGCGGCAAAGCCAGCAUGGACGGCCAG");
        // X.push_back("AGAACGGGAACGCUCAACCCGAAGGCCAAAAAGGCCCCCCGACAAUCGAGGGCGGGGCGGGGACGAGGAGCGCCAAAAGGACGCCCCGGGGCCAGCACGAGGAAAUCCCGGCCGCUAGCGAAAACGACGAGCGAAGGACCCCCCACGCGAGCGCCGAGCGAUGCAGGGCAAAAGCCCGUACGCCCAGCCCCAAAGGGAGCAUGGACGGCCCG");
        // std::vector<std::string> X_samfeo = readLinesFromFile("design.txt");
        // X.insert(X.end(), X_samfeo.begin(), X_samfeo.end());
        std::set<std::string> refs_checked;
        refs_checked.insert(ref2);
        std::vector<Constraint> cs_vec;
        Constraint cs_ref2 = Constraint(&critical_positions, &X);
        cs_ref2.setStructure(ref2);
        cs_vec.push_back(cs_ref2);
        return alg_2_cs(ref1, refs_checked, cs_vec, verbose, dangle_model);
    }
    std::cout<<"intial y_prime too bad!"<<std::endl;
    return "intial y_prime too bad";
}

void alg_3_helper(std::string& ref, std::string& seq, bool verbose, int dangle){
    TreeNode* root = parseStringToTree(ref);
    std::vector<std::pair<std::string, std::string>> subrefs;
    root->printTree(ref, seq, subrefs);
    std::cout<<"size of sub refs: "<<subrefs.size()<<std::endl;
    std::sort(subrefs.begin(), subrefs.end(), compareByFirstStringLength);
    for(int i = 0; i < subrefs.size(); i++){
        printf("%2d: L=%d\n", i, subrefs[i].first.length());
        std::cout<<subrefs[i].second<<std::endl;
        std::cout<<subrefs[i].first<<std::endl;

        std::string constr( subrefs[i].first.length(), '?');
        constr[0] = '(';
        constr[subrefs[i].first.length()-1] = ')';
        std::vector<std::string> subopts_raw = cs_fold(subrefs[i].second, constr, 0, false, false, dangle);
        std::vector<std::string> subopts;
        for(std::string subopt: subopts_raw){
            if(subopt[subrefs[i].first.length()-1]==')' && subopt[0]=='(')
                subopts.push_back(subopt);
        }
        if (isUMFE(subopts,  subrefs[i].first)){
            printf("(local) UMFE found!");
            std::cout<<subrefs[i].second<<std::endl;
            std::cout<<subrefs[i].first<<std::endl;
        }else{
            std::string ref_mfe = subopts[0];
            if(subopts[0]==subrefs[i].first)
                ref_mfe = subopts[1];
            std::string designability = alg_2_cs_helper(subrefs[i].first, ref_mfe, subrefs[i].second, verbose, dangle);
            if (designability == "undesignable")
                break;
        }
    }
    return;
}

std::vector<std::string> cs_fold(std::string seq, std::string& constr, int beamsize, bool sharpturn, bool verbose, int dangle){
    bool consflag = true;
    std::vector<std::string> subopts;
    std::set<char> consSet {'?', '.', '(', ')'};
    if (seq.length() != constr.length()){
        printf("The lengths don't match between sequence and constraints: %s, %s\n", seq.c_str(), constr.c_str());
        return subopts;
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
        BeamCKYParser parser(beamsize, !sharpturn, verbose, true, true, 0.0, "", false, dangle);
        
        BeamCKYParser::DecoderResult result = parser.parse(seq, &cons, subopts);
        // BeamCKYParser::DecoderResult result = parser.parse(seq, &cons);

        #ifdef lv
                double printscore = (result.score / -100.0);
        #else
                double printscore = result.score;
        #endif
        // printf("%s (%.2f)\n", result.structure.c_str(), printscore);
        return subopts;
    }
    return subopts;
}


std::vector<std::string> fold(std::string seq, int beamsize, bool sharpturn, bool verbose, int dangle, float energy_delta = 0.){
        // lhuang: moved inside loop, fixing an obscure but crucial bug in initialization
        BeamCKYParser parser(beamsize, !sharpturn, verbose, false, true, energy_delta, "", false, dangle);
        std::vector<std::string> subopts;
        BeamCKYParser::DecoderResult result = parser.parse(seq, NULL, subopts);
        return subopts;
}


std::vector<std::string> alg1_helper(std::string& seq, std::string& ref1, std::string& ref2, bool is_verbose, int dangle_model) {
    std::cout << "seq: " << seq << std::endl;
    std::cout << "  y: " << ref1 << std::endl;
    std::cout << " y': " << ref2 << std::endl;
    std::set<int> critical_positions;
    std::vector<std::vector<int>> cr_loops = find_critical_plus(ref1, ref2, critical_positions, is_verbose);
    std::set<int> critical_positions_v2 = loops2positions(cr_loops, ref1.length());
    if(critical_positions != critical_positions_v2){
        std::cout<<critical_positions.size()<<std::endl;
        std::cout<<critical_positions_v2.size()<<std::endl;
        assert (critical_positions == critical_positions_v2);
    }
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
    std::vector<std::string> X;
    if(n_enum > 0 && n_enum < MAX_ENUM){
        std::cout<<"alg 1"<<std::endl;
        auto X = alg_1(ref1, ref2, cr_loops, pairs_diff, seq, is_verbose, dangle_model);
        printf("X size: %d\n", X.size());
        if (X.size()==0)
            printf("undesignable!\n");
    }
    return X;
}

void test_cs(std::string& seq, std::string& ref1, std::string& ref2, bool is_verbose, int dangle_model) {
    std::cout << "seq: " << seq << std::endl;
    std::cout << "  y: " << ref1 << std::endl;
    std::cout << " y': " << ref2 << std::endl;
    std::set<int> critical_positions;
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
    if(n_enum > 0 && n_enum < MAX_ENUM){
        std::cout<<"alg 1"<<std::endl;
        auto X = alg_1(ref1, ref2, cr_loops, pairs_diff, seq, is_verbose, dangle_model);
        printf("X size: %d\n", X.size());
        std::set<std::string> refs;
        std::string constr(ref1.length(), '?');
        constr[0] = '(';
        constr[ref1.length()-1] = ')';
        std::cout<<constr<<std::endl;
        for(int i; i < min<int>(1000, X.size()); i++){
            std::string ref_x = cs_fold(X[i], constr, 0, false, false, 2)[0];
            // std::cout<<X[i]<<std::endl;
            // std::cout<<ref_x<<std::endl;
            if(ref_x == ref1){
                std::cout<<"local solution: "<<X[i]<<std::endl;
                std::cout<<"constraints:    "<<constr<<std::endl;
                std::cout<<"ref 1:          "<<ref1<<std::endl;
                std::cout<<"constr fold:    "<<ref_x<<std::endl;
            }
            refs.insert(ref_x);
        }
        std::cout<<"refs size: "<<refs.size()<<std::endl;
        for(auto ref_new: refs){
            std::set<int> critical_positions_new;
            cr_loops = find_critical_plus(ref1, ref_new, critical_positions_new, is_verbose);
            pairs_diff = idx2pair(critical_positions_new, ref1);
            n_enum = count_enum(pairs_diff);
            std::cout<<ref_new<<"\t"<<n_enum<<std::endl;
        }
    }
}

void show_configuration(){
    #ifdef SPECIAL_HP
    printf("SPECIAL_HP   defined.\n");
    #endif
    #ifdef SPECIAL_HP_3
    printf("SPECIAL_HP_3 defined.\n");
    #endif
    #ifdef SPECIAL_HP_4
    printf("SPECIAL_HP_4 defined.\n");
    #endif
    #ifdef SPECIAL_HP_6
    printf("SPECIAL_HP_6 defined.\n");
    #endif
    return;
}

int main(int argc, char* argv[]) {
    show_configuration();
    char* alg = (argc < 2) ? nullptr : argv[1];
    std::cout << "alg: " << ((alg == nullptr) ? "None" : alg) << std::endl;
    if (alg == nullptr){
        std::cout<<"no alg was selected!"<<std::endl;
        return 0;
    }
    bool verbose = false;
    int dangle = 2;
    if ( alg != nullptr && strcmp(alg, "csfold") == 0 ){  /* constrained folding */
        std::cout << alg << std::endl;
        int beamsize = 0;
        bool sharpturn = false;
        std::string seq;
        std::string constr;

        while (std::getline(std::cin, seq))
        {
            std::getline(std::cin, constr);
            auto refs =  cs_fold(seq, constr, beamsize, sharpturn, verbose, dangle);
            std::cout<<"subopts size: "<<refs.size()<<std::endl;
            for(auto ref: refs)
                std::cout<<ref<<std::endl;
        }
        return 0;
    }else if (alg != nullptr && strcmp(alg, "fold") == 0){ /* fold */
        int beamsize = 0;
        bool sharpturn = false;
        std::string seq;
        while (std::getline(std::cin, seq))
        {
            std::vector<std::string> refs =  fold(seq, beamsize, sharpturn, verbose, dangle);
            std::cout<<"subopts size: "<<refs.size()<<std::endl;
            for(auto ref: refs)
                std::cout<<ref<<std::endl;
        }
        return 0;
    }else if (alg != nullptr && strcmp(alg, "alg3") == 0){ /* alg 3 */
        std::cout << alg << std::endl;
        int beamsize = 0;
        bool sharpturn = false;
        std::string seq;
        std::string ref;
        // std::string constr;
        while (std::getline(std::cin, seq))
        {
            std::getline(std::cin, ref);
            auto start_time = std::chrono::high_resolution_clock::now();
            alg_3_helper(ref, seq, verbose, dangle);
            auto end_time = std::chrono::high_resolution_clock::now();
            const std::chrono::duration<double, std::milli> time_ms = end_time - start_time;
            printf("alg3 time: %.4f seconds\n", time_ms/1000.f);
        }
        return 0;
    }else if (alg != nullptr && strcmp(alg, "alg1") == 0){ /* alg 1 */
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
            auto start_time = std::chrono::high_resolution_clock::now();
            alg1_helper(seq, ref1, ref2, verbose, dangle);
            auto end_time = std::chrono::high_resolution_clock::now();
            const std::chrono::duration<double, std::milli> time_ms = end_time - start_time;
            printf("alg1 time: %.4f seconds\n", time_ms/1000.f);
        }
    }else if(alg != nullptr && strcmp(alg, "alg1v2") == 0){ /* alg 1 version 2 */
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
            auto start_time = std::chrono::high_resolution_clock::now();
            std::vector<std::string> X = alg_1_v2(ref1, ref2, seq, verbose, dangle);
            auto end_time = std::chrono::high_resolution_clock::now();
            const std::chrono::duration<double, std::milli> time_ms = end_time - start_time;
            printf("X size: %d\n", X.size());
            if (X.size()==0)
                printf("undesignable!\n");
            printf("alg1(v2) time: %.4f seconds\n", time_ms/1000.f);
        }
    }else if (alg != nullptr && strcmp(alg, "alg2") == 0){ /* alg 2 */
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
            auto start_time = std::chrono::high_resolution_clock::now();
            alg_2_helper(ref1, ref2, seq, verbose, dangle);
            auto end_time = std::chrono::high_resolution_clock::now();
            const std::chrono::duration<double, std::milli> time_ms = end_time - start_time;
            printf("alg2 time: %.4f seconds\n", time_ms/1000.f);
        }
    }else if (alg != nullptr && strcmp(alg, "alg2cs") == 0){ /* constrained alg2 */
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
            alg_2_cs_helper(ref1, ref2, seq, verbose, dangle);
        }
    }else if (alg != nullptr && strcmp(alg, "test_diff") == 0){ /* call the function test_diff  */
        verbose = true;
        std::string seq;
        std::string ref1;
        std::string ref2;
        // Read input line by line until EOF (end of file) is reached
        while (std::getline(std::cin, seq)) {
            getline(std::cin, ref1);
            getline(std::cin, ref2);
            test_diff(seq, ref1, ref2, verbose, dangle);
        }
    }else if (strcmp(argv[1], "critical") == 0){  /* critical positions */ 
        std::string ref1;
        std::string ref2;
        while(std::getline(std::cin, ref1)){
            getline(std::cin, ref2);
            std::set<int> critical_positions;
            auto cr_loops = find_critical_plus(ref1, ref2, critical_positions, true);
            auto sp_y = loops2specialhp(cr_loops, ref1.length());
            for(auto sp: sp_y){
                    std::cout<<sp.first<<"\t"<<sp.second<<std::endl;
            }
        }
    }else if (alg != nullptr && strcmp(alg, "eval") == 0){ /* energy evaluation  */
        verbose = true;
        std::string seq;
        std::string ref;
        // Read input line by line until EOF (end of file) is reached
        while (std::getline(std::cin, seq)) {
            getline(std::cin, ref);
            long energy = linear_eval(seq, ref, verbose, dangle);
            printf("total energy: %.2f\n", energy/-100.0);
        }
    }
    return 0;
}