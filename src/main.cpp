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
#include <iomanip>
#include <sstream>
#include <bits/stdc++.h> 

#include "eval.cpp"
#include "csv.cpp"
#include "cxxopts.hpp"
// #include "utils.hpp"
#include "utils.h"
#include "comps.h"
using namespace std;

#define MAX_ENUM 10000000000
#define MAX_CONSTRAINT 100000
#define MAX_SEQ 500
// global variables
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

bool verbose = false;
int dangle = 2;

int y_sub_start;
int y_sub_end;
std::string y_sub;
std::string seq_init;
std::vector<std::string> y_rivals;

std::vector<std::pair<int, int>> pairs_outside;
std::vector<std::pair<int, int>> pairs_inside;

std::vector<std::string> cs_fold(std::string seq, std::string& constr, int beamsize, bool sharpturn, bool verbose, int dangle);
std::vector<std::string> fold(std::string seq, int beamsize, bool sharpturn, bool verbose, int dangle, float energy_delta);
ulong count_enum(std::vector<std::tuple<int, int>>& pairs_diff);

// std::string removeNodeFromTree(TreeNode* node, std::string ref);
// std::string removeMNodeFromTree(TreeNode* node, std::string ref);
// std::string removeEdgeFromTree(TreeNode* node, std::string ref);
// std::string removeTwoNeighbors(TreeNode* node, std::string ref, int n1, int n2);
// std::string removeThreeNeighbors(TreeNode* node, std::string ref, std::vector<int>& powset);

void TreeNode::printTreeEnum(std::string& seq, std::string& y){
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
            node->parent = nodeStack.top();
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

std::string getCurrentTimestamp() {
    // Get the current time point
    auto currentTime = std::chrono::system_clock::now();

    // Convert the time point to a time_t object
    std::time_t currentTime_t = std::chrono::system_clock::to_time_t(currentTime);

    // Convert the time_t to a std::tm structure
    std::tm* timeInfo = std::localtime(&currentTime_t);

    // Create a stringstream to format the timestamp
    std::stringstream timestampStream;
    
    // Use strftime to format the timestamp
    timestampStream << std::put_time(timeInfo, "%Y%m%d%H%M%S");

    // Convert the stringstream to a string
    std::string timestamp = timestampStream.str();

    return timestamp;
}

template <typename T>
std::set<T> setIntersection(const std::set<T>& set1, const std::set<T>& set2) {
    std::set<T> intersection;
    std::set_intersection(set1.begin(), set1.end(), set2.begin(), set2.end(), std::inserter(intersection, intersection.begin()));
    return intersection;
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
        if(X.size() > MAX_SEQ)
            break;
    }
    if (X.size() > MAX_SEQ)
        X.resize(MAX_SEQ);
    std::cout<<"X.size: "<<X.size()<<std::endl;
    for(auto x: X){
        assert (check_compatible(x, ref1));
        std::vector<std::string> subopts = fold(x, 0, false, verbose, dangle_model, 0.0);
        if (isUMFE(subopts, ref1)){
            printf("UMFE found!\n");
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
            y_rivals.clear();
            y_rivals.push_back(y_prime.second.first);
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
                    std::cout<<"constraint counts:\t"<<cs_vec.size()+1<<std::endl; // the latest one was not added to cs_vec
                    for(int i = 0; i<cs_vec.size(); i++)
                        std::cout<<cs_vec[i].seqs->size()<<"\t";
                    std::cout<<cs_new.seqs->size()<<"\t";
                    std::cout<<std::endl;
                    y_rivals.clear();
                    for(int i = 0; i<cs_vec.size(); i++){
                        std::cout<<cs_vec[i].structure<<std::endl;
                        y_rivals.push_back(cs_vec[i].structure);
                    }
                    std::cout<<cs_new.structure<<std::endl;
                    y_rivals.push_back(cs_new.structure);
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
            std::cout<<"constraint counts:\t"<<cs_vec.size()<<std::endl;
            for(int i = 0; i<cs_vec.size(); i++)
                std::cout<<cs_vec[i].seqs->size()<<"\t";
            std::cout<<std::endl;
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
        if(X.size() > MAX_SEQ)
            break;
    }
    if (X.size() > MAX_SEQ)
        X.resize(MAX_SEQ);
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
        std::cout<<x<<std::endl;
        for(std::string ref: subopts){
            std::cout<<ref<<std::endl;
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
            y_sub = ref1;
            y_rivals.clear();
            y_rivals.push_back(y_prime.second.first);
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
                    y_sub = ref1;
                    y_rivals.clear();
                    for(int i = 0; i<cs_vec.size(); i++){
                        std::cout<<cs_vec[i].structure<<std::endl;
                        y_rivals.push_back(cs_vec[i].structure);
                    }
                    std::cout<<cs_new.structure<<std::endl;
                    y_rivals.push_back(cs_new.structure);
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
        std::set<std::string> refs_checked;
        std::vector<Constraint> cs_vec;
        if (X.size() == 0){
            std::cout<<"y :"<<ref1<<std::endl;
            std::cout<<"y':"<<ref2<<std::endl;
            std::cout<<"undesignable!"<<std::endl;
            y_rivals.clear();
            y_rivals.push_back(ref2);
            return "undesignable";
        }else if (X.size() > MAX_CONSTRAINT){
            std::cout<<"too many constraints: "<<X.size()<<"\t"<<"out of "<<ref2<<std::endl;
        }
        else{
            Constraint cs_ref2 = Constraint(&critical_positions, &X);
            cs_ref2.setStructure(ref2);
            cs_vec.push_back(cs_ref2);
        }
        refs_checked.insert(ref2);
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
        std::set<std::string> refs_checked;
        std::vector<Constraint> cs_vec;
        if (X.size() == 0){
            std::cout<<"undesignable!"<<std::endl;
            y_sub = ref1;
            y_rivals.clear();
            y_rivals.push_back(ref2);
            return "undesignable";
        }
        // (IMPORTANT) need to be optimized later
        // else if (X.size() > MAX_CONSTRAINT){  
        //     std::cout<<"too many constraints: "<<X.size()<<"\t"<<"out of "<<ref2<<std::endl;
        // }
        else{
            Constraint cs_ref2 = Constraint(&critical_positions, &X);
            cs_ref2.setStructure(ref2);
            cs_vec.push_back(cs_ref2);
        }
        refs_checked.insert(ref2);
        // std::set<std::string> refs_checked;
        // refs_checked.insert(ref2);
        // std::vector<Constraint> cs_vec;
        // Constraint cs_ref2 = Constraint(&critical_positions, &X);
        // cs_ref2.setStructure(ref2);
        // cs_vec.push_back(cs_ref2);
        return alg_2_cs(ref1, refs_checked, cs_vec, verbose, dangle_model);
    }
    std::cout<<"intial y_prime too bad!"<<std::endl;
    return "intial y_prime too bad";
}

std::string alg_3_helper(std::string& ref, std::string& seq, bool verbose, int dangle){
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
            if (designability == "undesignable"){
                size_t found = ref.find(subrefs[i].first);
                if (found != std::string::npos) {
                    printf("y*: %s\n", ref.c_str());
                    printf("context-constrained undesignable structure y*[%d: %d]\n", found, found+subrefs[i].first.length());
                    std::cout<<subrefs[i].first<<std::endl;
                    std::cout << "Found at position: " << found << std::endl;
                    std::cout<<"substr indices:"<<found<<","<<found+subrefs[i].first.length()<<std::endl;
                } else {
                    std::cout << "substr indices not found." << std::endl;
                }
                return "undesignable";
            }
        }
    }
    return "unknown";
}

std::string alg_3_span_helper(std::string& ref, std::string& seq, std::vector<std::pair<std::string, std::string>>& subrefs, int i, bool verbose, int dangle){
    // TreeNode* root = parseStringToTree(ref);
    // std::vector<std::pair<std::string, std::string>> subrefs;
    // root->printTree(ref, seq, subrefs);
    // std::cout<<"size of sub refs: "<<subrefs.size()<<std::endl;
    // std::sort(subrefs.begin(), subrefs.end(), compareByFirstStringLength);
    // for(int i = 0; i < subrefs.size(); i++){
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
            if (designability == "undesignable"){
                size_t found = ref.find(subrefs[i].first);
                if (found != std::string::npos) {
                    printf("y*: %s\n", ref.c_str());
                    printf("context-constrained undesignable structure y*[%d: %d]\n", found, found+subrefs[i].first.length());
                    std::cout<<subrefs[i].first<<std::endl;
                    std::cout << "Found at position: " << found << std::endl;
                    std::cout<<"substr indices:"<<found<<","<<found+subrefs[i].first.length()<<std::endl;
                } else {
                    std::cout << "substr indices not found." << std::endl;
                }
                return "undesignable";
            }
        }
    // }
    return "unknown";
}

std::string alg_5_cs(std::string& ref1, std::set<std::string>& refs_checked, std::vector<Constraint>& cs_vec, std::string& constr, bool verbose, int dangle_model){ // ref1, ref2, X, is_verbose, dangle_model
    int count_cs = cs_vec.size();
    std::vector<std::pair<ulong, std::pair<std::string, std::string>>> y_primes;
    std::cout<<"inside alg5cs"<<std::endl;
    std::vector<std::string> X;
    // X.push_back(seq_init);
    for(auto cs: cs_vec){
        std::cout<<"cs.seqs size: "<<cs.seqs->size()<<std::endl;
        X.insert(X.end(), cs.seqs->begin(), cs.seqs->end());
        if(X.size() > MAX_SEQ)
            break;
    }
    if (X.size() > MAX_SEQ)
        X.resize(MAX_SEQ);
    std::cout<<"X.size: "<<X.size()<<std::endl;
    // std::string constr(ref1.length(), '?');
    // constr[0] = '(';
    // constr[ref1.length()-1] = ')';
    for(auto x: X){
        assert (check_compatible(x, ref1));
        std::vector<std::string> subopts_raw = cs_fold(x, constr, 0, false, verbose, dangle_model);
        std::vector<std::string> subopts;
        for(std::string subopt: subopts_raw){
            bool flag = true;
            for(int idx = 0; idx < subopt.length(); idx++){
                if (constr[idx]!='?'&&constr[idx]!=subopt[idx]){
                    flag = false;
                    break;
                }
            }
            if(flag)
                subopts.push_back(subopt);
        }
        assert (subopts.size() > 0);
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
            y_sub = ref1;
            y_rivals.clear();
            y_rivals.push_back(y_prime.second.first);
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
                    y_sub = ref1;
                    y_rivals.clear();
                    for(int i = 0; i<cs_vec.size(); i++){
                        std::cout<<cs_vec[i].structure<<std::endl;
                        y_rivals.push_back(cs_vec[i].structure);
                    }
                    std::cout<<cs_new.structure<<std::endl;
                    y_rivals.push_back(cs_new.structure);
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
        return alg_5_cs(ref1, refs_checked, cs_vec, constr, verbose, dangle_model);
    else{
        std::cout<<"no conclusion!"<<std::endl;
        return "no conclusion";
    }
}


std::string alg_5_helper(std::string& ref1, std::string& ref2, std::string&constr, std::string& seq, bool verbose, int dangle_model){
    std::cout << "seq: " << seq << std::endl;
    std::cout << "  y: " << ref1 << std::endl;
    std::cout << " y': " << ref2 << std::endl;
    std::cout << "cst: " << constr << std::endl;
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
        std::set<std::string> refs_checked;
        std::vector<Constraint> cs_vec;
        if (X.size() == 0){
            std::cout<<"undesignable!"<<std::endl;
            y_sub = ref1;
            y_rivals.clear();
            y_rivals.push_back(ref2);
            return "undesignable";
        }else if (X.size() > MAX_CONSTRAINT){
            std::cout<<"too many constraints: "<<X.size()<<"\t"<<"out of "<<ref2<<std::endl;
        }
        else{
            Constraint cs_ref2 = Constraint(&critical_positions, &X);
            cs_ref2.setStructure(ref2);
            cs_vec.push_back(cs_ref2);
        }
        refs_checked.insert(ref2);
        return alg_5_cs(ref1, refs_checked, cs_vec, constr, verbose, dangle_model);
    }
    std::cout<<"intial y_prime too bad!"<<std::endl;
    return "intial y_prime too bad";
}


std::string alg_5_helper_v2(std::string& ref1, std::string& ref2, std::string&constr, std::string& seq, bool verbose, int dangle_model){
    std::cout<< "inside alg_5_helper_v2"<<std::endl;
    std::cout << "seq: " << seq << std::endl;
    std::cout << "  y: " << ref1 << std::endl;
    seq_init = seq;
    std::vector<std::string> mfes = cs_fold(seq, constr, 0, false, verbose, dangle);
    std::string ref_mfe;
    if(isUMFE(mfes, ref1)){
        std::cout<<seq<<std::endl;
        std::cout<<"designable!"<<std::endl;
        return "UMFE";
    }else if(isMFE(mfes, ref1)){
        std::cout<<"mfe designable"<<std::endl;
        if(mfes[0] == ref1)
            ref_mfe = mfes[1];
        else
            ref_mfe = mfes[0];
    }else{
        ref_mfe = mfes[0];
    }
    std::cout << " y': " << ref2 << std::endl;
    std::cout << "cst: " << constr << std::endl;
    std::cout << "mfe: " << ref_mfe << std::endl;

    std::set<int> critical_positions;
    std::set<std::string> refs_checked;
    std::vector<Constraint> cs_vec;

    std::vector<std::vector<int>> cr_loops = find_critical_plus(ref1, ref2, critical_positions, verbose);
    // long delta_energy = diff_eval(seq, cr_loops, verbose, dangle_model);
    std::vector<std::tuple<int, int>> pairs_diff = idx2pair(critical_positions, ref1);
    ulong n_enum = count_enum(pairs_diff);
    std::cout<<"enumeration count: "<<n_enum<<std::endl;
    // check the ref from removing internal pairs
    std::vector<std::string> X;
    Constraint cs_ref2 = Constraint(&critical_positions, &X);
    if(n_enum > 0 && n_enum < MAX_ENUM){
        std::cout<<"alg 1"<<std::endl;
        X = alg_1_v2(ref1, ref2, seq, verbose, dangle_model);
        printf("X size: %d\n", X.size());
        if (X.size() == 0){
            std::cout<<"undesignable!"<<std::endl;
            y_sub = ref1;
            y_rivals.clear();
            y_rivals.push_back(ref2);
            return "undesignable";
        }else if (X.size() > MAX_CONSTRAINT){
            std::cout<<"too many constraints: "<<X.size()<<"\t"<<"out of "<<ref2<<std::endl;
            refs_checked.insert(ref2);
        }
        else{
            cs_ref2 = Constraint(&critical_positions, &X);
            cs_ref2.setStructure(ref2);
            // cs_vec.push_back(cs_ref2);
        }
        // refs_checked.insert(ref2);
        // return alg_5_cs(ref1, refs_checked, cs_vec, constr, verbose, dangle_model);
    }
    // finish the first check

    // check the ref from mfe
    std::set<int> critical_positions_mfe;
    cr_loops = find_critical_plus(ref1, ref_mfe, critical_positions_mfe, verbose);
    // delta_energy = diff_eval(seq, cr_loops, verbose, dangle_model);
    pairs_diff = idx2pair(critical_positions_mfe, ref1);
    n_enum = count_enum(pairs_diff);
    std::cout<<"enumeration count: "<<n_enum<<std::endl;
    std::vector<std::string> X_mfe;
    if(n_enum > 0 && n_enum < MAX_ENUM){
        std::cout<<"alg 1"<<std::endl;
        // auto X = alg_1(ref1, ref2, cr_loops, pairs_diff, seq, verbose, dangle_model);
        X_mfe = alg_1_v2(ref1, ref_mfe, seq, verbose, dangle_model);
        printf("X_mfe size: %d\n", X_mfe.size());
        // std::set<std::string> refs_checked;
        // std::vector<Constraint> cs_vec;
        if (X_mfe.size() == 0){
            std::cout<<"undesignable!"<<std::endl;
            y_sub = ref1;
            y_rivals.clear();
            y_rivals.push_back(ref_mfe);
            return "undesignable";
        }else if (X_mfe.size() > MAX_CONSTRAINT){
            std::cout<<"too many constraints: "<<X_mfe.size()<<"\t"<<"out of "<<ref_mfe<<std::endl;
            if(cs_ref2.seqs->size()){
                cs_vec.push_back(cs_ref2);
                refs_checked.insert(ref2);
            }
        }
        else{
            Constraint cs_ref_mfe = Constraint(&critical_positions_mfe, &X_mfe);
            cs_ref_mfe.setStructure(ref_mfe);
            cs_vec.push_back(cs_ref_mfe);
        }
        refs_checked.insert(ref_mfe);
        return alg_5_cs(ref1, refs_checked, cs_vec, constr, verbose, dangle_model);
    }
    std::cout<<"intial y_prime too bad!"<<std::endl;
    return "intial y_prime too bad";
}


std::vector<std::string> cs_fold(std::string seq, std::string& constr, int beamsize, bool sharpturn, bool verbose, int dangle){
    return subopt(seq, constr);
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

        // Use std::find to search for the value
        if (std::find(subopts.begin(), subopts.end(), result.structure) == subopts.end()){
            subopts.push_back(result.structure);
        } 
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


std::string alg1_helper(std::string& seq, std::string& ref1, std::string& ref2, bool is_verbose, int dangle_model) {
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
        // auto X = alg_1(ref1, ref2, cr_loops, pairs_diff, seq, is_verbose, dangle_model);
        auto X = alg_1_v2(ref1, ref2, seq, verbose, dangle_model);
        printf("X size: %d\n", X.size());
        if (X.size()==0){
            printf("undesignable!\n");
            return "undesignable";
        }
        else{
            Constraint cs(&critical_positions, &X);
            printf("constraint size: %d\n", cs.seqs->size());
            return "unknown";
        }
    }
    std::cout<<"intial y_prime too bad!"<<std::endl;
    return "intial y_prime too bad";
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

void csv_process(std::string csv, std::string alg){
    auto df = read_csv(csv.c_str());
    printf("df shape: %d, %d\n", df.size(), df[0].size());
    std::vector<std::string> records;
    // Specify the file name
    std::string fileName = csv + "." + alg + ".log."+getCurrentTimestamp()+".txt";
    std::string file_m1 = replaceFileExtension(csv, "m1");
    std::map<std::string, std::set<std::string>> id2m1 = readMotif(file_m1.c_str());
    // Open the file for writing
    std::ofstream outputFile(fileName);

    // Check if the file is open
    if (!outputFile.is_open()) {
        std::cerr << "Error opening the file: " << fileName << std::endl;
        return;
    }
    
    for(int i = 1; i < df.size(); i++){
        auto row = df[i];
        if (row[7] != "None"){
            std::cout<<"Puzzle ID: "<<row[1]<<std::endl;
            std::cout<<"Puzzle name: "<<row[2]<<std::endl;
            std::cout<<row[8]<<std::endl;
            std::cout<<row[4]<<std::endl;
            std::cout<<row[7]<<std::endl;
            std::cout<<std::endl;

            std::string puzzle_id = row[1];
            std::string seq = row[8];
            std::string y_star = row[4];
            std::string y_prim = row[7];

            if (alg == "1"){
                auto start_time = std::chrono::high_resolution_clock::now();
                std::string result = alg1_helper(seq, y_star, y_prim, verbose, dangle);
                auto end_time = std::chrono::high_resolution_clock::now();
                const std::chrono::duration<double, std::milli> time_ms = end_time - start_time;
                printf("alg1(v2) time: %.4f seconds\n", time_ms/1000.f);
                if (result == "undesignable"){
                    std::string r = row[1]+","+y_star+","+y_prim;
                    std::set<std::pair<int, int>> pairset_star = ref2pairset(y_star);
                    std::set<std::pair<int, int>> pairset_prim = ref2pairset(y_prim);
                    for (std::pair<int, int> p: pairset_prim){
                        assert (pairset_star.find(p) != pairset_star.end());
                    }

                    std::set<std::pair<int, int>> pairset_diff;
                    std::set_difference(pairset_star.begin(), pairset_star.end(), pairset_prim.begin(), pairset_prim.end(),
                    std::inserter(pairset_diff, pairset_diff.begin()));
                    r += ","+std::to_string(pairset_diff.size());
                    for (auto pair: pairset_diff)
                        r += ","+std::to_string(pair.first)+","+std::to_string(pair.second);
                    records.push_back(r);
                } 
            }
            if (alg == "2"){ // && (puzzle_id == "60" || puzzle_id == "88")
                auto start_time = std::chrono::high_resolution_clock::now();
                std::string result = alg_2_helper(y_star, y_prim, seq, verbose, dangle);
                auto end_time = std::chrono::high_resolution_clock::now();
                const std::chrono::duration<double, std::milli> time_ms = end_time - start_time;
                printf("alg2 time: %.4f seconds\n", time_ms/1000.f);
                if (result == "undesignable"){
                    std::string r = row[1]+","+y_star+","+std::to_string(y_rivals.size());
                    for(auto rival: y_rivals)
                        r += ","+rival;
                    records.push_back(r);
                }
            }
            if (alg == "3"){
                auto start_time = std::chrono::high_resolution_clock::now();
                std::string result = alg_3_helper(y_star, seq, verbose, dangle);
                auto end_time = std::chrono::high_resolution_clock::now();
                const std::chrono::duration<double, std::milli> time_ms = end_time - start_time;
                printf("alg3 time: %.4f seconds\n", time_ms/1000.f);
                if (result == "undesignable"){
                    size_t found = y_star.find(y_sub);
                    assert (found != std::string::npos);
                    printf("y*: %s\n", y_star.c_str());
                    printf("context-constrained undesignable structure y*[%d: %d]\n", found, found+y_sub.length());
                    std::cout<<y_sub<<std::endl;
                    std::cout << "Found at position: " << found << std::endl;
                    std::cout<<"substr indices:"<<found<<","<<found+y_sub.length()<<std::endl;

                    std::string r = row[1]+","+y_star+","+std::to_string(found)+","+std::to_string(found+y_sub.length())+","+y_sub+","+std::to_string(y_rivals.size());
                    for(auto rival: y_rivals)
                        r += ","+rival;
                    records.push_back(r);
                }
            }
            if (alg == "3_v2"){
                TreeNode* root = parseStringToTree(y_star);
                std::vector<std::pair<std::string, std::string>> subrefs;
                root->printTree(y_star, seq, subrefs);
                std::cout<<"size of sub refs: "<<subrefs.size()<<std::endl;
                std::sort(subrefs.begin(), subrefs.end(), compareByFirstStringLength);
                std::vector<std::pair<int, int>> span_vec;
                for(int i = 0; i < subrefs.size(); i++){
                    auto start_time = std::chrono::high_resolution_clock::now();
                    std::string result = alg_3_span_helper(y_star, seq, subrefs, i, verbose, dangle);
                    auto end_time = std::chrono::high_resolution_clock::now();
                    const std::chrono::duration<double, std::milli> time_ms = end_time - start_time;
                    printf("alg3(v2) time: %.4f seconds\n", time_ms/1000.f);
                    if (result == "undesignable"){
                        std::vector<int> pos_vec = findAllOccurrences(y_star, y_sub);
                        for (auto found: pos_vec){
                            // size_t found = y_star.find(y_sub);
                            assert (found != std::string::npos);
                            int found_end = found+y_sub.length();
                            bool duplicate = false;
                            for (int j = 0; j < span_vec.size(); j++){
                                if ( found <= span_vec[j].first && found_end >= span_vec[j].second ){  // the equality may be worth a double-check
                                    duplicate = true;
                                    break;
                                }else if ( found >= span_vec[j].first && found_end <= span_vec[j].second ){  // the equality may be worth a double-check
                                    assert (false); // subref is ordered by length from short to long
                                }
                            }
                            if (!duplicate){
                                printf("y*: %s\n", y_star.c_str());
                                printf("context-constrained undesignable structure y*[%d: %d]\n", found, found+y_sub.length());
                                std::cout<<y_sub<<std::endl;
                                std::cout << "Found at position: " << found << std::endl;
                                std::cout<<"substr indices:"<<found<<","<<found+y_sub.length()<<std::endl;

                                std::string r = row[1]+","+y_star+","+std::to_string(found)+","+std::to_string(found+y_sub.length())+","+y_sub+","+std::to_string(y_rivals.size());
                                for(auto rival: y_rivals)
                                    r += ","+rival;
                                records.push_back(r);
                                outputFile << r << std::endl;
                                span_vec.push_back(std::make_pair(found, found_end));
                            }
                        }
                    }
                }
            }
            if (alg == "edge"){
                auto start_time = std::chrono::high_resolution_clock::now();
                std::vector<LoopComplex> lc_list;
                TreeNode* root = parseStringToTree(y_star);
                tree2Edges(root, y_star, lc_list);
                printf("lc_list size: %d\n", lc_list.size());
                // Sort the vector using a lambda expression
                std::sort(lc_list.begin(), lc_list.end(), [](const LoopComplex &a, const LoopComplex &b) {
                    return a.count_uk < b.count_uk;});
                for (auto lc: lc_list){
                    std::string target = y_star.substr(lc.start, lc.end-lc.start+1);
                    std::string subseq = seq.substr(lc.start, lc.end-lc.start+1);
                    printf(" count: %d\n", lc.count_uk);
                    printf("target: %s\n", target.c_str());
                    printf("   ref: %s\n", lc.ref.c_str());
                    printf("constr: %s\n", lc.constr.c_str());

                    std::string result = alg_5_helper_v2(target, lc.ref, lc.constr, subseq, verbose, dangle);
                    if (result == "undesignable"){
                        std::cout<<"undesignable!"<<std::endl;
                        auto end_time = std::chrono::high_resolution_clock::now();
                        const std::chrono::duration<double, std::milli> time_ms = end_time - start_time;
                        printf("time cost: %.4f seconds\n", time_ms/1000.f);
                        size_t found = y_star.find(y_sub);
                        assert (found != std::string::npos);
                        int found_end = found+y_sub.length();
                        std::string r = row[1]+","+y_star+",1,"+std::to_string(lc.node->first)+","+std::to_string(lc.node->second)+","+y_sub+","+std::to_string(y_rivals.size());
                        for(auto rival: y_rivals)
                            r += ","+rival;
                        // Convert duration to seconds and then cast to float
                        float time_seconds = std::chrono::duration_cast<std::chrono::duration<float>>(time_ms).count();
                        r += " time:" + fl2str(time_seconds) + ";";
                        std::cout<<r<<std::endl;
                        records.push_back(r);
                        std::string id = row[1] + "_" + alg;
                        std::string args4plot = compose_args4plot(id, y_star, lc.ps_outside, lc.ps_inside);
                        outputFile << r << std::endl;
                        outputFile << args4plot << std::endl;
                        std::string pairstring = id + ":" + compose_pairstr(lc.ps_inside, lc.ps_outside);
                        outputFile << pairstring <<endl;
                        // break;
                    }
                    printf("\n");
                }
            }
            if (alg == "loop"){
                seq = tg_init(y_star);
                std::vector<LoopComplex> lc_list;
                TreeNode* root = parseStringToTree(y_star);
                tree2Loops(root, y_star, lc_list);
                printf("lc_list size: %d\n", lc_list.size());

                // Sort the vector using a lambda expression
                std::sort(lc_list.begin(), lc_list.end(), [](const LoopComplex &a, const LoopComplex &b) {
                    return a.count_uk < b.count_uk;});
                for (auto lc: lc_list){
                    std::string target = y_star.substr(lc.start, lc.end-lc.start+1);
                    std::string subseq = seq.substr(lc.start, lc.end-lc.start+1);
                    printf(" count: %d\n", lc.count_uk);
                    printf("target: %s\n", target.c_str());
                    printf("   ref: %s\n", lc.ref.c_str());
                    printf("constr: %s\n", lc.constr.c_str());

                    std::string result = alg_5_helper_v2(target, lc.ref, lc.constr, subseq, verbose, dangle);
                    if (result == "undesignable"){
                        std::cout<<"undesignable!"<<std::endl;
                        int count_pairs = lc.node->children.size() + 1;
                        std::string r = row[1]+","+y_star+","+y_prim+","+std::to_string(count_pairs)+","+std::to_string(lc.node->first)+","+std::to_string(lc.node->second);
                        for(auto child: lc.node->children)
                            r += ","+std::to_string(child->first)+","+std::to_string(child->second);
                        r = r + ","+y_sub+","+std::to_string(y_rivals.size());
                        for(auto rival: y_rivals)
                            r += ","+rival;
                        std::cout<<r<<std::endl;
                        records.push_back(r);
                        outputFile << r << std::endl;
                        // break;
                        // return;
                    }
                    printf("\n");
                }
            }
            if (alg == "mloop"){
                std::vector<LoopComplex> lc_list;
                TreeNode* root = parseStringToTree(y_star);
                tree2MLoops(root, y_star, lc_list);
                printf("lc_list size: %d\n", lc_list.size());

                // Sort the vector using a lambda expression
                std::sort(lc_list.begin(), lc_list.end(), [](const LoopComplex &a, const LoopComplex &b) {
                    return a.count_uk < b.count_uk;});
                for (auto lc: lc_list){
                    std::string target = y_star.substr(lc.start, lc.end-lc.start+1);
                    std::string subseq = seq.substr(lc.start, lc.end-lc.start+1);
                    printf(" count: %d\n", lc.count_uk);
                    printf("target: %s\n", target.c_str());
                    printf("   ref: %s\n", lc.ref.c_str());
                    printf("constr: %s\n", lc.constr.c_str());

                    std::string result = alg_5_helper_v2(target, lc.ref, lc.constr, subseq, verbose, dangle);
                    if (result == "undesignable"){
                        std::cout<<"undesignable!"<<std::endl;
                        int count_pairs = lc.node->children.size() + 1;
                        std::string r = row[1]+","+y_star+","+y_prim+","+std::to_string(count_pairs)+","+std::to_string(lc.node->first)+","+std::to_string(lc.node->second);
                        for(auto child: lc.node->children)
                            r += ","+std::to_string(child->first)+","+std::to_string(child->second);
                        r = r + ","+y_sub+","+std::to_string(y_rivals.size());
                        for(auto rival: y_rivals)
                            r += ","+rival;
                        std::cout<<r<<std::endl;
                        records.push_back(r);
                        outputFile << r << std::endl;
                        // break;
                        // return;
                    }
                    printf("\n");
                }
            }
            if (alg == "neighbor2" || alg == "neighbor3" ){
                auto start_time = std::chrono::high_resolution_clock::now();
                std::vector<LoopComplex> lc_list;
                TreeNode* root = parseStringToTree(y_star);
                int max_internal = max_single(root);
                if(max_internal > 30){
                    std::string r = puzzle_id+","+std::to_string(max_internal);
                    records.push_back(r);
                    outputFile << r << std::endl;
                    continue;
                }
                if (alg == "neighbor2")
                    tree2TwoNeighbor(root, y_star, lc_list);
                else
                    tree2ThreeNeighbor(root, y_star, lc_list);
                printf("lc_list size: %d\n", lc_list.size());

                // Sort the vector using a lambda expression
                std::sort(lc_list.begin(), lc_list.end(), [](const LoopComplex &a, const LoopComplex &b) {
                    return a.count_uk < b.count_uk;});
                for (auto lc: lc_list){
                    bool isSubset = false;
                    if(id2m1.find(puzzle_id) != id2m1.end()){
                        for(auto pair: lc.ps_inside){
                            std::string pstr = std::to_string(pair.first) + "," + std::to_string(pair.second);
                            if(id2m1[puzzle_id].count(pstr) > 0)
                                isSubset = true;
                        }
                    }
                    if(isSubset){
                        std::cout<<"the current motif contains an undesignable motif!"<<std::endl;
                        continue;
                    }
                    std::string target = y_star.substr(lc.start, lc.end-lc.start+1);
                    std::string subseq = seq.substr(lc.start, lc.end-lc.start+1);
                    printf(" count: %d\n", lc.count_uk);
                    printf("target: %s\n", target.c_str());
                    printf("   ref: %s\n", lc.ref.c_str());
                    printf("constr: %s\n", lc.constr.c_str());

                    std::string result = alg_5_helper_v2(target, lc.ref, lc.constr, subseq, verbose, dangle);
                    if (result == "undesignable"){
                        std::cout<<"undesignable!"<<std::endl;
                        auto end_time = std::chrono::high_resolution_clock::now();
                        const std::chrono::duration<double, std::milli> time_ms = end_time - start_time;
                        printf("time cost: %.4f seconds\n", time_ms/1000.f);
                        int count_pairs = lc.node->children.size() + 1;
                        std::string r = puzzle_id+","+y_star+","+std::to_string(count_pairs)+","+std::to_string(lc.node->first)+","+std::to_string(lc.node->second);
                        for(auto child: lc.node->children)
                            r += ","+std::to_string(child->first)+","+std::to_string(child->second);
                        r = r + ","+y_sub+","+std::to_string(y_rivals.size());
                        for(auto rival: y_rivals)
                            r += ","+rival;
                        // Convert duration to seconds and then cast to float
                        float time_seconds = std::chrono::duration_cast<std::chrono::duration<float>>(time_ms).count();
                        r += " time:" + fl2str(time_seconds) + ";";
                        std::cout<<r<<std::endl;
                        records.push_back(r);
                        std::string id = puzzle_id + "_" + alg;
                        std::string args4plot = compose_args4plot(id, y_star, lc.ps_outside, lc.ps_inside);
                        outputFile << r << std::endl;
                        outputFile << args4plot <<std::endl;
                        std::string pairstring = id + ":" + compose_pairstr(lc.ps_inside, lc.ps_outside);
                        outputFile << pairstring <<endl;
                        // break;
                        // return;
                    }
                    printf("\n");
                }
            }
            if (alg == "nb2plot" || alg == "nb3plot"){
                auto start_time = std::chrono::high_resolution_clock::now();
                std::vector<LoopComplex> lc_list;
                TreeNode* root = parseStringToTree(y_star);
                int max_internal = max_single(root);
                if(max_internal > 30){
                    std::string r = puzzle_id+","+std::to_string(max_internal);
                    records.push_back(r);
                    outputFile << r << std::endl;
                    continue;
                }
                printf("alg: %s\n", alg.c_str());
                if (alg == "nb2plot")
                    tree2TwoNeighbor(root, y_star, lc_list);
                else
                    tree2ThreeNeighbor(root, y_star, lc_list);
                printf("lc_list size: %d\n", lc_list.size());

                // Sort the vector using a lambda expression
                std::sort(lc_list.begin(), lc_list.end(), [](const LoopComplex &a, const LoopComplex &b) {
                    return a.count_uk < b.count_uk;});
                for (auto lc: lc_list){
                    printf("alg: %s\n", alg.c_str());
                    bool isSubset = false;
                    if(id2m1.find(puzzle_id) != id2m1.end()){
                        for(auto pair: lc.ps_inside){
                            std::string pstr = std::to_string(pair.first) + "," + std::to_string(pair.second);
                            if(id2m1[puzzle_id].count(pstr) > 0)
                                isSubset = true;
                        }
                    }
                    if(isSubset){
                        std::cout<<"the current motif contains an undesignable motif!"<<std::endl;
                        continue;
                    }
                    std::string target = y_star.substr(lc.start, lc.end-lc.start+1);
                    std::string subseq = seq.substr(lc.start, lc.end-lc.start+1);
                    printf(" count: %d\n", lc.count_uk);
                    printf("target: %s\n", target.c_str());
                    printf("   ref: %s\n", lc.ref.c_str());
                    printf("constr: %s\n", lc.constr.c_str());
                    {
                        std::string id = puzzle_id + "_" + alg;
                        std::string args4plot = compose_args4plot(id, y_star, lc.ps_outside, lc.ps_inside);
                        std::cout<< args4plot <<std::endl;
                        outputFile << args4plot <<std::endl;
                        std::string pairstring = id + ":" + compose_pairstr(lc.ps_inside, lc.ps_outside);
                        outputFile << pairstring << std::endl;
                        std::cout<< pairstring <<std::endl;
                    }
                    printf("\n");
                }
            }
            if (alg == "dsedge"){
                std::vector<LoopComplex> lc_list;
                TreeNode* root = parseStringToTree(y_star);
                int max_internal = max_single(root);
                if(max_internal > 30){
                    std::string r = puzzle_id+","+std::to_string(max_internal);
                    records.push_back(r);
                    outputFile << r << std::endl;
                    continue;
                }
                tree2Edges(root, y_star, lc_list);
                printf("lc_list size: %d\n", lc_list.size());
                // Sort the vector using a lambda expression
                std::sort(lc_list.begin(), lc_list.end(), [](const LoopComplex &a, const LoopComplex &b) {
                    return a.count_uk < b.count_uk;});
                std::string r = puzzle_id+";"+y_star+";";
                std::vector<std::pair<int, int>> pairs_ds;
                std::vector<std::pair<int, int>> pairs_ud;
                for (auto lc: lc_list){
                    std::string target = y_star.substr(lc.start, lc.end-lc.start+1);
                    std::string subseq = seq.substr(lc.start, lc.end-lc.start+1);
                    printf(" count: %d\n", lc.count_uk);
                    printf("target: %s\n", target.c_str());
                    printf("   ref: %s\n", lc.ref.c_str());
                    printf("constr: %s\n", lc.constr.c_str());
                    std::string result = alg_5_helper_v2(target, lc.ref, lc.constr, subseq, verbose, dangle);
                    if (result == "UMFE"){
                        std::cout<<"UMFE!"<<std::endl;
                        r += std::to_string(lc.left) + "," + std::to_string(lc.right) + ";";
                        pairs_ds.push_back(std::make_pair(lc.left, lc.right));
                        // break;
                    }else if (result == "undesignable"){
                        pairs_ud.push_back(std::make_pair(lc.left, lc.right));
                    }
                    printf("\n");
                }
                std::cout<<r<<std::endl;
                records.push_back(r);
                std::string id = puzzle_id + "_" + "dsedge";
                std::string pairsplot = compose_pairsplot(id, y_star, pairs_ds, pairs_ud);
                outputFile << r << std::endl;
                outputFile << pairsplot << std::endl;
            }
            if (alg == "csgen"){
                std::vector<LoopComplex> lc_list;
                TreeNode* root = parseStringToTree(y_star);
                int max_internal = max_single(root);
                if(max_internal > 30){
                    std::string r = puzzle_id+","+std::to_string(max_internal);
                    records.push_back(r);
                    outputFile << r << std::endl;
                    continue;
                }
                tree2Loops(root, y_star, lc_list);
                printf("lc_list size: %d\n", lc_list.size());

                // Sort the vector using a lambda expression
                std::sort(lc_list.begin(), lc_list.end(), [](const LoopComplex &a, const LoopComplex &b) {
                    return a.count_uk < b.count_uk;});
                for (auto lc: lc_list){
                    std::string target = y_star.substr(lc.start, lc.end-lc.start+1);
                    std::string subseq = seq.substr(lc.start, lc.end-lc.start+1);
                    printf(" count: %d\n", lc.count_uk);
                    printf("target: %s\n", target.c_str());
                    printf("   ref: %s\n", lc.ref.c_str());
                    printf("constr: %s\n", lc.constr.c_str());
                    outputFile << subseq << std::endl;
                    outputFile << lc.constr << std::endl;
                    printf("\n");
                }
            }
        }
    }
    for (auto r: records)
        std::cout<<r<<std::endl;
    
    // Close the file
    outputFile.close();
    std::cout << "Strings written to file: " << fileName << std::endl;
}

void txt_process(std::string txt, std::string alg){
    std::vector<std::string> lines = readLinesFromFile(txt);
    printf("line count: %d\n", lines.size());
    // Specify the file name
    std::string fileName = txt + "." + alg + ".log."+getCurrentTimestamp()+".txt";
    // Open the file for writing
    std::ofstream outputFile(fileName);

    // Check if the file is open
    if (!outputFile.is_open()) {
        std::cerr << "Error opening the file: " << fileName << std::endl;
        return;
    }
    std::vector<std::string> records;
    for(int i = 0; i < lines.size(); i++){
        if (true){
            std::cout<<"Puzzle ID: "<<i<<std::endl;
            std::cout<<"Puzzle name: "<<i<<std::endl;

            std::string puzzle_id = std::to_string(i);
            std::string y_star = lines[i];
            std::string seq = tg_init(y_star);
            // std::string y_prim = row[7];

            if (alg == "inspect"){
                 TreeNode* root = parseStringToTree(y_star);
                 int max_internal = max_single(root);
                 std::string r = puzzle_id+","+std::to_string(max_internal);
                 records.push_back(r);
                 outputFile << r << std::endl;
            }
            if (alg == "mdegree"){
                 TreeNode* root = parseStringToTree(y_star);
                 std::string mloops = ml_degree(root);
                 std::string r = puzzle_id + "," + y_star + ";" + mloops;
                 records.push_back(r);
                 outputFile << r << std::endl;
            }
            // if (alg == "edge"){
            //     std::vector<LoopComplex> lc_list;
            //     TreeNode* root = parseStringToTree(y_star);
            //     int max_internal = max_single(root);
            //     if(max_internal > 30){
            //         std::string r = puzzle_id+","+std::to_string(max_internal);
            //         records.push_back(r);
            //         outputFile << r << std::endl;
            //         continue;
            //     }
            //     tree2Edges(root, y_star, lc_list);
            //     printf("lc_list size: %d\n", lc_list.size());
            //     // Sort the vector using a lambda expression
            //     std::sort(lc_list.begin(), lc_list.end(), [](const LoopComplex &a, const LoopComplex &b) {
            //         return a.count_uk < b.count_uk;});
            //     for (auto lc: lc_list){
            //         std::string target = y_star.substr(lc.start, lc.end-lc.start+1);
            //         std::string subseq = seq.substr(lc.start, lc.end-lc.start+1);
            //         printf(" count: %d\n", lc.count_uk);
            //         printf("target: %s\n", target.c_str());
            //         printf("   ref: %s\n", lc.ref.c_str());
            //         printf("constr: %s\n", lc.constr.c_str());

            //         std::string result = alg_5_helper_v2(target, lc.ref, lc.constr, subseq, verbose, dangle);
            //         if (result == "undesignable"){
            //             std::cout<<"undesignable!"<<std::endl;
            //             size_t found = y_star.find(y_sub);
            //             assert (found != std::string::npos);
            //             int found_end = found+y_sub.length();
            //             std::string r = puzzle_id+","+y_star+",1,"+std::to_string(lc.node->first)+","+std::to_string(lc.node->second)+","+y_sub+","+std::to_string(y_rivals.size());
            //             for(auto rival: y_rivals)
            //                 r += ","+rival;
            //             std::cout<<r<<std::endl;
            //             records.push_back(r);
            //             outputFile << r << std::endl;
            //             // break;
            //         }
            //         printf("\n");
            //     }
            // }
            if (alg == "neighbor2"){
                auto start_time = std::chrono::high_resolution_clock::now();
                std::vector<LoopComplex> lc_list;
                TreeNode* root = parseStringToTree(y_star);
                int max_internal = max_single(root);
                if(max_internal > 30){
                    std::string r = puzzle_id+","+std::to_string(max_internal);
                    records.push_back(r);
                    outputFile << r << std::endl;
                    continue;
                }
                tree2TwoNeighbor(root, y_star, lc_list);
                printf("lc_list size: %d\n", lc_list.size());

                // Sort the vector using a lambda expression
                std::sort(lc_list.begin(), lc_list.end(), [](const LoopComplex &a, const LoopComplex &b) {
                    return a.count_uk < b.count_uk;});
                for (auto lc: lc_list){
                    std::string target = y_star.substr(lc.start, lc.end-lc.start+1);
                    std::string subseq = seq.substr(lc.start, lc.end-lc.start+1);
                    printf(" count: %d\n", lc.count_uk);
                    printf("target: %s\n", target.c_str());
                    printf("   ref: %s\n", lc.ref.c_str());
                    printf("constr: %s\n", lc.constr.c_str());

                    std::string result = alg_5_helper_v2(target, lc.ref, lc.constr, subseq, verbose, dangle);
                    if (result == "undesignable"){
                        std::cout<<"undesignable!"<<std::endl;
                        auto end_time = std::chrono::high_resolution_clock::now();
                        const std::chrono::duration<double, std::milli> time_ms = end_time - start_time;
                        printf("time cost: %.4f seconds\n", time_ms/1000.f);
                        int count_pairs = lc.node->children.size() + 1;
                        std::string r = puzzle_id+","+y_star+","+std::to_string(count_pairs)+","+std::to_string(lc.node->first)+","+std::to_string(lc.node->second);
                        for(auto child: lc.node->children)
                            r += ","+std::to_string(child->first)+","+std::to_string(child->second);
                        r = r + ","+y_sub+","+std::to_string(y_rivals.size());
                        for(auto rival: y_rivals)
                            r += ","+rival;
                        // Convert duration to seconds and then cast to float
                        float time_seconds = std::chrono::duration_cast<std::chrono::duration<float>>(time_ms).count();
                        r += " time:" + fl2str(time_seconds) + ";";
                        std::cout<<r<<std::endl;
                        records.push_back(r);
                        std::string id = puzzle_id + "_" + alg;
                        std::string args4plot = compose_args4plot(id, y_star, lc.ps_outside, lc.ps_inside);
                        outputFile << r << std::endl;
                        outputFile << args4plot <<std::endl;
                        // break;
                        // return;
                    }
                    printf("\n");
                }
            }
            if (alg == "dsedge"){
                std::vector<LoopComplex> lc_list;
                TreeNode* root = parseStringToTree(y_star);
                int max_internal = max_single(root);
                if(max_internal > 30){
                    std::string r = puzzle_id+","+std::to_string(max_internal);
                    records.push_back(r);
                    outputFile << r << std::endl;
                    continue;
                }
                tree2Edges(root, y_star, lc_list);
                printf("lc_list size: %d\n", lc_list.size());
                // Sort the vector using a lambda expression
                std::sort(lc_list.begin(), lc_list.end(), [](const LoopComplex &a, const LoopComplex &b) {
                    return a.count_uk < b.count_uk;});
                std::string r = puzzle_id+";"+y_star+";";
                std::vector<std::pair<int, int>> pairs_ds;
                std::vector<std::pair<int, int>> pairs_ud;
                for (auto lc: lc_list){
                    std::string target = y_star.substr(lc.start, lc.end-lc.start+1);
                    std::string subseq = seq.substr(lc.start, lc.end-lc.start+1);
                    printf(" count: %d\n", lc.count_uk);
                    printf("target: %s\n", target.c_str());
                    printf("   ref: %s\n", lc.ref.c_str());
                    printf("constr: %s\n", lc.constr.c_str());
                    std::string result = alg_5_helper_v2(target, lc.ref, lc.constr, subseq, verbose, dangle);
                    if (result == "UMFE"){
                        std::cout<<"UMFE!"<<std::endl;
                        r += std::to_string(lc.left) + "," + std::to_string(lc.right) + ";";
                        pairs_ds.push_back(std::make_pair(lc.left, lc.right));
                        // break;
                    }else if (result == "undesignable"){
                        pairs_ud.push_back(std::make_pair(lc.left, lc.right));
                    }
                    printf("\n");
                }
                std::cout<<r<<std::endl;
                records.push_back(r);
                std::string id = puzzle_id + "_" + "dsedge";
                std::string pairsplot = compose_pairsplot(id, y_star, pairs_ds, pairs_ud);
                outputFile << r << std::endl;
                outputFile << pairsplot << std::endl;
            }
            if (alg == "loop"){
                std::vector<LoopComplex> lc_list;
                TreeNode* root = parseStringToTree(y_star);
                int max_internal = max_single(root);
                if(max_internal > 30){
                    std::string r = puzzle_id+","+std::to_string(max_internal);
                    records.push_back(r);
                    outputFile << r << std::endl;
                    continue;
                }
                tree2Loops(root, y_star, lc_list);
                printf("lc_list size: %d\n", lc_list.size());

                // Sort the vector using a lambda expression
                std::sort(lc_list.begin(), lc_list.end(), [](const LoopComplex &a, const LoopComplex &b) {
                    return a.count_uk < b.count_uk;});
                for (auto lc: lc_list){
                    std::string target = y_star.substr(lc.start, lc.end-lc.start+1);
                    std::string subseq = seq.substr(lc.start, lc.end-lc.start+1);
                    printf(" count: %d\n", lc.count_uk);
                    printf("target: %s\n", target.c_str());
                    printf("   ref: %s\n", lc.ref.c_str());
                    printf("constr: %s\n", lc.constr.c_str());

                    std::string result = alg_5_helper_v2(target, lc.ref, lc.constr, subseq, verbose, dangle);
                    if (result == "undesignable"){
                        std::cout<<"undesignable!"<<std::endl;
                        int count_pairs = lc.node->children.size() + 1;
                        std::string r = puzzle_id+","+y_star+","+std::to_string(count_pairs)+","+std::to_string(lc.node->first)+","+std::to_string(lc.node->second);
                        for(auto child: lc.node->children)
                            r += ","+std::to_string(child->first)+","+std::to_string(child->second);
                        r = r + ","+y_sub+","+std::to_string(y_rivals.size());
                        for(auto rival: y_rivals)
                            r += ","+rival;
                        std::cout<<r<<std::endl;
                        records.push_back(r);
                        outputFile << r << std::endl;
                        // break;
                        // return;
                    }
                    printf("\n");
                }
            }
            if (alg == "mloop"){
                std::vector<LoopComplex> lc_list;
                TreeNode* root = parseStringToTree(y_star);
                int max_internal = max_single(root);
                if(max_internal > 30){
                    std::string r = puzzle_id+","+std::to_string(max_internal);
                    records.push_back(r);
                    outputFile << r << std::endl;
                    continue;
                }
                tree2MLoops(root, y_star, lc_list);
                printf("lc_list size: %d\n", lc_list.size());

                // Sort the vector using a lambda expression
                std::sort(lc_list.begin(), lc_list.end(), [](const LoopComplex &a, const LoopComplex &b) {
                    return a.count_uk < b.count_uk;});
                for (auto lc: lc_list){
                    std::string target = y_star.substr(lc.start, lc.end-lc.start+1);
                    std::string subseq = seq.substr(lc.start, lc.end-lc.start+1);
                    printf(" count: %d\n", lc.count_uk);
                    printf("target: %s\n", target.c_str());
                    printf("   ref: %s\n", lc.ref.c_str());
                    printf("constr: %s\n", lc.constr.c_str());

                    std::string result = alg_5_helper_v2(target, lc.ref, lc.constr, subseq, verbose, dangle);
                    if (result == "undesignable"){
                        std::cout<<"undesignable!"<<std::endl;
                        int count_pairs = lc.node->children.size() + 1;
                        std::string r = puzzle_id+","+y_star+","+std::to_string(count_pairs)+","+std::to_string(lc.node->first)+","+std::to_string(lc.node->second);
                        for(auto child: lc.node->children)
                            r += ","+std::to_string(child->first)+","+std::to_string(child->second);
                        r = r + ","+y_sub+","+std::to_string(y_rivals.size());
                        for(auto rival: y_rivals)
                            r += ","+rival;
                        std::cout<<r<<std::endl;
                        records.push_back(r);
                        outputFile << r << std::endl;
                        // break;
                        // return;
                    }
                    printf("\n");
                }
            }
            if (alg == "neighbor2"){
                std::vector<LoopComplex> lc_list;
                TreeNode* root = parseStringToTree(y_star);
                int max_internal = max_single(root);
                if(max_internal > 30){
                    std::string r = puzzle_id+","+std::to_string(max_internal);
                    records.push_back(r);
                    outputFile << r << std::endl;
                    continue;
                }
                tree2TwoNeighbor(root, y_star, lc_list);
                printf("lc_list size: %d\n", lc_list.size());

                // Sort the vector using a lambda expression
                std::sort(lc_list.begin(), lc_list.end(), [](const LoopComplex &a, const LoopComplex &b) {
                    return a.count_uk < b.count_uk;});
                for (auto lc: lc_list){
                    std::string target = y_star.substr(lc.start, lc.end-lc.start+1);
                    std::string subseq = seq.substr(lc.start, lc.end-lc.start+1);
                    printf(" count: %d\n", lc.count_uk);
                    printf("target: %s\n", target.c_str());
                    printf("   ref: %s\n", lc.ref.c_str());
                    printf("constr: %s\n", lc.constr.c_str());

                    std::string result = alg_5_helper_v2(target, lc.ref, lc.constr, subseq, verbose, dangle);
                    if (result == "undesignable"){
                        std::cout<<"undesignable!"<<std::endl;
                        int count_pairs = lc.node->children.size() + 1;
                        std::string r = puzzle_id+","+y_star+","+std::to_string(count_pairs)+","+std::to_string(lc.node->first)+","+std::to_string(lc.node->second);
                        for(auto child: lc.node->children)
                            r += ","+std::to_string(child->first)+","+std::to_string(child->second);
                        r = r + ","+y_sub+","+std::to_string(y_rivals.size());
                        for(auto rival: y_rivals)
                            r += ","+rival;
                        std::cout<<r<<std::endl;
                        records.push_back(r);
                        outputFile << r << std::endl;
                        // break;
                        // return;
                    }
                    printf("\n");
                }
            }
            if (alg == "csgen"){
                std::vector<LoopComplex> lc_list;
                TreeNode* root = parseStringToTree(y_star);
                int max_internal = max_single(root);
                if(max_internal > 30){
                    std::string r = puzzle_id+","+std::to_string(max_internal);
                    records.push_back(r);
                    outputFile << r << std::endl;
                    continue;
                }
                tree2Loops(root, y_star, lc_list);
                printf("lc_list size: %d\n", lc_list.size());

                // Sort the vector using a lambda expression
                std::sort(lc_list.begin(), lc_list.end(), [](const LoopComplex &a, const LoopComplex &b) {
                    return a.count_uk < b.count_uk;});
                for (auto lc: lc_list){
                    std::string target = y_star.substr(lc.start, lc.end-lc.start+1);
                    std::string subseq = seq.substr(lc.start, lc.end-lc.start+1);
                    printf(" count: %d\n", lc.count_uk);
                    printf("target: %s\n", target.c_str());
                    printf("   ref: %s\n", lc.ref.c_str());
                    printf("constr: %s\n", lc.constr.c_str());
                    outputFile << subseq << std::endl;
                    outputFile << lc.constr << std::endl;
                    printf("\n");
                }
            }
        }
    }
    for (auto r: records)
        std::cout<<r<<std::endl;
    
    // Close the file
    outputFile.close();
    std::cout << "Strings written to file: " << fileName << std::endl;
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
    cxxopts::Options options("MyProgram", "One line description of MyProgram");
    options.add_options()
    ("a,alg", "Algorithm", cxxopts::value<std::string>()->default_value("0"))
    ("c,csv", "csv file", cxxopts::value<std::string>()->default_value(""))
    ("t,txt", "txt file", cxxopts::value<std::string>()->default_value(""))
    ("d,dangle", "Dangle mode", cxxopts::value<int>()->default_value("2"))
    ("v,verbose", "Verbose output", cxxopts::value<bool>()->default_value("false"))
    ;

    auto result = options.parse(argc, argv);
    std::string alg = result["alg"].as<std::string>();
    std::string csv = result["csv"].as<std::string>();
    std::string txt = result["txt"].as<std::string>();
    verbose = result["verbose"].as<bool>();
    dangle = result["dangle"].as<int>();
    printf("alg: %s, verbose: %d, dangle: %d\n", alg.c_str(), verbose, dangle);
    show_configuration();

    if (alg == "0"){
        std::cout<<"no alg was selected!"<<std::endl;
        return 0;
    }

    if (!csv.empty()){
        csv_process(csv, alg);
        return 0;
    }

    if (!txt.empty()){
        txt_process(txt, alg);
        return 0;
    }


    if ( alg == "csfold" ){  /* constrained folding */
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
    }else if (alg == "fold"){ /* fold */
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
    }else if (alg == "3"){ /* alg 3 */
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
    }else if (alg == "-1"){ /* alg 1 (deprecated)*/
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
    }else if(alg == "1"){ /* alg 1 version 2 */
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
            // std::vector<std::string> X = alg_1_v2(ref1, ref2, seq, verbose, dangle);
            alg1_helper(seq, ref1, ref2, verbose, dangle);
            auto end_time = std::chrono::high_resolution_clock::now();
            const std::chrono::duration<double, std::milli> time_ms = end_time - start_time;
            // printf("X size: %d\n", X.size());
            // if (X.size()==0)
            //     printf("undesignable!\n");
            printf("alg1(v2) time: %.4f seconds\n", time_ms/1000.f);
        }
    }else if (alg == "2"){ /* alg 2 */
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
    }else if (alg == "2c"){ /* constrained alg2 */
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
    }else if (alg == "ed"){ /* energy difference: call the function test_diff  */
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
    }else if (alg == "dp"){  /* differential positions */ 
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
    }else if (alg == "eval"){ /* energy evaluation  */
        std::string seq;
        std::string ref;
        // Read input line by line until EOF (end of file) is reached
        while (std::getline(std::cin, seq)) {
            getline(std::cin, ref);
            long energy = linear_eval(seq, ref, verbose, dangle);
            printf("total energy: %.2f\n", energy/-100.0);
        }
    }else if (alg == "loop"){ /* loops evaluation  */
        std::string seq;
        std::string ref;
        // Read input line by line until EOF (end of file) is reached
        while (std::getline(std::cin, seq)) {
            getline(std::cin, ref);
            std::vector<LoopComplex> lc_list;
            TreeNode* root = parseStringToTree(ref);
            tree2Loops(root, ref, lc_list);
            printf("lc_list size: %d\n", lc_list.size());
            // Sort the vector using a lambda expression
            std::sort(lc_list.begin(), lc_list.end(), [](const LoopComplex &a, const LoopComplex &b) {
                return a.count_uk < b.count_uk;});
            for (auto lc: lc_list){
                std::string target = ref.substr(lc.start, lc.end-lc.start+1);
                std::string subseq = seq.substr(lc.start, lc.end-lc.start+1);
                printf(" count: %d\n", lc.count_uk);
                printf("target: %s\n", target.c_str());
                printf("   ref: %s\n", lc.ref.c_str());
                printf("constr: %s\n", lc.constr.c_str());

                if (true){
                    std::string result = alg_5_helper_v2(target, lc.ref, lc.constr, subseq, verbose, dangle);
                    if (result == "undesignable"){
                        printf("undesignable span: %d\t", lc.node->first);
                        for(auto child: lc.node->children){
                            printf("%d\t%d\t", child->first, child->second);
                        }
                        printf("%d\n", lc.node->second);
                        break;
                    }
                }

                printf("\n");
            }
        }
    }else if (alg == "edge"){ /* edges evaluation  */
        std::string seq;
        std::string ref;
        // Read input line by line until EOF (end of file) is reached
        while (std::getline(std::cin, seq)) {
            getline(std::cin, ref);
            std::vector<LoopComplex> lc_list;
            TreeNode* root = parseStringToTree(ref);
            tree2Edges(root, ref, lc_list);
            printf("lc_list size: %d\n", lc_list.size());
            // Sort the vector using a lambda expression
            std::sort(lc_list.begin(), lc_list.end(), [](const LoopComplex &a, const LoopComplex &b) {
                return a.count_uk < b.count_uk;});
            for (auto lc: lc_list){
                std::string target = ref.substr(lc.start, lc.end-lc.start+1);
                std::string subseq = seq.substr(lc.start, lc.end-lc.start+1);
                printf(" count: %d\n", lc.count_uk);
                printf("target: %s\n", target.c_str());
                printf("   ref: %s\n", lc.ref.c_str());
                printf("constr: %s\n", lc.constr.c_str());

                if (true){
                    std::string result = alg_5_helper_v2(target, lc.ref, lc.constr, subseq, verbose, dangle);
                    if (result == "undesignable")
                        break;
                }

                printf("\n");
            }
        }
    }else if (alg == "mloop"){ /* multi-loops evaluation  */
        std::string seq;
        std::string ref;
        // Read input line by line until EOF (end of file) is reached
        while (std::getline(std::cin, seq)) {
            getline(std::cin, ref);
            std::vector<LoopComplex> lc_list;
            TreeNode* root = parseStringToTree(ref);
            tree2MLoops(root, ref, lc_list);
            printf("lc_list size: %d\n", lc_list.size());
            // Sort the vector using a lambda expression
            std::sort(lc_list.begin(), lc_list.end(), [](const LoopComplex &a, const LoopComplex &b) {
                return a.count_uk < b.count_uk;});
            for (auto lc: lc_list){
                std::string target = ref.substr(lc.start, lc.end-lc.start+1);
                std::string subseq = seq.substr(lc.start, lc.end-lc.start+1);
                printf(" count: %d\n", lc.count_uk);
                printf("target: %s\n", target.c_str());
                printf("   ref: %s\n", lc.ref.c_str());
                printf("constr: %s\n", lc.constr.c_str());

                if (true){
                    std::string result = alg_5_helper_v2(target, lc.ref, lc.constr, subseq, verbose, dangle);
                    if (result == "undesignable"){
                        printf("undesignable span: %d\t", lc.node->first);
                        for(auto child: lc.node->children){
                            printf("%d\t%d\t", child->first, child->second);
                        }
                        printf("%d\n", lc.node->second);
                        break;
                    }
                }

                printf("\n");
            }
        }
    }else if (alg == "motif"){ /* motif evaluation  */
        std::string seq;
        std::string target;
        std::string ref;
        std::string cst;
        // Read input line by line until EOF (end of file) is reached
        while (std::getline(std::cin, seq)) {
            getline(std::cin, target);
            getline(std::cin, cst);
            // ref = cst;
            // std::replace(ref.begin(), ref.end(), '?', '.');
            // std::vector<std::string> mfes = cs_fold(seq, cst, 0, false, verbose, dangle);
            // if(isUMFE(mfes, target)){
            //     // std::count<<seq<<std::endl;
            //     std::cout<<seq<<std::endl;
            //     std::cout<<"designable! (umfe) "<<std::endl;
            //     // std::count<<"designable! (umfe)"<<std::endl;
            //     continue;
            // }else if(isMFE(mfes, target)){
            //     if(mfes[0] == target)
            //         ref = mfes[1];
            //     else
            //         ref = mfes[0];
            // }else{
            //     ref = mfes[0];
            // }
            
            // std::vector<LoopComplex> lc_list;
            // TreeNode* root = parseStringToTree(ref);
            // tree2Loops(root, ref, lc_list);
            // printf("lc_list size: %d\n", lc_list.size());
            // // Sort the vector using a lambda expression
            // std::sort(lc_list.begin(), lc_list.end(), [](const LoopComplex &a, const LoopComplex &b) {
            //     return a.count_uk < b.count_uk;});
            // for (auto lc: lc_list){
                // std::string target = ref.substr(lc.start, lc.end-lc.start+1);
                // std::string subseq = seq.substr(lc.start, lc.end-lc.start+1);
                printf(" count: %d\n", countOccurrences(cst, '?'));
                printf("target: %s\n", target.c_str());
                printf("constr: %s\n", cst.c_str());
                // printf("   ref: %s\n", ref.c_str());
                if (true){
                    std::string result = alg_5_helper_v2(target, ref, cst, seq, verbose, dangle);
                    std::cout<<"result: "<<result<<std::endl;
                    // if (result == "undesignable"){
                    //     // printf("undesignable span: %d\t", lc.node->first);
                    //     for(auto child: lc.node->children){
                    //         printf("%d\t%d\t", child->first, child->second);
                    //     }
                    //     printf("%d\n", lc.node->second);
                    //     break;
                    // }
                }
                printf("\n");
            // }
        }
    }else if (alg == "neighbor2"){ /* edges evaluation  */
        std::string seq;
        std::string ref;
        // Read input line by line until EOF (end of file) is reached
        while (std::getline(std::cin, seq)) {
            getline(std::cin, ref);
            std::vector<LoopComplex> lc_list;
            TreeNode* root = parseStringToTree(ref);
            tree2TwoNeighbor(root, ref, lc_list);
            printf("lc_list size: %d\n", lc_list.size());
            // Sort the vector using a lambda expression
            std::sort(lc_list.begin(), lc_list.end(), [](const LoopComplex &a, const LoopComplex &b) {
                return a.count_uk < b.count_uk;});
            for (auto lc: lc_list){
                std::string target = ref.substr(lc.start, lc.end-lc.start+1);
                std::string subseq = seq.substr(lc.start, lc.end-lc.start+1);
                printf(" count: %d\n", lc.count_uk);
                printf("target: %s\n", target.c_str());
                printf("   ref: %s\n", lc.ref.c_str());
                printf("constr: %s\n", lc.constr.c_str());
                auto start_time = std::chrono::high_resolution_clock::now();
                std::string result = alg_5_helper_v2(target, lc.ref, lc.constr, subseq, verbose, dangle);
                if (result == "undesignable"){
                    std::cout<<"undesignable!"<<std::endl;
                    auto end_time = std::chrono::high_resolution_clock::now();
                    const std::chrono::duration<double, std::milli> time_ms = end_time - start_time;
                    printf("time cost: %.4f seconds\n", time_ms/1000.f);
                    int count_pairs = lc.node->children.size() + 1;
                    std::string r = seq+","+ref+","+std::to_string(count_pairs)+","+std::to_string(lc.node->first)+","+std::to_string(lc.node->second);
                    for(auto child: lc.node->children)
                        r += ","+std::to_string(child->first)+","+std::to_string(child->second);
                    r = r + ","+y_sub+","+std::to_string(y_rivals.size());
                    for(auto rival: y_rivals)
                        r += ","+rival;
                    // Convert duration to seconds and then cast to float
                    float time_seconds = std::chrono::duration_cast<std::chrono::duration<float>>(time_ms).count();
                    r += " time:" + fl2str(time_seconds) + ";";
                    std::cout<<r<<std::endl;
                    // std::string id = puzzle_id + "_" + alg;
                    std::string args4plot = compose_args4plot("id", ref, lc.ps_outside, lc.ps_inside);
                    std::cout << args4plot <<std::endl;
                    std::string pairstring = compose_pairstr(lc.ps_inside, lc.ps_outside);
                    std::cout << pairstring << std::endl;
                }
                printf("\n");
            }
        }
    }else if (alg == "showtree"){
        std::string ref;
        // Read input line by line until EOF (end of file) is reached
        while (std::getline(std::cin, ref)) {
            TreeNode* root = parseStringToTree(ref);
            root->printTree();
        }
    }
    return 0;
}