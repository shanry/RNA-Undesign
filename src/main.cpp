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
#include <cassert> 

// #include "eval.cpp"
#include "cxxopts.hpp"
#include "utils.h"
#include "comps.h"
#include "eval.h"
using namespace std;

#define MAX_ENUM 10000000000
#define MAX_CONSTRAINT 100000
#define MAX_SEQ 500


/* Old compatibility names for C types.  */
// typedef unsigned long int ulong;
// typedef unsigned short int ushort;
// typedef unsigned int uint;

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

int SEED_RAND = 2024;

// std::string PATH_DESIGNABLE_LIB = "lib_designable.txt";
// std::string PATH_UNDESIGNABLE_LIB = "lib_undesignable.txt";

// load lib from file
std::set<std::string> loadLib(std::string path){
    std::set<std::string> lib;
    std::ifstream file(path);
    std::string line;
    while (std::getline(file, line)){
        json js_record = json::parse(line);
        assert(!js_record["is_duplicated"]);
        Node* tree = new Node(js_record["motif"]);
        std::string treestr = tree->toDotBracket();
        lib.insert(treestr);
        for(Node* rotree : tree->rotated(0)){
            std::string rotreestr = rotree->toDotBracket();
            lib.insert(rotreestr);
            delete rotree;
        }
        delete tree;
    }
    return lib;
}

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
            node->child_id = nodeStack.top()->children.size(); // 1-based child id
            nodeStack.push(node);
        } else if (c == ')') {
            nodeStack.top()->second = i;
            nodeStack.top()->setLoop();
            nodeStack.pop();

        }
    }
    assert(nodeStack.size() == 1);
    nodeStack.top()->setLoop();
    nodeStack.pop();

    return root;
}

std::string getCurrentTimestamp() {
    // Get the current time point
    auto currentTime = std::chrono::system_clock::now();

    // Convert the time point to a time_t object
    std::time_t currentTime_t = std::chrono::system_clock::to_time_t(currentTime);

    // Convert the time_t to a std::tm structure
    std::tm* timeInfo = std::localtime(&currentTime_t);
    
    // Use strftime to format the timestamp
    char buffer[15]; // Adjust the size as needed based on your format
    std::strftime(buffer, sizeof(buffer), "%Y%m%d%H%M%S", timeInfo);

    // Convert the stringstream to a string
    std::ostringstream timestampStream;
    timestampStream << buffer;

    // Return the formatted timestamp as a string
    return timestampStream.str();
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
    std::cout << "    x: " << seq << std::endl;
    std::cout << "ystar: " << y << std::endl;
    std::cout << "yprim: " << y_prime << std::endl;
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
        printf("alg1.X size: %d\n", X.size());
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
    printf("alg1.X size: %d\n", X.size());
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
    std::cout<<"inside alg5cs"<<std::endl;
    int count_cs = cs_vec.size();
    printf("design constraints: %d\n", count_cs);
    for(int idx_cs; idx_cs < count_cs; idx_cs++){
        std::cout<<cs_vec[idx_cs].structure<<std::endl;
        std::cout<<cs_vec[idx_cs].seqs->size()<<std::endl;
    }
    std::vector<std::pair<ulong, std::pair<std::string, std::string>>> y_primes;
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
            std::cout<<"UMFE log: "<<ref1<<std::endl;
            std::cout<<"UMFE log: "<<constr<<std::endl;
            std::cout<<"UMFE log: "<<x<<std::endl;
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

std::string alg_5_cs_plus(std::string& ref1, std::set<std::string>& refs_checked, std::vector<Constraint>& cs_vec, std::vector<std::string> X, std::string& constr, bool verbose, int dangle_model){ // ref1, ref2, X, is_verbose, dangle_model
    std::cout<<"inside alg5cs"<<std::endl;
    int count_cs = cs_vec.size();
    printf("design constraints: %d\n", count_cs);
    for(int idx_cs; idx_cs < count_cs; idx_cs++){
        std::cout<<cs_vec[idx_cs].structure<<std::endl;
        std::cout<<cs_vec[idx_cs].seqs->size()<<std::endl;
    }
    std::vector<std::pair<ulong, std::pair<std::string, std::string>>> y_primes;
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
            std::cout<<"UMFE log: "<<ref1<<std::endl;
            std::cout<<"UMFE log: "<<constr<<std::endl;
            std::cout<<"UMFE log: "<<x<<std::endl;
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

std::string alg_5_helper_v2(std::string& ref1, std::string& ref2, std::string&constr, std::string& seq, bool verbose, int dangle_model){
    std::cout<< "inside alg_5_helper_v2"<<std::endl;
    std::cout << "seq: " << seq << std::endl;
    std::cout << "  y: " << ref1 << std::endl;
    seq_init = seq;
    std::vector<std::string> mfes = cs_fold(seq, constr, 0, false, verbose, dangle);
    std::cout<<"mfes.size: "<<mfes.size()<<std::endl;
    std::string ref_mfe;
    if(isUMFE(mfes, ref1)){
        std::cout<<seq<<std::endl;
        std::cout<<"designable!"<<std::endl;
        std::cout<<"UMFE log: "<<ref1<<std::endl;
        std::cout<<"UMFE log: "<<constr<<std::endl;
        std::cout<<"UMFE log: "<<seq<<std::endl;
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
    std::cout<<"critical positions: ";
    for(int critical_position: critical_positions)
        std::cout<<critical_position<<"\t";
    std::cout<<std::endl;
    // long delta_energy = diff_eval(seq, cr_loops, verbose, dangle_model);
    std::vector<std::tuple<int, int>> pairs_diff = idx2pair(critical_positions, ref1);
    std::cout<<"pairs_diff: ";
    for (auto pair: pairs_diff)
        std::cout<<std::get<0>(pair)<<"\t"<<std::get<1>(pair)<<";\t";
    std::cout<<std::endl;
    ulong n_enum = count_enum(pairs_diff);
    std::cout<<"enumeration count: "<<n_enum<<std::endl;
    // check the ref from removing internal pairs
    std::vector<std::string> X;
    bool cs_flag = false;
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
        }
        else{            
            cs_flag = true;
            // cs_vec.push_back(cs_ref2);
            Constraint cs_ref2 = Constraint(&critical_positions, &X);
            cs_ref2.setStructure(ref2);
            cs_vec.push_back(cs_ref2);
            refs_checked.insert(ref2);
        }
    }
    std::vector<std::string> Xseeds;
    Xseeds.push_back(seq);
    for(int i=0; i<10; i++)
        Xseeds.push_back(tg_init(ref1)); // can be deduplicated here
    return alg_5_cs_plus(ref1, refs_checked, cs_vec, Xseeds, constr, verbose, dangle_model);
    

    // finish the first check
    // if(ref_mfe == ref2){
    //     if(cs_flag){
    //         Constraint cs_ref2 = Constraint(&critical_positions, &X);
    //         cs_ref2.setStructure(ref2);
    //         cs_vec.push_back(cs_ref2);
    //         refs_checked.insert(ref2);
    //     }
    //     if(cs_vec.size())
    //         return alg_5_cs(ref1, refs_checked, cs_vec, constr, verbose, dangle_model);
    //     else{
    //         std::cout<<"intial y' has too many constraints!"<<std::endl;
    //         assert (X.size() > 0);
    //         return alg_5_cs_plus(ref1, refs_checked, cs_vec, X, constr, verbose, dangle_model);
    //     }
    // }

    // check the ref from mfe
    // std::set<int> critical_positions_mfe;
    // cr_loops = find_critical_plus(ref1, ref_mfe, critical_positions_mfe, verbose);
    // // delta_energy = diff_eval(seq, cr_loops, verbose, dangle_model);
    // pairs_diff = idx2pair(critical_positions_mfe, ref1);
    // n_enum = count_enum(pairs_diff);
    // std::cout<<"enumeration count: "<<n_enum<<std::endl;
    // std::vector<std::string> X_mfe;
    // if(n_enum > 0 && n_enum < MAX_ENUM){
    //     std::cout<<"alg 1"<<std::endl;
    //     // auto X = alg_1(ref1, ref2, cr_loops, pairs_diff, seq, verbose, dangle_model);
    //     X_mfe = alg_1_v2(ref1, ref_mfe, seq, verbose, dangle_model);
    //     printf("X_mfe size: %d\n", X_mfe.size());
    //     // std::vector<Constraint> cs_vec;
    //     if (X_mfe.size() == 0){
    //         std::cout<<"undesignable!"<<std::endl;
    //         y_sub = ref1;
    //         y_rivals.clear();
    //         y_rivals.push_back(ref_mfe);
    //         return "undesignable";
    //     }else if (X_mfe.size() > MAX_CONSTRAINT){
    //         std::cout<<"too many constraints: "<<X_mfe.size()<<"\t"<<"out of "<<ref_mfe<<std::endl;
    //         if(cs_flag){
    //             Constraint cs_ref2 = Constraint(&critical_positions, &X);
    //             cs_ref2.setStructure(ref2);
    //             cs_vec.push_back(cs_ref2);
    //             refs_checked.insert(ref2);
    //         }
    //     }else{
    //         Constraint cs_ref_mfe = Constraint(&critical_positions_mfe, &X_mfe);
    //         cs_ref_mfe.setStructure(ref_mfe);
    //         cs_vec.push_back(cs_ref_mfe);
    //         refs_checked.insert(ref_mfe);
    //     }
    //     if(cs_vec.size())
    //         return alg_5_cs(ref1, refs_checked, cs_vec, constr, verbose, dangle_model);
    //     else{
    //         std::cout<<"both intial ys have too many constraints!"<<std::endl;
    //         assert (X_mfe.size() > 0);
    //         return alg_5_cs_plus(ref1, refs_checked, cs_vec, X_mfe, constr, verbose, dangle_model);
    //     }
    // }
    // std::cout<<"intial y_prime too bad!"<<std::endl;
    // return "intial y_prime too bad";
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
    std::vector<std::string> records_designable;
    int count_designable = 0;
    // Specify the file name
    std::string fileName = csv + "." + alg + ".log."+getCurrentTimestamp()+".txt";
    #ifdef SPECIAL_HP
    #else
        fileName = csv + "." + alg + ".log.nosh."+getCurrentTimestamp()+".txt";
    #endif

    // Output file stream
    std::ofstream outputFile(fileName);
    // Check if the file is open
    if (!outputFile.is_open()) {
        std::cerr << "Error opening the file: " << fileName << std::endl;
        return;
    }
    // Library files for designable and undesignable motifs
    std::ofstream designableLibFile("lib_designable.txt", std::ios::app);
    std::ofstream undesignableLibFile("lib_undesignable.txt", std::ios::app);
    if (!designableLibFile.is_open()) {
        std::cerr << "Error opening the file: " << "lib_designable.txt" << std::endl;
        return;
    }
    if (!undesignableLibFile.is_open()) {
        std::cerr << "Error opening the file: " << "lib_undesignable.txt" << std::endl;
        return;
    }

    std::string fileTime = csv + "." + alg + ".time."+getCurrentTimestamp()+".csv";
    #ifdef SPECIAL_HP
    #else
        fileTime = csv + "." + alg + ".time.nosh."+getCurrentTimestamp()+".csv";
    #endif
    // Open the file for writing
    std::ofstream timeFile(fileTime);

    // Check if the file is open
    if (!timeFile.is_open()) {
        std::cerr << "Error opening the file: " << fileTime << std::endl;
        return;
    }else{
        timeFile << "ID,Time(s)" << std::endl;
    }
    std::unordered_map<std::string, GroupY> constr2groupy;
    // std::unordered_map<std::string, std::string> uniq_ud;
    const char*  var_undesignable_lib = std::getenv("PATH_UNDESIGNABLE_LIB");
    const char*  var_designable_lib = std::getenv("PATH_DESIGNABLE_LIB");
    if(var_undesignable_lib == NULL){
        std::cerr << "Error: PATH_UNDESIGNABLE_LIB is not set" << std::endl;
        return;
    }
    if(var_designable_lib == NULL){
        std::cerr << "Error: PATH_DESIGNABLE_LIB is not set" << std::endl;
        return;
    }
    std::string path_undesignable_lib(var_undesignable_lib);
    std::string path_designable_lib(var_designable_lib);
    std::set<std::string> uniq_ud = loadLib(path_undesignable_lib); // undesignable motifs
    std::set<std::string> uniq_ds = loadLib(path_designable_lib);   // designable   motifs
    for(int i = 1; i < df.size(); i++){
        auto row = df[i];
        {
            string sharph0 = "()";
            string sharph1 = "(.)";
            string sharph2 = "(..)";

            std::string puzzle_id = row[0];
            std::string seq;
            std::string y_star = row[1];
            if(y_star.find(sharph0) != std::string::npos || y_star.find(sharph1) != std::string::npos || y_star.find(sharph2) != std::string::npos){
                std::cout<<"sharp turn!"<<std::endl;
                std::cout<<y_star<<std::endl;
                continue;
            }
            if(row.size() > 2)
                seq = row[2];
            else
                seq = tg_init(y_star);
            std::string y_prim; // currently dummy

            std::cout<<"Puzzle ID: "<<puzzle_id<<std::endl;

            if (alg == "1"){
                auto start_time = std::chrono::high_resolution_clock::now();
                std::string result = alg1_helper(seq, y_star, y_prim, verbose, dangle);
                auto end_time = std::chrono::high_resolution_clock::now();
                const std::chrono::duration<double, std::milli> time_ms = end_time - start_time;
                printf("alg1(v2) time: %.4f seconds\n", time_ms/1000.f);
                if (result == "undesignable"){
                    std::string r = puzzle_id+","+y_star+","+y_prim;
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
                    std::string r = puzzle_id+","+y_star+","+std::to_string(y_rivals.size());
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

                    std::string r = puzzle_id+","+y_star+","+std::to_string(found)+","+std::to_string(found+y_sub.length())+","+y_sub+","+std::to_string(y_rivals.size());
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

                                std::string r = puzzle_id+","+y_star+","+std::to_string(found)+","+std::to_string(found+y_sub.length())+","+y_sub+","+std::to_string(y_rivals.size());
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
                // int max_internal = max_single(root);
                // if(max_internal > 30){
                //     std::string r = puzzle_id+","+std::to_string(max_internal);
                //     records.push_back(r);
                //     outputFile << r << std::endl;
                //     continue;
                // }
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
                    std::string result;
                    try{
                        result = alg_5_helper_v2(target, lc.ref, lc.constr, subseq, verbose, dangle);
                    }catch(...){
                        // Code to handle any exception
                        std::cerr << "An exception occurred" << std::endl;
                        result = "exception";
                    }
                    if (result == "exception")
                        continue;
                    if (result == "undesignable"){
                        std::cout<<"undesignable!"<<std::endl;
                        auto end_time = std::chrono::high_resolution_clock::now();
                        const std::chrono::duration<double, std::milli> time_ms = end_time - start_time;
                        float time_seconds = std::chrono::duration_cast<std::chrono::duration<float>>(time_ms).count();
                        printf("time cost: %.4f seconds\n", time_ms/1000.f);
                        auto js = jsrecords(lc, y_star, y_sub, y_rivals, puzzle_id);
                        js["time"] = time_seconds;
                        std::string jstring = js.dump();
                        outputFile << jstring << std::endl;
                        records.push_back(jstring);
                        // size_t found = y_star.find(y_sub);
                        // assert (found != std::string::npos);
                        // int found_end = found+y_sub.length();
                        // std::string r = puzzle_id+","+y_star+",1,"+std::to_string(lc.node->first)+","+std::to_string(lc.node->second)+","+y_sub+","+std::to_string(y_rivals.size());
                        // for(auto rival: y_rivals)
                        //     r += ","+rival;
                        // // Convert duration to seconds and then cast to float
                        // r += " time:" + fl2str(time_seconds, 4) + ";";
                        // std::cout<<r<<std::endl;
                        // records.push_back(r);
                        // std::string id = puzzle_id + "_" + alg;
                        // std::string args4plot = compose_args4plot(id, y_star, lc.ps_outside, lc.ps_inside);
                        // outputFile << r << std::endl;
                        // outputFile << args4plot << std::endl;
                        // std::string pairstring = id + ":" + compose_pairstr(lc.ps_inside, lc.ps_outside);
                        // std::cout << pairstring << std::endl;
                        // outputFile << pairstring << std::endl;
                        // std::string jstr = lc.jsmotif(puzzle_id);
                        // std::cout << jstr << std::endl;
                        // outputFile << jstr << std::endl;
                        // break;
                    }
                    printf("\n");
                }
            }
            // std::string goal_test = "p [] (M [0, 5, 0] (p [] (), B [1, 0] (S [0, 0] (p [] ()))))";
            // power neighbor set search
            if (alg == "pn"){
                auto start_time = std::chrono::high_resolution_clock::now();
                TreeNode* root = parseStringToTree(y_star);
                std::set<string> ds_ipairs; // designable internal pairs
                std::set<string> ud_ipairs; // undesinable internal pairs

                // motif of 2 loops (1 edge/internal pair)
                std::vector<LoopComplex> lc_list;
                tree2Edges(root, y_star, lc_list);
                printf("lc_list size: %d\n", lc_list.size());
                // sort motifs by the number of unknown positions
                std::sort(lc_list.begin(), lc_list.end(), [](const LoopComplex &a, const LoopComplex &b) {
                    return a.count_uk < b.count_uk;});

                // motif of 3 loops (2 edges/internal pairs)
                std::vector<LoopComplex> lc_list_2nbs;
                tree2TwoNeighbor(root, y_star, lc_list_2nbs);
                printf("lc_list_2nbs size: %d\n", lc_list_2nbs.size());
                // sort motifs by the number of unknown positions
                std::sort(lc_list_2nbs.begin(), lc_list_2nbs.end(), [](const LoopComplex &a, const LoopComplex &b) {
                    return a.count_uk < b.count_uk;});
                lc_list.insert(lc_list.end(), lc_list_2nbs.begin(), lc_list_2nbs.end());

                // motif of 4 loops (3 edges/internal pairs)
                std::vector<LoopComplex> lc_list_3nbs;
                tree2ThreeNeighbor(root, y_star, lc_list_3nbs);
                printf("lc_list_3nbs size: %d\n", lc_list_3nbs.size());
                // sort motifs by the number of unknown positions
                std::sort(lc_list_3nbs.begin(), lc_list_3nbs.end(), [](const LoopComplex &a, const LoopComplex &b) {
                    return a.count_uk < b.count_uk;});
                lc_list.insert(lc_list.end(), lc_list_3nbs.begin(), lc_list_3nbs.end());

                for (auto lc: lc_list){
                    auto start_time_lc = std::chrono::high_resolution_clock::now();
                    lc.printLoopLens();
                    if (lc.hasLongLoop()){
                        printf("the loop exceeds length limit!");
                        continue;
                    }
                    std::string target = y_star.substr(lc.start, lc.end-lc.start+1);
                    std::string subseq = seq.substr(lc.start, lc.end-lc.start+1);
                    printf(" count: %d\n", lc.count_uk);
                    printf("target: %s\n", target.c_str());
                    printf("   ref: %s\n", lc.ref.c_str());
                    printf("constr: %s\n", lc.constr.c_str());
                    auto ipairs_subsets = pairSubSet(lc.ps_inside);
                    std::string result;
                    // check if the (rotated) motif is already found undesignable
                    json js_motif = json::parse(lc.jsmotif(puzzle_id));
                    std::cout<<"js_motif: "<<js_motif<<std::endl;
                    Node* tree = new Node(js_motif);
                    std::string treestr = tree->toDotBracket();
                    // if(!goal_test.empty() && treestr != goal_test){
                    //     delete tree;
                    //     continue;
                    // }
                    std::cout<<"treestr: "<<treestr<<std::endl;
                    if(uniq_ds.find(treestr) != uniq_ds.end()){
                        std::cout<<"already designable!"<<std::endl;
                        result = "designable";
                    }else if(uniq_ud.find(treestr) != uniq_ud.end()){
                        result = "undesignable";
                        std::cout<<"recur lc.constr: "<<lc.constr<<std::endl;
                        // std::cout<<"recur    groupy: "<<constr2groupy[lc.constr].constr<<std::endl;
                        std::cout<<"recur lc.constr: "<<target<<std::endl;
                        std::cout<<"recur   treestr: "<<treestr<<std::endl;
                        // std::cout<<"recur     ystar: "<<constr2groupy[lc.constr].star<<std::endl;
                    }else{
                        std::string ref_lc = lc.ref;
                        std::string constr_lc = lc.constr;
                        for( int ib = 1; ib < lc.ps_outside.size(); ib++){ // skip the most outside boudary pair
                            auto bpair = lc.ps_outside[ib];
                            std::cout<<"bpair: "<<bpair.first<<"\t"<<bpair.second<<std::endl;
                            int len_branch =  bpair.second - bpair.first + 1;
                            std::string y_branch = y_star.substr(bpair.first, len_branch);
                            std::cout<<"y_branch:"<<std::endl;
                            std::cout<<y_branch<<std::endl;
                            TreeNode* root_branch = parseStringToTree(y_branch);
                            if(max_single(root_branch) > SINGLE_MAX_LEN || max_multi(root_branch) > MULTIPLE_FIRST_MAX_LEN){
                                std::string helix_branch = genHelix(len_branch);
                                std::string seq_branch = tg_init(helix_branch);
                                target.replace(bpair.first - lc.start, len_branch, helix_branch);
                                ref_lc.replace(bpair.first - lc.start, len_branch, helix_branch);
                                constr_lc.replace(bpair.first - lc.start, len_branch, helix_branch);
                                subseq.replace(bpair.first - lc.start, len_branch, seq_branch);
                            }
                        }
                        printf("target: %s\n", target.c_str());
                        printf("   ref: %s\n", ref_lc.c_str());
                        printf("constr: %s\n", constr_lc.c_str());
                        bool ud = false; // check if the motif is already found undesignable
                        for(auto ipairs: ipairs_subsets){
                            if(ud_ipairs.find(pairs2string(ipairs))!=ud_ipairs.end()){
                                ud = true;
                                break;
                            }
                        }
                        if(ud)
                            continue;
                        try{
                            std::cout<<"UMFE log: before alg 5 "<<std::endl;
                            result = alg_5_helper_v2(target, ref_lc, constr_lc, subseq, verbose, dangle);
                            if (result == "UMFE"){
                                result = "designable";
                                std::cout<<"UMFE log: afer alg 5 "<<std::endl;
                                std::cout<<"UMFE log: "<<tree->toDotBracket()<<std::endl;
                                std::cout<<"UMFE log: "<<target<<std::endl;
                                std::cout<<"UMFE log: "<<constr_lc<<std::endl;
                            }
                        }catch(...){
                            // Code to handle any exception
                            std::cerr << "An exception occurred" << std::endl;
                            result = "exception";
                            continue;
                        }
                    }
                    if (result == "designable"){
                        std::cout<<"designable!"<<std::endl;
                        ds_ipairs.insert(pairs2string(lc.ps_inside));
                        count_designable++;
                        y_sub = y_star; // set y_sub as y_star
                        y_rivals.clear(); // clear y_rivals
                        bool found_ds = false;
                        if(uniq_ds.find(treestr) != uniq_ds.end()){
                            found_ds = true;
                            std::cout<<"already found designable!"<<std::endl;
                        }else{
                            uniq_ds.insert(treestr);
                            for(Node* rotree: tree->rotated(0)){
                                std::string rotreestr = rotree->toDotBracket();
                                uniq_ds.insert(rotreestr);
                                delete rotree;
                            }
                            auto end_time_lc = std::chrono::high_resolution_clock::now();
                            const std::chrono::duration<double, std::milli> time_ms = end_time_lc - start_time_lc;
                            float time_seconds = std::chrono::duration_cast<std::chrono::duration<float>>(time_ms).count();
                            printf("time cost: %.4f seconds\n", time_seconds);
                            auto js = jsrecords(lc, y_star, y_sub, y_rivals, puzzle_id);
                            js["time"] = time_seconds;
                            js["seed"] = SEED_RAND;
                            js["is_duplicated"] = found_ds;
                            std::string jstring = js.dump();
                            designableLibFile << jstring << std::endl;
                            records_designable.push_back(jstring);
                        }
                    }
                    if (result == "undesignable"){
                        std::cout<<"undesignable!"<<std::endl;
                        bool found_ud = false; // check if the motif is already found undesignable
                        // if(constr2groupy.find(lc.constr) != constr2groupy.end() && constr2groupy[lc.constr].star == target)
                        if(uniq_ud.find(treestr) != uniq_ud.end()){
                            found_ud = true;
                        }else{
                            uniq_ud.insert(treestr);
                            for(Node* rotree: tree->rotated(0)){
                                std::string rotreestr = rotree->toDotBracket();
                                uniq_ud.insert(rotreestr);
                                delete rotree;
                            }
                        }
                        auto end_time_lc = std::chrono::high_resolution_clock::now();
                        const std::chrono::duration<double, std::milli> time_ms = end_time_lc - start_time_lc;
                        float time_seconds = std::chrono::duration_cast<std::chrono::duration<float>>(time_ms).count();
                        printf("time cost: %.4f seconds\n", time_seconds);
                        auto js = jsrecords(lc, y_star, y_sub, y_rivals, puzzle_id);
                        js["time"] = time_seconds;
                        js["seed"] = SEED_RAND;
                        js["is_duplicated"] = found_ud;
                        ud_ipairs.insert(pairs2string(lc.ps_inside));
                        std::vector<std::vector<std::pair<int, int>>> uk_pairs;
                        for(auto ipairs: ipairs_subsets){
                            if(ds_ipairs.find(pairs2string(ipairs))==ds_ipairs.end()){
                                uk_pairs.push_back(ipairs);
                            }
                        }
                        if(uk_pairs.size()){
                            js["ismin"] = false;
                            js["uk_ipairs"] = uk_pairs;
                        }else{
                            js["ismin"] = true;
                        }
                        std::string jstring = js.dump();
                        outputFile << jstring << std::endl;
                        records.push_back(jstring);
                        if (!found_ud){
                            undesignableLibFile << jstring << std::endl;
                        }
                    }
                    for(auto pair: ds_ipairs)
                        std::cout<<pair<<"  ";
                    delete tree;
                    printf("\n");
                }
                auto end_time = std::chrono::high_resolution_clock::now();
                const std::chrono::duration<double, std::milli> time_ms_y = end_time - start_time;
                float time_seconds_lc = std::chrono::duration_cast<std::chrono::duration<float>>(time_ms_y).count();
                printf("time cost for whole structure: %.4f seconds\n", time_seconds_lc);
                timeFile << puzzle_id << "," << time_seconds_lc <<std::endl;
            }
            if (alg == "scan"){
                auto start_time = std::chrono::high_resolution_clock::now();
                TreeNode* root = parseStringToTree(y_star);
                std::set<string> ds_ipairs; // designable internal pairs
                std::set<string> ud_ipairs; // undesinable internal pairs

                // motif of 2 loops (1 edge/internal pair)
                std::vector<LoopComplex> lc_list;
                tree2Edges(root, y_star, lc_list);
                printf("lc_list size: %d\n", lc_list.size());
                // sort motifs by the number of unknown positions
                std::sort(lc_list.begin(), lc_list.end(), [](const LoopComplex &a, const LoopComplex &b) {
                    return a.count_uk < b.count_uk;});

                // motif of 3 loops (2 edges/internal pairs)
                std::vector<LoopComplex> lc_list_2nbs;
                tree2TwoNeighbor(root, y_star, lc_list_2nbs);
                printf("lc_list_2nbs size: %d\n", lc_list_2nbs.size());
                // sort motifs by the number of unknown positions
                std::sort(lc_list_2nbs.begin(), lc_list_2nbs.end(), [](const LoopComplex &a, const LoopComplex &b) {
                    return a.count_uk < b.count_uk;});
                lc_list.insert(lc_list.end(), lc_list_2nbs.begin(), lc_list_2nbs.end());

                // motif of 4 loops (3 edges/internal pairs)
                std::vector<LoopComplex> lc_list_3nbs;
                tree2ThreeNeighbor(root, y_star, lc_list_3nbs);
                printf("lc_list_3nbs size: %d\n", lc_list_3nbs.size());
                // sort motifs by the number of unknown positions
                std::sort(lc_list_3nbs.begin(), lc_list_3nbs.end(), [](const LoopComplex &a, const LoopComplex &b) {
                    return a.count_uk < b.count_uk;});
                lc_list.insert(lc_list.end(), lc_list_3nbs.begin(), lc_list_3nbs.end());

                for (auto lc: lc_list){
                    auto start_time_lc = std::chrono::high_resolution_clock::now();
                    lc.printLoopLens();
                    if (lc.hasLongLoop()){
                        printf("the loop exceeds length limit!");
                        continue;
                    }
                    std::string target = y_star.substr(lc.start, lc.end-lc.start+1);
                    std::string subseq = seq.substr(lc.start, lc.end-lc.start+1);
                    printf(" count: %d\n", lc.count_uk);
                    printf("target: %s\n", target.c_str());
                    printf("   ref: %s\n", lc.ref.c_str());
                    printf("constr: %s\n", lc.constr.c_str());
                    auto ipairs_subsets = pairSubSet(lc.ps_inside);
                    std::string result;
                    // check if the (rotated) motif is already found undesignable
                    json js_motif = json::parse(lc.jsmotif(puzzle_id));
                    std::cout<<"js_motif: "<<js_motif<<std::endl;
                    Node* tree = new Node(js_motif);
                    std::string treestr = tree->toDotBracket();
                    // if(!goal_test.empty() && treestr != goal_test){
                    //     delete tree;
                    //     continue;
                    // }
                    std::cout<<"treestr: "<<treestr<<std::endl;
                    if(uniq_ds.find(treestr) != uniq_ds.end()){
                        std::cout<<"already designable!"<<std::endl;
                        result = "designable";
                    }else if(uniq_ud.find(treestr) != uniq_ud.end()){
                        result = "undesignable";
                        std::cout<<"recur lc.constr: "<<lc.constr<<std::endl;
                        // std::cout<<"recur    groupy: "<<constr2groupy[lc.constr].constr<<std::endl;
                        std::cout<<"recur lc.constr: "<<target<<std::endl;
                        std::cout<<"recur   treestr: "<<treestr<<std::endl;
                        // std::cout<<"recur     ystar: "<<constr2groupy[lc.constr].star<<std::endl;
                    }else {
                        result = "unknown";
                    }
                    if (result == "designable"){
                        std::cout<<"designable!"<<std::endl;
                        ds_ipairs.insert(pairs2string(lc.ps_inside));
                        count_designable++;
                        y_sub = y_star; // set y_sub as y_star
                        y_rivals.clear(); // clear y_rivals
                        bool found_ds = false;
                        if(uniq_ds.find(treestr) != uniq_ds.end()){
                            found_ds = true;
                            std::cout<<"already found designable!"<<std::endl;
                        }else{
                            uniq_ds.insert(treestr);
                            for(Node* rotree: tree->rotated(0)){
                                std::string rotreestr = rotree->toDotBracket();
                                uniq_ds.insert(rotreestr);
                                delete rotree;
                            }
                            auto end_time_lc = std::chrono::high_resolution_clock::now();
                            const std::chrono::duration<double, std::milli> time_ms = end_time_lc - start_time_lc;
                            float time_seconds = std::chrono::duration_cast<std::chrono::duration<float>>(time_ms).count();
                            printf("time cost: %.4f seconds\n", time_seconds);
                            auto js = jsrecords(lc, y_star, y_sub, y_rivals, puzzle_id);
                            js["time"] = time_seconds;
                            js["seed"] = SEED_RAND;
                            js["is_duplicated"] = found_ds;
                            std::string jstring = js.dump();
                            designableLibFile << jstring << std::endl;
                            records_designable.push_back(jstring);
                        }
                    }
                    if (result == "undesignable"){
                        std::cout<<"undesignable!"<<std::endl;
                        bool found_ud = false; // check if the motif is already found undesignable
                        // if(constr2groupy.find(lc.constr) != constr2groupy.end() && constr2groupy[lc.constr].star == target)
                        if(uniq_ud.find(treestr) != uniq_ud.end()){
                            found_ud = true;
                        }else{
                            uniq_ud.insert(treestr);
                            for(Node* rotree: tree->rotated(0)){
                                std::string rotreestr = rotree->toDotBracket();
                                uniq_ud.insert(rotreestr);
                                delete rotree;
                            }
                        }
                        auto end_time_lc = std::chrono::high_resolution_clock::now();
                        const std::chrono::duration<double, std::milli> time_ms = end_time_lc - start_time_lc;
                        float time_seconds = std::chrono::duration_cast<std::chrono::duration<float>>(time_ms).count();
                        printf("time cost: %.4f seconds\n", time_seconds);
                        auto js = jsrecords(lc, y_star, y_sub, y_rivals, puzzle_id);
                        js["time"] = time_seconds;
                        js["seed"] = SEED_RAND;
                        js["is_duplicated"] = found_ud;
                        ud_ipairs.insert(pairs2string(lc.ps_inside));
                        std::vector<std::vector<std::pair<int, int>>> uk_pairs;
                        for(auto ipairs: ipairs_subsets){
                            if(ds_ipairs.find(pairs2string(ipairs))==ds_ipairs.end()){
                                uk_pairs.push_back(ipairs);
                            }
                        }
                        if(uk_pairs.size()){
                            js["ismin"] = false;
                            js["uk_ipairs"] = uk_pairs;
                        }else{
                            js["ismin"] = true;
                        }
                        std::string jstring = js.dump();
                        outputFile << jstring << std::endl;
                        records.push_back(jstring);
                        if (!found_ud){
                            undesignableLibFile << jstring << std::endl;
                        }
                    }
                    for(auto pair: ds_ipairs)
                        std::cout<<pair<<"  ";
                    delete tree;
                    printf("\n");
                }
                auto end_time = std::chrono::high_resolution_clock::now();
                const std::chrono::duration<double, std::milli> time_ms_y = end_time - start_time;
                float time_seconds_lc = std::chrono::duration_cast<std::chrono::duration<float>>(time_ms_y).count();
                printf("time cost for whole structure: %.4f seconds\n", time_seconds_lc);
                timeFile << puzzle_id << "," << time_seconds_lc <<std::endl;
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
                        std::string r = puzzle_id+","+y_star+","+y_prim+","+std::to_string(count_pairs)+","+std::to_string(lc.node->first)+","+std::to_string(lc.node->second);
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
                        std::string r = puzzle_id+","+y_star+","+y_prim+","+std::to_string(count_pairs)+","+std::to_string(lc.node->first)+","+std::to_string(lc.node->second);
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
                std::string file_m1 = replaceFileExtension(csv, "m1");
                std::map<std::string, std::set<std::string>> id2m1 = readMotif(file_m1.c_str());
                auto start_time = std::chrono::high_resolution_clock::now();
                std::vector<LoopComplex> lc_list;
                TreeNode* root = parseStringToTree(y_star);
                // int max_internal = max_single(root);
                // if(max_internal > 30){
                //     std::string r = puzzle_id+","+std::to_string(max_internal);
                //     records.push_back(r);
                //     outputFile << r << std::endl;
                //     continue;
                // }
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
                    std::string result;
                    try{
                        result = alg_5_helper_v2(target, lc.ref, lc.constr, subseq, verbose, dangle);
                    }catch(...){
                        // Code to handle any exception
                        std::cerr << "An exception occurred" << std::endl;
                        result = "exception";
                    }
                    if (result == "exception")
                        continue;
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
                        std::string jstr = lc.jsmotif(puzzle_id);
                        std::cout << jstr << std::endl;
                        outputFile << jstr << std::endl;
                        // break;
                        // return;
                    }
                    printf("\n");
                }
            }
            if (alg == "nb1plot" || alg == "nb2plot" || alg == "nb3plot"){
                std::string file_m1 = replaceFileExtension(csv, "m1");
                std::map<std::string, std::set<std::string>> id2m1 = readMotif(file_m1.c_str());
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
                if (alg == "nb1plot")
                    tree2Edges(root, y_star, lc_list);
                else if (alg == "nb2plot")
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
    std::cout<< "count designable: " << count_designable << std::endl;
    std::cout << "designable records: " << records_designable.size() << std::endl;
    // Close the file
    outputFile.close();
    designableLibFile.close();
    std::cout << "Strings written to file: " << fileName << std::endl;
}

void online_process(std::string y, std::string path_prefix=""){
    // auto df = read_csv(csv.c_str());
    // printf("df shape: %d, %d\n", df.size(), df[0].size());
    std::string time_stamp = getCurrentTimestamp();
    if(path_prefix.empty())
        path_prefix = time_stamp;
    std::string alg = "pn";
    std::vector<std::string> y_list;
    y_list.push_back(y);
    std::cout<<"y_list size: "<<y_list.size()<<std::endl;
    std::vector<std::string> records;
    std::vector<std::string> records_designable;
    int count_designable = 0;
    // Specify the file name
    std::string fileName = path_prefix + "." + alg + ".log."+time_stamp+".txt";
    #ifdef SPECIAL_HP
    #else
        fileName = path_prefix + "." + alg + ".log.nosh."+time_stamp+".txt";
    #endif

    // Output file stream
    std::ofstream outputFile(fileName);
    // Check if the file is open
    if (!outputFile.is_open()) {
        std::cerr << "Error opening the file: " << fileName << std::endl;
        return;
    }
    // Library files for designable and undesignable motifs
    std::ofstream designableLibFile("lib_designable.txt", std::ios::app);
    std::ofstream undesignableLibFile("lib_undesignable.txt", std::ios::app);
    if (!designableLibFile.is_open()) {
        std::cerr << "Error opening the file: " << "lib_designable.txt" << std::endl;
        return;
    }
    if (!undesignableLibFile.is_open()) {
        std::cerr << "Error opening the file: " << "lib_undesignable.txt" << std::endl;
        return;
    }

    std::string fileTime = path_prefix + "." + alg + ".time."+getCurrentTimestamp()+".csv";
    #ifdef SPECIAL_HP
    #else
        fileTime = path_prefix + "." + alg + ".time.nosh."+getCurrentTimestamp()+".csv";
    #endif
    // Open the file for writing
    std::ofstream timeFile(fileTime);

    // Check if the file is open
    if (!timeFile.is_open()) {
        std::cerr << "Error opening the file: " << fileTime << std::endl;
        return;
    }else{
        timeFile << "ID,Time(s)" << std::endl;
    }
    std::unordered_map<std::string, GroupY> constr2groupy;
    // std::unordered_map<std::string, std::string> uniq_ud;
    const char*  var_undesignable_lib = std::getenv("PATH_UNDESIGNABLE_LIB");
    const char*  var_designable_lib = std::getenv("PATH_DESIGNABLE_LIB");
    if(var_undesignable_lib == NULL){
        std::cerr << "Error: PATH_UNDESIGNABLE_LIB is not set" << std::endl;
        return;
    }
    if(var_designable_lib == NULL){
        std::cerr << "Error: PATH_DESIGNABLE_LIB is not set" << std::endl;
        return;
    }
    std::string path_undesignable_lib(var_undesignable_lib);
    std::string path_designable_lib(var_designable_lib);
    std::set<std::string> uniq_ud = loadLib(path_undesignable_lib); // undesignable motifs
    std::set<std::string> uniq_ds = loadLib(path_designable_lib);   // designable   motifs
    for(int i = 0; i < y_list.size(); i++){
        // auto row = df[i];
        {
            string sharph0 = "()";
            string sharph1 = "(.)";
            string sharph2 = "(..)";

            std::string puzzle_id = path_prefix;
            std::string seq;
            std::string y_star = y_list[i];
            if(y_star.find(sharph0) != std::string::npos || y_star.find(sharph1) != std::string::npos || y_star.find(sharph2) != std::string::npos){
                std::cout<<"sharp turn!"<<std::endl;
                std::cout<<y_star<<std::endl;
                continue;
            }
            // if(row.size() > 2)
            //     seq = row[2];
            // else
            seq = tg_init(y_star);
            std::string y_prim; // currently dummy

            std::cout<<"Puzzle ID: "<<puzzle_id<<std::endl;

            // std::string goal_test = "p [] (M [0, 5, 0] (p [] (), B [1, 0] (S [0, 0] (p [] ()))))";
            // power neighbor set search
            if (alg == "pn"){
                auto start_time = std::chrono::high_resolution_clock::now();
                TreeNode* root = parseStringToTree(y_star);
                std::set<string> ds_ipairs; // designable internal pairs
                std::set<string> ud_ipairs; // undesinable internal pairs

                // motif of 2 loops (1 edge/internal pair)
                std::vector<LoopComplex> lc_list;
                tree2Edges(root, y_star, lc_list);
                printf("lc_list size: %d\n", lc_list.size());
                // sort motifs by the number of unknown positions
                std::sort(lc_list.begin(), lc_list.end(), [](const LoopComplex &a, const LoopComplex &b) {
                    return a.count_uk < b.count_uk;});

                // motif of 3 loops (2 edges/internal pairs)
                std::vector<LoopComplex> lc_list_2nbs;
                tree2TwoNeighbor(root, y_star, lc_list_2nbs);
                printf("lc_list_2nbs size: %d\n", lc_list_2nbs.size());
                // sort motifs by the number of unknown positions
                std::sort(lc_list_2nbs.begin(), lc_list_2nbs.end(), [](const LoopComplex &a, const LoopComplex &b) {
                    return a.count_uk < b.count_uk;});
                lc_list.insert(lc_list.end(), lc_list_2nbs.begin(), lc_list_2nbs.end());

                // motif of 4 loops (3 edges/internal pairs)
                std::vector<LoopComplex> lc_list_3nbs;
                tree2ThreeNeighbor(root, y_star, lc_list_3nbs);
                printf("lc_list_3nbs size: %d\n", lc_list_3nbs.size());
                // sort motifs by the number of unknown positions
                std::sort(lc_list_3nbs.begin(), lc_list_3nbs.end(), [](const LoopComplex &a, const LoopComplex &b) {
                    return a.count_uk < b.count_uk;});
                lc_list.insert(lc_list.end(), lc_list_3nbs.begin(), lc_list_3nbs.end());

                for (auto lc: lc_list){
                    auto start_time_lc = std::chrono::high_resolution_clock::now();
                    lc.printLoopLens();
                    if (lc.hasLongLoop()){
                        printf("the loop exceeds length limit!");
                        continue;
                    }
                    std::string target = y_star.substr(lc.start, lc.end-lc.start+1);
                    std::string subseq = seq.substr(lc.start, lc.end-lc.start+1);
                    printf(" count: %d\n", lc.count_uk);
                    printf("target: %s\n", target.c_str());
                    printf("   ref: %s\n", lc.ref.c_str());
                    printf("constr: %s\n", lc.constr.c_str());
                    auto ipairs_subsets = pairSubSet(lc.ps_inside);
                    std::string result;
                    // check if the (rotated) motif is already found undesignable
                    json js_motif = json::parse(lc.jsmotif(puzzle_id));
                    std::cout<<"js_motif: "<<js_motif<<std::endl;
                    Node* tree = new Node(js_motif);
                    std::string treestr = tree->toDotBracket();
                    // if(!goal_test.empty() && treestr != goal_test){
                    //     delete tree;
                    //     continue;
                    // }
                    std::cout<<"treestr: "<<treestr<<std::endl;
                    if(uniq_ds.find(treestr) != uniq_ds.end()){
                        std::cout<<"already designable!"<<std::endl;
                        result = "designable";
                    }else if(uniq_ud.find(treestr) != uniq_ud.end()){
                        result = "undesignable";
                        std::cout<<"recur lc.constr: "<<lc.constr<<std::endl;
                        // std::cout<<"recur    groupy: "<<constr2groupy[lc.constr].constr<<std::endl;
                        std::cout<<"recur lc.constr: "<<target<<std::endl;
                        std::cout<<"recur   treestr: "<<treestr<<std::endl;
                        // std::cout<<"recur     ystar: "<<constr2groupy[lc.constr].star<<std::endl;
                    }else{
                        std::string ref_lc = lc.ref;
                        std::string constr_lc = lc.constr;
                        for( int ib = 1; ib < lc.ps_outside.size(); ib++){ // skip the most outside boudary pair
                            auto bpair = lc.ps_outside[ib];
                            std::cout<<"bpair: "<<bpair.first<<"\t"<<bpair.second<<std::endl;
                            int len_branch =  bpair.second - bpair.first + 1;
                            std::string y_branch = y_star.substr(bpair.first, len_branch);
                            std::cout<<"y_branch:"<<std::endl;
                            std::cout<<y_branch<<std::endl;
                            TreeNode* root_branch = parseStringToTree(y_branch);
                            if(max_single(root_branch) > SINGLE_MAX_LEN || max_multi(root_branch) > MULTIPLE_FIRST_MAX_LEN){
                                std::string helix_branch = genHelix(len_branch);
                                std::string seq_branch = tg_init(helix_branch);
                                target.replace(bpair.first - lc.start, len_branch, helix_branch);
                                ref_lc.replace(bpair.first - lc.start, len_branch, helix_branch);
                                constr_lc.replace(bpair.first - lc.start, len_branch, helix_branch);
                                subseq.replace(bpair.first - lc.start, len_branch, seq_branch);
                            }
                        }
                        printf("target: %s\n", target.c_str());
                        printf("   ref: %s\n", ref_lc.c_str());
                        printf("constr: %s\n", constr_lc.c_str());
                        bool ud = false; // check if the motif is already found undesignable
                        for(auto ipairs: ipairs_subsets){
                            if(ud_ipairs.find(pairs2string(ipairs))!=ud_ipairs.end()){
                                ud = true;
                                break;
                            }
                        }
                        if(ud)
                            continue;
                        try{
                            std::cout<<"UMFE log: before alg 5 "<<std::endl;
                            result = alg_5_helper_v2(target, ref_lc, constr_lc, subseq, verbose, dangle);
                            if (result == "UMFE"){
                                result = "designable";
                                std::cout<<"UMFE log: afer alg 5 "<<std::endl;
                                std::cout<<"UMFE log: "<<tree->toDotBracket()<<std::endl;
                                std::cout<<"UMFE log: "<<target<<std::endl;
                                std::cout<<"UMFE log: "<<constr_lc<<std::endl;
                            }
                        }catch(...){
                            // Code to handle any exception
                            std::cerr << "An exception occurred" << std::endl;
                            result = "exception";
                            continue;
                        }
                    }
                    if (result == "designable"){
                        std::cout<<"designable!"<<std::endl;
                        ds_ipairs.insert(pairs2string(lc.ps_inside));
                        count_designable++;
                        y_sub = y_star; // set y_sub as y_star
                        y_rivals.clear(); // clear y_rivals
                        bool found_ds = false;
                        if(uniq_ds.find(treestr) != uniq_ds.end()){
                            found_ds = true;
                            std::cout<<"already found designable!"<<std::endl;
                        }else{
                            uniq_ds.insert(treestr);
                            for(Node* rotree: tree->rotated(0)){
                                std::string rotreestr = rotree->toDotBracket();
                                uniq_ds.insert(rotreestr);
                                delete rotree;
                            }
                            auto end_time_lc = std::chrono::high_resolution_clock::now();
                            const std::chrono::duration<double, std::milli> time_ms = end_time_lc - start_time_lc;
                            float time_seconds = std::chrono::duration_cast<std::chrono::duration<float>>(time_ms).count();
                            printf("time cost: %.4f seconds\n", time_seconds);
                            auto js = jsrecords(lc, y_star, y_sub, y_rivals, puzzle_id);
                            js["time"] = time_seconds;
                            js["seed"] = SEED_RAND;
                            js["is_duplicated"] = found_ds;
                            std::string jstring = js.dump();
                            designableLibFile << jstring << std::endl;
                            records_designable.push_back(jstring);
                        }
                    }
                    if (result == "undesignable"){
                        std::cout<<"undesignable!"<<std::endl;
                        bool found_ud = false; // check if the motif is already found undesignable
                        // if(constr2groupy.find(lc.constr) != constr2groupy.end() && constr2groupy[lc.constr].star == target)
                        if(uniq_ud.find(treestr) != uniq_ud.end()){
                            found_ud = true;
                        }else{
                            uniq_ud.insert(treestr);
                            for(Node* rotree: tree->rotated(0)){
                                std::string rotreestr = rotree->toDotBracket();
                                uniq_ud.insert(rotreestr);
                                delete rotree;
                            }
                        }
                        auto end_time_lc = std::chrono::high_resolution_clock::now();
                        const std::chrono::duration<double, std::milli> time_ms = end_time_lc - start_time_lc;
                        float time_seconds = std::chrono::duration_cast<std::chrono::duration<float>>(time_ms).count();
                        printf("time cost: %.4f seconds\n", time_seconds);
                        auto js = jsrecords(lc, y_star, y_sub, y_rivals, puzzle_id);
                        js["time"] = time_seconds;
                        js["seed"] = SEED_RAND;
                        js["is_duplicated"] = found_ud;
                        ud_ipairs.insert(pairs2string(lc.ps_inside));
                        std::vector<std::vector<std::pair<int, int>>> uk_pairs;
                        for(auto ipairs: ipairs_subsets){
                            if(ds_ipairs.find(pairs2string(ipairs))==ds_ipairs.end()){
                                uk_pairs.push_back(ipairs);
                            }
                        }
                        if(uk_pairs.size()){
                            js["ismin"] = false;
                            js["uk_ipairs"] = uk_pairs;
                        }else{
                            js["ismin"] = true;
                        }
                        std::string jstring = js.dump();
                        outputFile << jstring << std::endl;
                        records.push_back(jstring);
                        if (!found_ud){
                            undesignableLibFile << jstring << std::endl;
                        }
                    }
                    for(auto pair: ds_ipairs)
                        std::cout<<pair<<"  ";
                    delete tree;
                    printf("\n");
                }
                auto end_time = std::chrono::high_resolution_clock::now();
                const std::chrono::duration<double, std::milli> time_ms_y = end_time - start_time;
                float time_seconds_lc = std::chrono::duration_cast<std::chrono::duration<float>>(time_ms_y).count();
                printf("time cost for whole structure: %.4f seconds\n", time_seconds_lc);
                timeFile << puzzle_id << "," << time_seconds_lc <<std::endl;
            }
        }
    }
    for (auto r: records)
        std::cout<<r<<std::endl;
    std::cout<< "count designable: " << count_designable << std::endl;
    std::cout << "designable records: " << records_designable.size() << std::endl;
    // Close the file
    outputFile.close();
    designableLibFile.close();
    std::cout << "Strings written to file: " << fileName << std::endl;
}

void show_configuration(){
    #ifdef SPECIAL_HP
        printf("SPECIAL_HP   defined.\n");
        #define SPECIAL_HP_3
        #define SPECIAL_HP_4
        #define SPECIAL_HP_6
    #else
        printf("SPECIAL_HP undefined.\n");
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
    #ifdef SINGLE_MAX_LEN
        printf("SINGLE_MAX_LEN: %d\n", SINGLE_MAX_LEN);
    #endif
    #ifdef MULTIPLE_FIRST_MAX_LEN
        printf("MULTIPLE_FIRST_MAX_LEN: %d\n", MULTIPLE_FIRST_MAX_LEN);
    #endif
    return;
}

int main(int argc, char* argv[]) {
    cxxopts::Options options("MyProgram", "One line description of MyProgram");
    options.add_options()
    ("a,alg", "Algorithm", cxxopts::value<std::string>()->default_value("0"))
    ("c,csv", "csv file", cxxopts::value<std::string>()->default_value(""))
    ("t,txt", "txt file", cxxopts::value<std::string>()->default_value(""))
    ("s,seed", "random seed", cxxopts::value<int>()->default_value("2024"))
    ("d,dangle", "Dangle mode", cxxopts::value<int>()->default_value("2"))
    ("vrna,vienna", "Use ViennaRNA to fold, the environment variable VRNABIN has to be set", cxxopts::value<bool>()->default_value("false"))
    ("v,verbose", "Verbose output", cxxopts::value<bool>()->default_value("false"));

    auto result = options.parse(argc, argv);
    std::string alg = result["alg"].as<std::string>();
    std::string csv = result["csv"].as<std::string>();
    std::string txt = result["txt"].as<std::string>();
    bool vrna = result["vienna"].as<bool>();
    verbose = result["verbose"].as<bool>();
    dangle = result["dangle"].as<int>();
    SEED_RAND = result["seed"].as<int>();
    printf("alg: %s, vienna: %d, verbose: %d, dangle: %d\n", alg.c_str(), vrna, verbose, dangle);
    printf("random seed for target initialization: %d\n", SEED_RAND);
    show_configuration();

    std::unordered_map<std::string, std::string> struct2seq = loadlib_eterna("data/eterna_umfe_unsolved.csv");

    if (alg == "0"){
        std::cout<<"no alg was selected!"<<std::endl;
        return 0;
    }

    if (!csv.empty()){
        csv_process(csv, alg);
        return 0;
    }


    if ( alg == "csfold" || alg == "cf" ){  /* constrained folding */
        std::cout << alg << std::endl;
        int beamsize = 0;
        bool sharpturn = false;
        std::string seq;
        std::string constr;

        while (std::getline(std::cin, seq))
        {
            std::getline(std::cin, constr);
            std::vector<std::string> refs;
            if (alg == "csfold")
                refs =  cs_fold(seq, constr, beamsize, sharpturn, verbose, dangle);
            else
                refs =  cs_fold_vienna(seq, constr, beamsize, sharpturn, verbose, dangle);
            std::cout<<"subopts size: "<<refs.size()<<std::endl;
            for(auto ref: refs)
                std::cout<<ref<<std::endl;
        }
        return 0;
    }else if (alg == "fold"){ /* fold */
        int beamsize = 0;
        bool sharpturn = false;
        std::string seq;
        std::vector<std::string> refs;
        while (std::getline(std::cin, seq))
        {   
            if (!vrna)
                refs =  fold(seq, beamsize, sharpturn, verbose, dangle, 0.);
            else{
                std::string constr(seq.size(), '?');
                std::cout<<constr<<std::endl;
                refs = subopt_vienna(seq, constr);
            }
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
        while (std::getline(std::cin, ref1)) {
            // Process the line as needed
            // std::cout<<"got seq"<<std::endl;
            // getline(std::cin, ref1);
            getline(std::cin, ref2);
            auto start_time = std::chrono::high_resolution_clock::now();
            // std::vector<std::string> X = alg_1_v2(ref1, ref2, seq, verbose, dangle);
            seq = tg_init(ref1);
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
    }
    else if (alg == "loop"){ /* loops evaluation  */
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
    }
    else if (alg == "foldv"){ /* loops evaluation  */
        std::string seq;
        // Read input line by line until EOF (end of file) is reached
        while (std::getline(std::cin, seq)) {
                std::string mfe = fold_vienna(seq);
                std::cout<<seq<<std::endl;
                std::cout<<mfe<<std::endl;
            }
    }
    else if (alg == "online"){ /* loops evaluation  */
        std::string y;
        // Read input line by line until EOF (end of file) is reached
        while (std::getline(std::cin, y)) {
                std::cout<<y<<std::endl;
                online_process(y, "online");
            }
    }
    else if (alg == "mloop"){ /* multi-loops evaluation  */
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
        while (std::getline(std::cin, target)) {
            // getline(std::cin, target);
            getline(std::cin, cst);
            // getline(std::cin, ref);
            seq = tg_init(target);
            ref = cst;
            std::replace(ref.begin(), ref.end(), '?', '.');
            printf(" count: %d\n", countOccurrences(cst, '?'));
            printf("target: %s\n", target.c_str());
            printf("constr: %s\n", cst.c_str());
            printf("   ref: %s\n", ref.c_str());
            TreeNode* root = parseStringToTree(target);
            int max_internal = max_single(root);
            // if(max_internal > 30){        
            //     std::cout<<"the internal loop is too long: "<<max_internal<<std::endl;            
            //     continue;
            // }
            std::cout<<"max internal loop length: "<<max_internal<<std::endl;
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
    }else if (alg == "n1" || alg == "n2" || alg == "n3"){ /* edges evaluation  */
        std::string ref;
        // Read input line by line until EOF (end of file) is reached
        while (std::getline(std::cin, ref)){
            auto start_time = std::chrono::high_resolution_clock::now();
            std::string seq;
            std::string puzzle_id = "id";
            if(struct2seq.find(ref) != struct2seq.end()){
                seq = struct2seq[ref];
                std::cout<<"seq found in design lib: "<<seq<<std::endl;
            }
            else{
                seq = tg_init(ref);
                std::cout<<"seq via targeted initialization: "<<seq<<std::endl;
            }
            std::vector<LoopComplex> lc_list;
            TreeNode* root = parseStringToTree(ref);
            if(alg == "n1")
                tree2Edges(root, ref, lc_list);
            else if(alg == "n2")
                tree2TwoNeighbor(root, ref, lc_list);
            else{
                assert (alg == "n3");
                tree2ThreeNeighbor(root, ref, lc_list);
            }
            printf("lc_list size: %d\n", lc_list.size());
            // Sort the vector using a lambda expression
            std::sort(lc_list.begin(), lc_list.end(), [](const LoopComplex &a, const LoopComplex &b) {
                return a.count_uk < b.count_uk;});
            std::vector<std::string> motif_records;
            for (auto lc: lc_list){
                 auto start_time_lc = std::chrono::high_resolution_clock::now();
                lc.printLoopLens();
                if (lc.hasLongLoop()){
                    printf("the loop exceeds length limit!");
                    continue;
                }
                std::string target = ref.substr(lc.start, lc.end-lc.start+1);

                std::string subseq = seq.substr(lc.start, lc.end-lc.start+1);
                printf(" count: %d\n", lc.count_uk);
                printf("target: %s\n", target.c_str());
                printf("   ref: %s\n", lc.ref.c_str());
                printf("constr: %s\n", lc.constr.c_str());
                for(auto bp: lc.ps_outside)
                    printf("bpair: %d, %d\n", bp.first, bp.second);
                for(auto ip: lc.ps_inside)
                    printf("ipair: %d, %d\n", ip.first, ip.second);
                std::string ref_lc = lc.ref;
                std::string constr_lc = lc.constr;
                for( int ib = 1; ib < lc.ps_outside.size(); ib++){ // skip the most outside boudary pair
                    auto bpair = lc.ps_outside[ib];
                    // std::cout<<"bpair: "<<bpair.first<<"\t"<<bpair.second<<std::endl;
                    int len_branch =  bpair.second - bpair.first + 1;
                    std::string y_branch = ref.substr(bpair.first, len_branch);
                    // std::cout<<"y_branch:"<<std::endl;
                    std::cout<<y_branch<<std::endl;
                    TreeNode* root_branch = parseStringToTree(y_branch);
                    if(max_single(root_branch) > SINGLE_MAX_LEN || max_multi(root_branch) > MULTIPLE_FIRST_MAX_LEN){
                        std::string helix_branch = genHelix(len_branch);
                        std::string seq_branch = tg_init(helix_branch);
                        target.replace(bpair.first - lc.start, len_branch, helix_branch);
                        ref_lc.replace(bpair.first - lc.start, len_branch, helix_branch);
                        constr_lc.replace(bpair.first - lc.start, len_branch, helix_branch);
                        subseq.replace(bpair.first - lc.start, len_branch, seq_branch);
                    }
                }

                std::string result = alg_5_helper_v2(target, ref_lc, constr_lc, subseq, verbose, dangle);
                if (result == "undesignable"){
                    std::cout<<"undesignable!"<<std::endl;
                    auto end_time_lc = std::chrono::high_resolution_clock::now();
                    const std::chrono::duration<double, std::milli> time_ms_lc = end_time_lc - start_time_lc;
                    float time_seconds_lc = std::chrono::duration_cast<std::chrono::duration<float>>(time_ms_lc).count();
                    printf("time cost for the motif: %.4f seconds\n", time_seconds_lc);
                    auto js = jsrecords(lc, ref, y_sub, y_rivals, puzzle_id);
                    js["time"] = time_seconds_lc;
                    js["ismin"] = true;
                    std::string jstring = js.dump();
                    motif_records.push_back(jstring);
                }
            }
            printf("----------------------------------------------------------------------------------\n");
            if(motif_records.size()){
                std::cout<<"found "+std::to_string(motif_records.size())+ " undesignable motif(s)."<<std::endl;
                for(int i=0; i < motif_records.size(); i++){
                    std::cout<<"motif: "<<std::to_string(i+1)<<std::endl;
                    std::cout<<motif_records[i]<<std::endl;
                }
                auto end_time = std::chrono::high_resolution_clock::now();
                const std::chrono::duration<double, std::milli> time_ms_y = end_time - start_time;
                float time_seconds_y = std::chrono::duration_cast<std::chrono::duration<float>>(time_ms_y).count();
                for(int i=0; i < motif_records.size(); i++){
                    float jtime = json::parse(motif_records[i])["time"];
                    printf("time cost for the motif: %.4f seconds\n", jtime);
                }
                printf("time cost for the structure: %.4f seconds\n", time_seconds_y);
            }else{
                std::cout<<"no undesignable motifs found."<<std::endl;
            }
        }
    }else if (alg == "showtree"){
        std::string ref;
        // Read input line by line until EOF (end of file) is reached
        while (std::getline(std::cin, ref)) {
            TreeNode* root = parseStringToTree(ref);
            root->printTree();
        }
    }else if (alg == "subset1" || alg == "subset2" || alg == "subset3"){
        std::string ref;
        // Read input line by line until EOF (end of file) is reached
        while (std::getline(std::cin, ref)){
            std::string seq;
            if(struct2seq.find(ref) != struct2seq.end()){
                seq = struct2seq[ref];
                std::cout<<"seq found in design lib: "<<seq<<std::endl;
            }
            else{
                seq = tg_init(ref);
                std::cout<<"seq via targeted initialization: "<<seq<<std::endl;
            }
            std::vector<LoopComplex> lc_list;
            TreeNode* root = parseStringToTree(ref);
            if(alg == "subset1")
                tree2Edges(root, ref, lc_list);
            else if(alg == "subset2")
                tree2TwoNeighbor(root, ref, lc_list);
            else
                tree2ThreeNeighbor(root, ref, lc_list);
            printf("lc_list size: %d\n", lc_list.size());
            // Sort the vector using a lambda expression
            std::sort(lc_list.begin(), lc_list.end(), [](const LoopComplex &a, const LoopComplex &b) {
                return a.count_uk < b.count_uk;});
            std::vector<std::vector<std::string>> motif_records;
            for (auto lc: lc_list){
                std::cout<<lc.constr<<std::endl;
                std::cout<<lc.ps_inside.size()<<std::endl;
                for(auto pair: lc.ps_inside){
                    std::cout<<pair.first<<","<<pair.second<<";";
                }
                std::cout<<std::endl;
                auto subsets = pairSubSet(lc.ps_inside);
                for(auto subset: subsets){
                    for(auto pair: subset){
                        std::cout<<pair.first<<","<<pair.second<<";";
                    }
                    std::cout<<std::endl;
                }
            }
        }
    }else if (alg == "imax"){
        std::string y;
        // Read input line by line until EOF (end of file) is reached
        while (std::getline(std::cin, y)){
            TreeNode* root = parseStringToTree(y);
            int max_internal = max_single(root);
            printf("maximum internal loop length: %d\n", max_internal);
        }
    }else if (alg == "mmax"){
        std::string y;
        // Read input line by line until EOF (end of file) is reached
        while (std::getline(std::cin, y)){
            TreeNode* root = parseStringToTree(y);
            int max_mul = max_multi(root);
            printf("maximum multi-loop length: %d\n", max_mul);
        }
    }else if (alg == "l1" || alg == "l2" || alg == "l3"){ /* edges evaluation  */
        std::string ref;
        // Read input line by line until EOF (end of file) is reached
        while (std::getline(std::cin, ref)){
            std::string seq;
            std::string puzzle_id = "id";
            if(struct2seq.find(ref) != struct2seq.end()){
                seq = struct2seq[ref];
                std::cout<<"seq found in design lib: "<<seq<<std::endl;
            }
            else{
                seq = tg_init(ref);
                std::cout<<"seq via targeted initialization: "<<seq<<std::endl;
            }
            std::vector<LoopComplex> lc_list;
            TreeNode* root = parseStringToTree(ref);
            if(alg == "l1")
                tree2Edges(root, ref, lc_list);
            else if(alg == "l2")
                tree2TwoNeighbor(root, ref, lc_list);
            else
                tree2ThreeNeighbor(root, ref, lc_list);
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
                lc.printLoopLens();
                if (lc.hasLongLoop()){
                    printf("the loop exceeds length limit!\n");
                    continue;
                }
                std::string ref_lc = lc.ref;
                std::string constr_lc = lc.constr;
                for( int ib = 1; ib < lc.ps_outside.size(); ib++){ // skip the most outside boudary pair
                    auto bpair = lc.ps_outside[ib];
                    // std::cout<<"bpair: "<<bpair.first<<"\t"<<bpair.second<<std::endl;
                    int len_branch =  bpair.second - bpair.first + 1;
                    std::string y_branch = ref.substr(bpair.first, len_branch);
                    // std::cout<<"y_branch:"<<std::endl;
                    std::cout<<y_branch<<std::endl;
                    TreeNode* root_branch = parseStringToTree(y_branch);
                    if(max_single(root_branch) > SINGLE_MAX_LEN || max_multi(root_branch) > MULTIPLE_FIRST_MAX_LEN){
                        std::string helix_branch = genHelix(len_branch);
                        std::string seq_branch = tg_init(helix_branch);
                        target.replace(bpair.first - lc.start, len_branch, helix_branch);
                        ref_lc.replace(bpair.first - lc.start, len_branch, helix_branch);
                        constr_lc.replace(bpair.first - lc.start, len_branch, helix_branch);
                        subseq.replace(bpair.first - lc.start, len_branch, seq_branch);
                    }
                }
                std::cout<<subseq<<std::endl;
                std::cout<<constr_lc<<std::endl;
                std::vector<std::string>  refs =  cs_fold(subseq, constr_lc, 0, false, verbose, dangle);
                std::cout<<"subopts size: "<<refs.size()<<std::endl;
                for(auto ref: refs)
                    std::cout<<ref<<std::endl;
            }
        }
    }else if (alg == "helix"){
        std::string lenstr;
        while(std::getline(std::cin, lenstr)){
            int len = std::stoi(lenstr);
            std::string h = genHelix(len);
            std::string x = tg_init(h);
            std::cout<<h<<std::endl;
            std::cout<<x<<std::endl;
        }
    }
    return 0;
}