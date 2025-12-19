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
#include <cstdio>
// #include <random>

// #include "eval.cpp"
#include "cxxopts.hpp"
#include "utils.h"
#include "comps.h"
#include "eval.h"
#include <thread>
using namespace std;

double GASCONST = 1.98717;
double TEMPERATURE = 310.15;

long int MAX_ENUM = 10000000000; // default 10000000000
long int MAX_CONSTRAINT = 100000; // default 100000
long int N_SAMPLE = 500;  // default 500
long int MAX_RIVAL = 100; // default 100

const char*  var_undesignable_lib = std::getenv("PATH_UNDESIGNABLE_LIB");
const char*  var_designable_lib = std::getenv("PATH_DESIGNABLE_LIB");
const char*  var_unknown_lib = std::getenv("PATH_UNKNOWN_LIB");

std::ofstream designableLibFile; // ("lib_designable.txt", std::ios::app);
std::ofstream undesignableLibFile; // ("lib_undesignable.txt", std::ios::app);

std::set<std::string> uniq_ud; // = loadLib(path_undesignable_lib); // undesignable motifs
std::set<std::string> uniq_ds; // = loadLib(path_designable_lib);   // designable   motifs
std::set<std::string> uniq_unknown; // = loadLib(var_unknown_lib); // unknown motifs


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
bool PLOT = false;
int dangle = 2;

int y_sub_start;
int y_sub_end;
std::string y_sub;
std::string seq_init;
std::vector<std::string> y_rivals;
std::string x_umfe;

// cache mfe solution for alg2
std::string mfe_seq_for_alg2;

std::vector<std::pair<int, int>> pairs_outside;
std::vector<std::pair<int, int>> pairs_inside;

int SEED_RAND = 2024;

class MotifLib{
    public:
        std::map<std::string, int> dotbracket2id;
        int count_unique = 0;
        
        MotifLib(std::string path);

        int getID(std::string dotbracket){
            return dotbracket2id[dotbracket];
        }

        int getCount(){
            return count_unique;
        }

        // int insert(std::string dotbracket){
        //     if (dotbracket2id.count(dotbracket)==0){
        //         dotbracket2id[dotbracket] = count_unique+1;
        //         count_unique += 1;
        //         return 1;
        //     }
        //     return 0;
        // }

        int insert(Node* tree){
            std::string treestr = tree->toDotBracket();
            if (dotbracket2id.count(treestr)==0){
                std::set<std::string> rostrings;
                dotbracket2id[treestr] = count_unique+1;
                for(Node* rotree : tree->rotated(0)){
                    std::string rotreestr = rotree->toDotBracket();
                    if (dotbracket2id.count(rotreestr)==0)
                        dotbracket2id[rotreestr] = count_unique + 1;
                    else
                        assert (rostrings.count(rotreestr)==1);
                    delete rotree;
                }
                count_unique += 1;
                return 1;
            }
            return 0;
        }

        bool count(std::string dotbracket){
            return dotbracket2id.count(dotbracket);
        }
};

// load lib from file
std::set<std::string> loadLib(std::string path){
    std::set<std::string> lib;
    std::ifstream file(path);
    std::string line;
    printf("Loading lib from file: %s\n", path.c_str());
    size_t line_number = 0;
    // start time
    auto start_time = std::chrono::high_resolution_clock::now();
    while (std::getline(file, line)){
        line_number += 1;
        // printf("path %s, line %d\n", path.c_str(), line_number);
        json js_record = json::parse(line);
        std::set<std::string> rostrings;
        // assert(!js_record["is_duplicated"]);
        Node* tree = new Node(js_record["motif"]);
        std::string treestr = tree->toDotBracket();
        rostrings.insert(treestr);
        // if (lib.count(treestr) != 0)
        //     std::cout<<"duplicated motif: "<<treestr<<std::endl;
        // assert (lib.count(treestr)==0);
        lib.insert(treestr);
        for(Node* rotree : tree->rotated(0)){
            std::string rotreestr = rotree->toDotBracket();
            if (lib.count(rotreestr)==0)
                lib.insert(rotreestr);
            // else
            //     assert (rostrings.count(rotreestr)==1);
            rostrings.insert(rotreestr);
            delete rotree;
        }
        delete tree;
    }
    // end time
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    printf("Loaded %d motifs from %s in %.2f seconds.\n", lib.size(), path.c_str(), elapsed.count());
    return lib;
}


void set_globals(){
    std::string path_undesignable_lib(var_undesignable_lib);
    std::string path_designable_lib(var_designable_lib);
    designableLibFile.open(path_designable_lib, std::ios::app);
    undesignableLibFile.open(path_undesignable_lib, std::ios::app);
    uniq_ud = loadLib(path_undesignable_lib); // undesignable motifs
    uniq_ds = loadLib(path_designable_lib);   // designable   motifs
    if (var_unknown_lib != NULL)
        uniq_unknown = loadLib(var_unknown_lib); // load unknown motifs if the library is set
    return;
}


// enrich lib with new motifs
std::set<std::string> enrichLib(std::string path, std::set<std::string>& lib){
    std::ifstream file(path);
    std::string line;
    while (std::getline(file, line)){
        json js_record = json::parse(line);
        std::vector<std::string> dotbrackets = js_record["dot-bracket"];
        for(std::string dotbracket : dotbrackets){
            if (lib.count(dotbracket)==0){
                lib.insert(dotbracket);
            }
        }
    }
    return lib;
}


static void loadDesignableEnum(const std::string &path, std::set<std::string> &uniq_ds){
    std::ifstream in(path);
    if (!in.is_open()) {
        std::cerr << "Error opening file: " << path << std::endl;
        return;
    }
    std::string line;
    while (std::getline(in, line)) {
        // strip whitespace
        // trim leading/trailing whitespace
        auto first = line.find_first_not_of(" \t\r\n");
        if (first == std::string::npos) 
            continue;              // blank line
        auto last = line.find_last_not_of(" \t\r\n");
        std::string raw = line.substr(first, last - first + 1);

        // insert raw JSON‐or‐dotbracket‐string
        uniq_ds.insert(raw);
        // if (!uniq_ds.insert(raw).second) { // already exists
        //     designableLibFile << raw << "\n";
        // }
    }
}


MotifLib::MotifLib(std::string path){
    std::ifstream file(path);
    std::string line;
    while (std::getline(file, line)){
        json js_record = json::parse(line);
        std::set<std::string> rostrings;
        assert(!js_record["is_duplicated"]);
        Node* tree = new Node(js_record["motif"]);
        std::string treestr = tree->toDotBracket();
        rostrings.insert(treestr);
        // assert (dotbracket2id.count(treestr)==0);
        if (js_record["id_uniq"].empty())
            js_record["id_uniq"] = count_unique + 1;
        dotbracket2id[treestr] = js_record["id_uniq"].get<int>();
        for(Node* rotree : tree->rotated(0)){
            std::string rotreestr = rotree->toDotBracket();
            if (dotbracket2id.count(rotreestr)==0)
                dotbracket2id[rotreestr] = js_record["id_uniq"];
            else
                assert (rostrings.count(rotreestr)==1);
            rostrings.insert(rotreestr);
            delete rotree;
        }
        count_unique += 1;
        delete tree;
    }
}

ulong count_enum(std::vector<std::tuple<int, int>>& pairs_diff);

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


double alg2_ensemble(std::string& y, std::vector<std::string>& y_rival_vector, std::vector<std::vector<std::vector<int>>>& cr_loops_vector, std::vector<std::tuple<int, int>>& pairs_diff, std::string& seq, bool is_verbose, int dangle_model){
    std::cout<<"inside alg2_ensemble"<<std::endl;
    assert (y_rival_vector.size() == cr_loops_vector.size());
    auto start = std::chrono::high_resolution_clock::now();
    ulong nEnum = count_enum(pairs_diff);
    std::cout<<"count enum: "<<nEnum<<std::endl;
    // get the upper bound \frac{1}{1 + \sum_{y' \in Y_{rival}} e^{-(\Delta G(\boldsymbol{x}, \boldsymbol{y'}) - \Delta G(\boldsymbol{x}, \boldsymbol{y^\star})) / R T}}
    // double prob_bound = 0.;
    int num_threads = std::thread::hardware_concurrency();
    double prob_bound_thread[num_threads];
    for (int t = 0; t < num_threads; t++)
        prob_bound_thread[t] = 0.;
    // long e_diff_min = 10000;
    #pragma omp parallel for
    for(ulong i=0; i < nEnum; i++){
        int thread_id = omp_get_thread_num();
        std::string seq_i = enumerate(pairs_diff, i, seq);
        if ((i+1)%1000000 == 0){
            auto pause = std::chrono::high_resolution_clock::now();
            printf("thread %2d, %8d, %.4f, %.2f seconds\n", thread_id, (i+1)/10000, prob_bound_thread[thread_id], std::chrono::duration<double, std::milli>(pause - start)/1000.f);
        }
        double denumerator = 1.;
        // long e_diff_max_local = 0;
        for (int j = 0; j < y_rival_vector.size(); j++){
            std::string y_prime = y_rival_vector[j];
            auto cr_loops = cr_loops_vector[j];
            if(check_compatible(seq_i, y_prime)){
                long e_diff = -diff_eval(seq_i, cr_loops, is_verbose, dangle_model); // not divided by 100, \delta G(x, y_star) - \delta G(x, y_prime)
                denumerator += std::exp(e_diff * 10 / (GASCONST * TEMPERATURE));
                // e_diff_max_local = std::max(e_diff_max_local, e_diff);
            }
        }
        // #pragma omp critical
        // {
        //     prob_bound = std::max(prob_bound, 1. / denumerator);
        //     e_diff_min = std::min(e_diff_min, e_diff_max_local);
        // }
        prob_bound_thread[thread_id] = std::max(prob_bound_thread[thread_id], 1. / denumerator);
    }
    // std::cout<<"e_diff_min: "<<e_diff_min<<std::endl;
    auto end = std::chrono::high_resolution_clock::now();
    // std::cout<<"time cost: "<<std::chrono::duration<double, std::milli>(end - start)/1000.f<<" seconds"<<std::endl;
    printf("time cost: %.4f seconds\n", std::chrono::duration<double, std::milli>(end - start)/1000.f);
    // double p_bound_local = 1 / (1. + std::exp(e_diff_min * 10 / (GASCONST * TEMPERATURE)));
    // std::cout<<"probability bound (local): "<<p_bound_local<<std::endl;
    double prob_bound = 0.;
    for (int t = 0; t < num_threads; t++)
        prob_bound = std::max(prob_bound, prob_bound_thread[t]);
    std::cout<<"probability bound: "<<prob_bound<<std::endl;
    return prob_bound;
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
    // long e_diff_min = 100;
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
            long e_diff = -diff_eval(seq_i, cr_loops, is_verbose, dangle_model); // not divided by 100, \delta G(y_star) - \delta G(y')
            if(e_diff < 0){
                #pragma omp critical
                idX.push_back({i, seq_i});
                // X.push_back(seq_i);
            }
            // #pragma omp critical
            // e_diff_min = std::min(e_diff_min, e_diff);
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
    // printf("alg1.e_diff_min: %ld\n", e_diff_min);
    return X;
}

std::string alg_2(std::string& ref1, std::set<std::string>& refs_checked, std::vector<Constraint>& cs_vec, bool verbose, int dangle_model){ // ref1, ref2, X, is_verbose, dangle_model
    int count_cs = cs_vec.size();
    std::vector<std::pair<ulong, std::pair<std::string, std::string>>> y_primes;
    std::cout<<"inside alg2"<<std::endl;
    std::vector<std::string> X;
    for(auto cs: cs_vec){
        X.insert(X.end(), cs.seqs->begin(), cs.seqs->end());
        if(X.size() > N_SAMPLE)
            break;
    }
    if (X.size() > N_SAMPLE)
        X.resize(N_SAMPLE);
    // if (X.size() > N_SAMPLE) {
    //     // static std::random_device rd;
    //     static std::mt19937 g(2025);
    //     std::shuffle(X.begin(), X.end(), g);
    //     X.resize(N_SAMPLE);
    // }
    std::cout<<"X.size: "<<X.size()<<std::endl;
    for(auto x: X){
        assert (check_compatible(x, ref1));
        std::vector<std::string> subopts = fold(x, 0, false, verbose, dangle_model, 0.0);
        if (isUMFE(subopts, ref1)){
            printf("UMFE found!\n");
            std::cout<<x<<std::endl;
            return "UMFE:" + x;
        }else if(isMFE(subopts, ref1)){
            printf("MFE found!\n");
            mfe_seq_for_alg2 = x;
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
            // update y_rivals
            y_rivals.push_back(y_prime.second.first);
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
    if (cs_vec.size() < MAX_RIVAL)
        return alg_2(ref1, refs_checked, cs_vec, verbose, dangle_model);
    else{
        std::cout<<"no conclusion!"<<std::endl;
        return "exceeded MAX_RIVAL: " + std::to_string(MAX_RIVAL);
    }
}

std::string alg_2_cs(std::string& ref1, std::set<std::string>& refs_checked, std::vector<Constraint>& cs_vec, bool verbose, int dangle_model){ // ref1, ref2, X, is_verbose, dangle_model
    int count_cs = cs_vec.size();
    std::vector<std::pair<ulong, std::pair<std::string, std::string>>> y_primes;
    std::cout<<"inside alg2cs"<<std::endl;
    std::vector<std::string> X;
    for(auto cs: cs_vec){
        X.insert(X.end(), cs.seqs->begin(), cs.seqs->end());
        if(X.size() > N_SAMPLE)
            break;
    }
    if (X.size() > N_SAMPLE)
        X.resize(N_SAMPLE);
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
    if (cs_vec.size() < MAX_RIVAL)
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
    y_rivals.clear();
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
            print(cs_ref2);
            cs_ref2.setStructure(ref2);
            cs_vec.push_back(cs_ref2);
            y_rivals.push_back(ref2);
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
        if(X.size() > N_SAMPLE)
            break;
    }
    if (X.size() > N_SAMPLE)
        X.resize(N_SAMPLE);
    std::cout<<"X.size: "<<X.size()<<std::endl;
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
    if (cs_vec.size() < MAX_RIVAL)
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
        if(X.size() > N_SAMPLE)
            break;
    }
    if (X.size() > N_SAMPLE)
        X.resize(N_SAMPLE);
    std::cout<<"X.size: "<<X.size()<<std::endl;
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
            x_umfe = x;
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
    if (cs_vec.size() < MAX_RIVAL)
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
        x_umfe = seq;
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

void csv_process(std::string csv, std::string alg){
    auto df = read_csv(csv.c_str());
    printf("df shape: %d, %d\n", df.size(), df[0].size());
    std::vector<std::string> records;
    std::vector<std::string> records_designable;
    int count_designable = 0;

    // Specify the file name
    std::string path_undesignable = csv + "." + alg + ".log."+getCurrentTimestamp()+".txt";
    std::string path_designable = csv + "." + alg + ".designable.log."+getCurrentTimestamp()+".txt";
    std::string path_unknown = csv + "." + alg + ".unknown.log."+getCurrentTimestamp()+".txt";
    
    #ifdef SPECIAL_HP
    #else
        path_undesignable = csv + "." + alg + ".log.nosh."+getCurrentTimestamp()+".txt";
        path_designable = csv + "." + alg + ".designable.log.nosh."+getCurrentTimestamp()+".txt";
        path_unknown = csv + "." + alg + ".unknown.log.nosh."+getCurrentTimestamp()+".txt";
    #endif

    // Output file stream
    std::ofstream outputFile(path_undesignable);
    std::ofstream outputFile_designable(path_designable);
    std::ofstream outputFile_unknown(path_unknown);
    // Check if the file is open
    if (!outputFile.is_open()) {
        std::cerr << "Error opening the file: " << path_undesignable << std::endl;
        return;
    }
    if (!outputFile_designable.is_open()) {
        std::cerr << "Error opening the file: " << path_designable << std::endl;
        return;
    }
    if (!outputFile_unknown.is_open()) {
        std::cerr << "Error opening the file: " << path_unknown << std::endl;
        return;
    }

    std::string path_time = csv + "." + alg + ".time."+getCurrentTimestamp()+".csv";
    #ifdef SPECIAL_HP
    #else
        path_time = csv + "." + alg + ".time.nosh."+getCurrentTimestamp()+".csv";
    #endif
    // Open the file for writing
    std::ofstream timeFile(path_time);

    // Check if the file is open
    if (!timeFile.is_open()) {
        std::cerr << "Error opening the file: " << path_time << std::endl;
        return;
    }else{
        timeFile << "ID,Time(s)" << std::endl;
    }
    std::unordered_map<std::string, GroupY> constr2groupy;
    // print global environment variables
    printf("PATH_UNDESIGNABLE_LIB: %s\n", var_undesignable_lib);
    printf("PATH_DESIGNABLE_LIB: %s\n", var_designable_lib);
    printf("PATH_UNKNOWN_LIB: %s\n", var_unknown_lib);
    if(var_undesignable_lib == NULL){
        std::cerr << "Error: PATH_UNDESIGNABLE_LIB is not set" << std::endl;
        return;
    }
    if(var_designable_lib == NULL){
        std::cerr << "Error: PATH_DESIGNABLE_LIB is not set" << std::endl;
        return;
    }
    if(var_unknown_lib == NULL){
        std::cerr << "Error: PATH_UNKNOWN_LIB is not set" << std::endl;
        return;
    }
    // std::set<std::string> uniq_ud = loadLib(path_undesignable_lib); // undesignable motifs
    // std::set<std::string> uniq_ds = loadLib(path_designable_lib);   // designable   motifs
    // std::set<std::string> uniq_unkown;
    // if (var_unknown_lib != NULL)
    //     uniq_unkown = loadLib(var_unknown_lib); // load unknown motifs if the library is set
    std::cout<<"size of designable motifs: "<<uniq_ds.size()<<std::endl;
    std::cout<<"size of undesignable motifs: "<<uniq_ud.size()<<std::endl;
    std::string path_undesignable_enum = "data/motifs_maxlen14_no_external/results.uniq.json";
    enrichLib(path_undesignable_enum, uniq_ud);
    std::cout<<"size of undesignable motifs after enrichment: "<<uniq_ud.size()<<std::endl;
    std::string path_undesignable_external = "data/undes_motifs_exloop.json";
    enrichLib(path_undesignable_external, uniq_ud);
    std::cout<<"size of undesignable motifs after enrichment (motifs with external loops): "<<uniq_ud.size()<<std::endl;

    std::string path_designable_external_enum = "data/motifs_maxlen14_with_external/short14_external_designable_dg0.txt";
    // add each line of path_designable_external_enum to the designable library
    loadDesignableEnum(path_designable_external_enum, uniq_ds);
    std::cout<<"size of designable motifs after enrichment (motifs with external loops): "<<uniq_ds.size()<<std::endl;
    std::string path_designable_noexternal_enum = "data/motifs_maxlen14_no_external/short14_designable_dg0.txt";
    loadDesignableEnum(path_designable_noexternal_enum, uniq_ds);
    std::cout<<"size of designable motifs after enrichment (motifs without external loops): "<<uniq_ds.size()<<std::endl;

    size_t count_motifs = 0;
    size_t count_motifs_designable = 0;
    size_t count_motifs_undesignable = 0;
    size_t count_motifs_unknown = 0;
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
                    std::cout<<"treestr: "<<treestr<<std::endl;
                    if(uniq_ds.find(treestr) != uniq_ds.end()){
                        std::cout<<"already designable!"<<std::endl;
                        result = "designable";
                    }else if(uniq_ud.find(treestr) != uniq_ud.end()){
                        result = "undesignable";
                        std::cout<<"recur lc.constr: "<<lc.constr<<std::endl;
                        std::cout<<"recur lc.constr: "<<target<<std::endl;
                        std::cout<<"recur   treestr: "<<treestr<<std::endl;
                    }else if(uniq_unknown.find(treestr) != uniq_unknown.end()){
                        result = "unknown";
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
                            // designableLibFile << jstring << std::endl;
                            outputFile_designable << jstring << std::endl;
                            records_designable.push_back(jstring);
                        }
                    }
                    if (result == "undesignable"){
                        std::cout<<"undesignable!"<<std::endl;
                        bool found_ud = false; // check if the motif is already found undesignable
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
                        // if (!found_ud){ // if new undesignable motif, output to lib file
                        //     undesignableLibFile << jstring << std::endl;
                        // }
                    }else{ // unknown case
                        // if new unknown motif, output to file
                        if(uniq_unknown.find(treestr) == uniq_unknown.end()){
                            auto js = jsrecords(lc, y_star, y_sub, y_rivals, puzzle_id);
                            std::string jstring = js.dump();
                            outputFile_unknown << jstring << std::endl;
                        }
                        // add to the unknown library
                        uniq_unknown.insert(treestr);
                        for(Node* rotree: tree->rotated(0)){
                            std::string rotreestr = rotree->toDotBracket();
                            uniq_unknown.insert(rotreestr);
                            delete rotree;
                        }
                    }
                    std::cout<<"designable internal pairs: "<<std::endl;
                    for(auto pair: ds_ipairs)
                        std::cout<<pair<<"  ";
                    printf("\n");
                    delete tree;
                }
                auto end_time = std::chrono::high_resolution_clock::now();
                const std::chrono::duration<double, std::milli> time_ms_y = end_time - start_time;
                float time_seconds_lc = std::chrono::duration_cast<std::chrono::duration<float>>(time_ms_y).count();
                printf("time cost for whole structure: %.4f seconds\n", time_seconds_lc);
                timeFile << puzzle_id << "," << time_seconds_lc <<std::endl;
            }
            if (alg == "getunknown"){
                // remove outputFile, outputFile_designable, and timeFile
                if (std::remove(path_undesignable.c_str()) == 0) {
                    std::cout << "File removed successfully.\n";
                } else {
                    std::perror("Error deleting file");
                }
                if (std::remove(path_designable.c_str()) == 0) {
                    std::cout << "File removed successfully.\n";
                } else {
                    std::perror("Error deleting file");
                }
                if (std::remove(path_time.c_str()) == 0) {
                    std::cout << "File removed successfully.\n";
                } else {
                    std::perror("Error deleting file");
                }
                
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
                // count_motifs += lc_list.size();
                for (auto lc: lc_list){
                    auto start_time_lc = std::chrono::high_resolution_clock::now();
                    lc.printLoopLens();
                    if (lc.hasLongLoop()){
                        printf("the loop exceeds length limit!");
                        continue;
                    }
                    count_motifs++;
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
                        count_motifs_designable++;
                    }else if(uniq_ud.find(treestr) != uniq_ud.end()){
                        result = "undesignable";
                        count_motifs_undesignable++;
                        std::cout<<"recur lc.constr: "<<lc.constr<<std::endl;
                        // std::cout<<"recur    groupy: "<<constr2groupy[lc.constr].constr<<std::endl;
                        std::cout<<"recur lc.constr: "<<target<<std::endl;
                        std::cout<<"recur   treestr: "<<treestr<<std::endl;
                        // std::cout<<"recur     ystar: "<<constr2groupy[lc.constr].star<<std::endl;
                    }else{
                        result = "unknown";
                    }
                    if (result == "unknown")
                    { // unknown case
                        count_motifs_unknown++;
                        // get js
                        auto js = jsrecords(lc, y_star, y_sub, y_rivals, puzzle_id);
                        std::string jstring = js.dump();
                        outputFile_unknown << jstring << std::endl;
                        // uniq_unkown.insert(treestr);
                        // for(Node* rotree: tree->rotated(0)){
                        //     std::string rotreestr = rotree->toDotBracket();
                        //     uniq_unkown.insert(rotreestr);
                        //     delete rotree;
                        // }
                    }
                    for(auto pair: ds_ipairs)
                        std::cout<<pair<<"  ";
                    delete tree;
                    printf("\n");
                }
            }
        }
    }
    // flush the standard output
    std::cout.flush();
    for (auto r: records)
        std::cout<<r<<std::endl;
    std::cout<< "count designable: " << count_designable << std::endl;
    std::cout << "designable records: " << records_designable.size() << std::endl;
    // Close the file
    outputFile.close();
    designableLibFile.close();
    outputFile_unknown.close();
    timeFile.close();
    std::cout << "Strings written to file: " << path_undesignable << std::endl;
    if (alg == "getunknown"){
        cout << "count all motifs: " << count_motifs << std::endl;
        cout << "count designable: " << count_motifs_designable << std::endl;
        cout << "count undesignable: " << count_motifs_undesignable << std::endl;
        cout << "count unknown: " << count_motifs_unknown << std::endl;
        // print the outputFile_unknown path
        std::cout << "Strings written to file: " << path_unknown << std::endl;
        outputFile_unknown.close();
        assert(count_motifs == count_motifs_designable + count_motifs_undesignable + count_motifs_unknown);
    }
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
    std::string path_undesignable = path_prefix + "." + alg + ".log."+time_stamp+".txt";
    #ifdef SPECIAL_HP
    #else
        path_undesignable = path_prefix + "." + alg + ".log.nosh."+time_stamp+".txt";
    #endif

    // Output file stream
    std::ofstream outputFile(path_undesignable);
    // Check if the file is open
    if (!outputFile.is_open()) {
        std::cerr << "Error opening the file: " << path_undesignable << std::endl;
        return;
    }
    // Library files for designable and undesignable motifs
    // if (!designableLibFile.is_open()) {
    //     std::cerr << "Error opening the file: " << "lib_designable.txt" << std::endl;
    //     return;
    // }
    // if (!undesignableLibFile.is_open()) {
    //     std::cerr << "Error opening the file: " << "lib_undesignable.txt" << std::endl;
    //     return;
    // }

    std::string path_time = path_prefix + "." + alg + ".time."+getCurrentTimestamp()+".csv";
    #ifdef SPECIAL_HP
    #else
        path_time = path_prefix + "." + alg + ".time.nosh."+getCurrentTimestamp()+".csv";
    #endif
    // Open the file for writing
    std::ofstream timeFile(path_time);

    // Check if the file is open
    if (!timeFile.is_open()) {
        std::cerr << "Error opening the file: " << path_time << std::endl;
        return;
    }else{
        timeFile << "ID,Time(s)" << std::endl;
    }
    std::unordered_map<std::string, GroupY> constr2groupy;
    std::cout << "[ProgressInfo] Start loading motif libs ..." << std::endl;
    auto time_start_load_lib = std::chrono::high_resolution_clock::now();
    // const char*  var_undesignable_lib = std::getenv("PATH_UNDESIGNABLE_LIB");
    // const char*  var_designable_lib = std::getenv("PATH_DESIGNABLE_LIB");
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
    std::cout<<"size of designable motifs: "<<uniq_ds.size()<<std::endl;
    std::cout<<"size of undesignable motifs: "<<uniq_ud.size()<<std::endl;
    // MotifLib motiflib_ud(path_undesignable_lib);
    // MotifLib motiflib_ds(path_designable_lib);
    // std::cout<<"undesignable motiflib size: "<<motiflib_ud.count_unique<<std::endl;
    // std::cout<<"designable motiflib size: "<<motiflib_ds.count_unique<<std::endl;
    auto time_end_load_lib = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double, std::milli> time_ms_load_lib = time_end_load_lib - time_start_load_lib;
    std::cout << "[ProgressInfo] Time cost for loading motif libs: " << time_ms_load_lib.count()/1000.f << " seconds" << std::endl;
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
                    if(uniq_ds.count(treestr)){
                        std::cout<<"already designable!"<<std::endl;
                        result = "designable";
                    }else if(uniq_ud.count(treestr)){
                        result = "undesignable";
                        std::cout<<"recur lc.constr: "<<lc.constr<<std::endl;
                        // std::cout<<"recur    groupy: "<<constr2groupy[lc.constr].constr<<std::endl;
                        std::cout<<"recur lc.constr: "<<target<<std::endl;
                        std::cout<<"recur   treestr: "<<treestr<<std::endl;
                        // std::cout<<"recur     ystar: "<<constr2groupy[lc.constr].star<<std::endl;
                    }else{ // check designability of the motif
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
                            if(ud_ipairs.find(pairs2string(ipairs))!=ud_ipairs.end()){ // check if any sub-motif is already found undesignable
                                ud = true;
                                break;
                            }
                        }
                        if(ud)  // skip the motif if it is already found undesignable by undesignable sub-motif
                            continue;
                        try{
                            std::cout<<"UMFE log: before alg 5 "<<std::endl;
                            result = alg_5_helper_v2(target, ref_lc, constr_lc, subseq, verbose, dangle); // attempting to identify rival motif(s)
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
                        if(uniq_ds.count(treestr)){
                            found_ds = true;
                            std::cout<<"already found designable!"<<std::endl;
                        }else{
                            uniq_ds.insert(treestr);
                            auto end_time_lc = std::chrono::high_resolution_clock::now();
                            const std::chrono::duration<double, std::milli> time_ms = end_time_lc - start_time_lc;
                            float time_seconds = std::chrono::duration_cast<std::chrono::duration<float>>(time_ms).count();
                            printf("time cost: %.4f seconds\n", time_seconds);
                            auto js = jsrecords(lc, y_star, y_sub, y_rivals, puzzle_id);
                            js["time"] = time_seconds;
                            js["seed"] = SEED_RAND;
                            js["is_duplicated"] = found_ds;
                            js["dot-bracket"] = treestr;
                            std::string jstring = js.dump();
                            designableLibFile << jstring << std::endl;
                            records_designable.push_back(jstring);
                        }
                    }
                    if (result == "undesignable"){
                        std::cout<<"undesignable!"<<std::endl;
                        auto end_time_lc = std::chrono::high_resolution_clock::now();
                        const std::chrono::duration<double, std::milli> time_ms = end_time_lc - start_time_lc;
                        float time_seconds = std::chrono::duration_cast<std::chrono::duration<float>>(time_ms).count();
                        printf("time cost: %.4f seconds\n", time_seconds);
                        bool found_ud = false; // check if the motif is already found undesignable
                        if(uniq_ud.count(treestr)){
                            found_ud = true;
                        }else{
                            uniq_ud.insert(treestr);
                        }
                        auto js = jsrecords(lc, y_star, y_sub, y_rivals, puzzle_id);
                        js["time"] = time_seconds;
                        js["seed"] = SEED_RAND;
                        js["is_duplicated"] = found_ud;
                        // js["id_uniq"] = motiflib_ud.getID(treestr);
                        js["length"] = countDotBrackets(treestr);
                        js["cardinality"] = lc.ps_inside.size() + 1;
                        bool has_external = false;
                        if (lc.ps_outside[0].first == -1){
                            has_external = true;
                        }
                        js["has_external"] = has_external;
                        std::vector<std::string> treestr_all;
                        treestr_all.push_back(treestr);
                        for(Node* rotree: tree->rotated(0)){
                            treestr_all.push_back(rotree->toDotBracket());
                        }
                        js["dot-bracket"] = treestr_all;
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
                        std::cout << "[ProgressInfo] " << "found undesignable motif: " + treestr << std::endl;
                        std::cout << "[ProgressInfo] " << "time cost: " + fl2str(time_seconds) << " seconds" << std::endl;
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
                printf("[ProgressInfo] time cost for whole structure: %.4f seconds\n", time_seconds_lc);
                timeFile << puzzle_id << "," << time_seconds_lc <<std::endl;
            }
        }
    }
    std::cout<<"count undesignable: "<<records.size()<<std::endl;
    for (auto r: records)
        std::cout<<r<<std::endl;
    outputFile.close();
    designableLibFile.close();
    std::cout << "Strings written to file: " << path_undesignable << std::endl;
    if (PLOT && records.size() > 0){
        const char* PATH_FASTMOTIF_PTR = std::getenv("PATH_FASTMOTIF");
        if(PATH_FASTMOTIF_PTR == NULL){
            std::cerr << "Error: PATH_FASTMOTIF is not set" << std::endl;
            return;
        }else{
            std::cout << "PATH_FASTMOTIF: " << PATH_FASTMOTIF_PTR << std::endl;
        }
        std::string PATH_FASTMOTIF(PATH_FASTMOTIF_PTR);
        std::string cmd_str = PATH_FASTMOTIF + "/scripts/parser.py -m y -p " + path_undesignable;
        std::cout<<"extracting plotstr: "<<cmd_str<<std::endl;
        const char* cmd_cstr = cmd_str.c_str();
        std::string path_plotstr = exec_command(cmd_cstr);
        path_plotstr = path_plotstr.substr(0, path_plotstr.size()-1);
        cout << path_plotstr << endl;
        std::ifstream plotFile(path_plotstr);
        std::string line;
        // Check if the file exists and is open
        int retry_count = 5;
        while (!plotFile.is_open() && retry_count > 0) {
            std::cerr << "Error opening the file: " << path_plotstr << ". Retrying..." << std::endl;
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
            plotFile.open(path_plotstr);
            retry_count--;
        }
        if (!plotFile.is_open()) {
            std::cerr << "Error opening the file after retries: " << path_plotstr << std::endl;
            return;
        }
        std::vector<std::string> path_plots;
        while (std::getline(plotFile, line)) {
            // std::cout << line << std::endl;
            std::string cmd_draw = PATH_FASTMOTIF + "/scripts/draw_motif.sh " + line;
            std::cout << "[ProgressInfo] drawing motif: " << line << std::endl;
            const char* cmd_draw_cstr = cmd_draw.c_str();
            std::string output_draw = exec_command(cmd_draw_cstr);
            std::istringstream iss(output_draw);
            std::string last_line;
            std::string line;
            while (std::getline(iss, line)) {
                last_line = line;
            }
            std::cout << "[ProgressInfo] " << last_line << std::endl;
            path_plots.push_back(last_line);
        }
        plotFile.close();
        std::cout<<"-----------------"<<std::endl;
        for(auto path_plot: path_plots){
            std::cout<<path_plot<<std::endl;
            // convert pdf to svg
            // cmd_str = "inkscape --without-gui --file=" + path_plot + " --export-plain-svg=" + path_plot + ".svg";
            // std::cout<<"[ProgressInfo] converting to svg: "<<cmd_str<<std::endl;
            // const char* cmd_cstr2 = cmd_str.c_str();
            // try{
            //     std::string output_svg = exec_command(cmd_cstr2);
            //     std::cout<< "[ProgressInfo] " << path_plot + ".svg" <<std::endl;
            // }catch(const std::exception& e)
            // {
            //     std::cerr << e.what() << '\n';
            //     std::cout << "inkscape error!" << std::endl;
            // }
        }
    }
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
    ("p,plot", "Plot motifs", cxxopts::value<bool>()->default_value("false"))
    ("v,verbose", "Verbose output", cxxopts::value<bool>()->default_value("false"))
    ("max_enum", "max number of enumeration", cxxopts::value<long int>()->default_value("10000000000"))
    ("max_constraint", "max number of constraints", cxxopts::value<long int>()->default_value("100000"))
    ("n_sample", "sampling size ", cxxopts::value<long int>()->default_value("500"))
    ("max_rival", "max number of rival structures or motifs", cxxopts::value<long int>()->default_value("100"));


    auto result = options.parse(argc, argv);
    std::string alg = result["alg"].as<std::string>();
    std::string csv = result["csv"].as<std::string>();
    std::string txt = result["txt"].as<std::string>();
    bool vrna = result["vienna"].as<bool>();
    verbose = result["verbose"].as<bool>();
    PLOT = result["plot"].as<bool>();
    dangle = result["dangle"].as<int>();
    SEED_RAND = result["seed"].as<int>();
    printf("alg: %s, vienna: %d, verbose: %d, dangle: %d\n", alg.c_str(), vrna, verbose, dangle);
    printf("random seed for target initialization: %d\n", SEED_RAND);
    show_configuration();

    MAX_ENUM = result["max_enum"].as<long int>();
    MAX_CONSTRAINT = result["max_constraint"].as<long int>();
    N_SAMPLE = result["n_sample"].as<long int>();
    MAX_RIVAL = result["max_rival"].as<long int>();

    int num_threads = std::thread::hardware_concurrency();
    std::cout<<"number of threads: "<<num_threads<<std::endl;

    if (alg == "0"){
        std::cout<<"no alg was selected!"<<std::endl;
        return 0;
    }

    if (!csv.empty()){
        set_globals();
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
    }else if (alg == "subopt"){ /* suboptimal folding */
        std::cout << alg << std::endl;
        int beamsize = 0;
        bool sharpturn = false;
        float energy_delta = 0.;
        std::string seq;
        std::string constr;
        std::string energy_delta_str;
        while (std::getline(std::cin, seq))
        {
            std::getline(std::cin, constr);
            // input energy delta
            std::getline(std::cin, energy_delta_str);
            // replace x with .; replace . with ?
            std::replace(constr.begin(), constr.end(), 'x', '.');
            std::replace(constr.begin(), constr.end(), '.', '?');
            energy_delta = std::stof(energy_delta_str);
            std::cout<<"constr: "<<constr<<std::endl;
            std::cout<<"energy_delta: "<<energy_delta<<std::endl;
            std::vector<std::string> refs;
            refs =  subopt_fold(seq, constr, energy_delta, beamsize, sharpturn, verbose, dangle);
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
            std::string result = alg_2_helper(ref1, ref2, seq, verbose, dangle);
            auto end_time = std::chrono::high_resolution_clock::now();
            const std::chrono::duration<double, std::milli> time_ms = end_time - start_time;
            printf("alg2 time: %.4f seconds\n", time_ms/1000.f);
            std::cout<< "result: " << result << std::endl;
            // create json object
            json js;
            js["target"] = ref1;
            js["result"] = "unkown";
            std::cout << "number of rival structures: " << y_rivals.size() << std::endl;
            // if UMFE is a substring of result, get the sequence x (after UMFE:)
            if (result.find("UMFE") != std::string::npos) {
                js["result"] = "designable";
                js["umfe"] = result.substr(result.find("UMFE:")+5);
            } else if (result.find("undesignable") != std::string::npos) {
                js["result"] = "undesignable";
            }
            if (mfe_seq_for_alg2 != ""){
                js["mfe"] = mfe_seq_for_alg2;
            }
            js["time"] = time_ms.count()/1000.f;
            std::cout << js.dump() << std::endl;
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
    }else if(alg == "prbound"){  /* get the differential positions between a target structure and multiple rival structures, and obtain a probability bound */
        std::string seq;
        std::string y_target;
        while(std::getline(std::cin, y_target)){
            seq = tg_init(y_target);
            int num_rivals = 1;
            // a vector of rival structures
            std::vector<std::string> y_rivals;
            // get target from input
            // std::getline(std::cin, y_target);
            std::cin >> num_rivals;
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            y_rivals.resize(num_rivals);
            for(int i = 0; i < num_rivals; ++i){
                std::cin >> y_rivals[i];
            }
            // deduplicate rival structures
            std::set<std::string> y_rival_set;
            for(const auto& rival: y_rivals){
                y_rival_set.insert(rival);
            }
            y_rivals.assign(y_rival_set.begin(), y_rival_set.end());
            std::cout<<"number of unique rival structures: "<<y_rivals.size()<<std::endl;
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            // get diffrential positions between the target and each rival, then obtain the union
            std::set<int> differential_positions;
            std::vector<std::vector<std::vector<int>>> cr_loops_vector;
            for(const auto& rival: y_rivals){
                std::set<int> temp_positions;
                auto cr_loops = find_critical_plus(y_target, rival, temp_positions, true);
                cr_loops_vector.push_back(cr_loops);
                differential_positions.insert(temp_positions.begin(), temp_positions.end());
            }
            std::cout<<"differential positions: "<<std::endl;
            for(auto pos: differential_positions){
                std::cout<<pos<<"\t";
            }
            std::cout<<std::endl;
            // calculate number of enumerations for these positions: 4^unpaired * 6^pairs
            std::vector<std::tuple<int, int>> pairs_diff = idx2pair(differential_positions, y_target);
            std::cout<<"differential pairs: "<<std::endl;
            for(auto& pair: pairs_diff)
                std::cout<<std::get<0>(pair)<<"\t"<<std::get<1>(pair)<<std::endl;
            ulong n_enum = count_enum(pairs_diff);
            std::cout<<"enumeration size: "<<n_enum<<std::endl;
            // call alg2_ensemble
            double prob_bound = alg2_ensemble(y_target, y_rivals, cr_loops_vector, pairs_diff, seq, verbose, dangle);
            std::cout<<"probability bound: "<<prob_bound<<std::endl;
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
    }else if (alg == "foldv"){ /* RNA folding via ViennaRNA */
        std::string seq;
        // Read input line by line until EOF (end of file) is reached
        while (std::getline(std::cin, seq)) {
                std::string mfe = fold_vienna(seq);
                std::cout<<seq<<std::endl;
                std::cout<<mfe<<std::endl;
            }
    }else if (alg == "fastmotif"){ /* loops evaluation  */
        var_unknown_lib = NULL;
        set_globals();
        std::string y;
        // Read input line by line until EOF (end of file) is reached
        while (std::getline(std::cin, y)) {
                std::cout<<"please input a structure ->: "<<y<<std::endl;
                std::cout<<y<<std::endl;
                // generate a time stamp used for id
                std::string puzzle_id = "online_" + getCurrentTimestamp() + getRandomIntString();
                std::cout<<"[ProgressInfo] generate a structure id: "<<puzzle_id<<std::endl;
                online_process(y, puzzle_id);
            }
    }else if (alg == "motif"){ /* motif evaluation  */
        MAX_ENUM = 20000000000;
        MAX_CONSTRAINT = 40000000;
        N_SAMPLE = 20000;
        std::string dotbracket;
        std::string seq;
        std::string target;
        std::string ref;
        std::string constr;
        // Read input line by line until EOF (end of file) is reached
        while (std::getline(std::cin, dotbracket)) {
            printf("input motif: %s\n", dotbracket.c_str());
            target = dotbracket2target(dotbracket);
            constr = dotbracket2constraint(dotbracket);
            seq = tg_init(target);
            ref = constr;
            std::replace(ref.begin(), ref.end(), '?', '.');
            printf(" count: %d\n", countOccurrences(constr, '?'));
            printf("target: %s\n", target.c_str());
            printf("constr: %s\n", constr.c_str());
            printf("   ref: %s\n", ref.c_str());
            assert(target.length() == constr.length());
            // TreeNode* root = parseStringToTree(target);
            // int max_internal = max_single(root);
            // if(max_internal > 30){        
            //     std::cout<<"the internal loop is too long: "<<max_internal<<std::endl;            
            //     continue;
            // }
            auto time_start = std::chrono::high_resolution_clock::now();
            std::string result = alg_5_helper_v2(target, ref, constr, seq, verbose, dangle);
            if (result == "UMFE"){
                result = "designable";
            }
            auto time_end = std::chrono::high_resolution_clock::now();
            const std::chrono::duration<double, std::milli> time_ms = time_end - time_start;
            float time_seconds = std::chrono::duration_cast<std::chrono::duration<float>>(time_ms).count();
            // printf("input motif: %s\n", dotbracket.c_str());
            std::cout<<"output  motif: "<<dotbracket<<std::endl;
            std::cout<<"output target: "<<target<<std::endl;
            std::cout<<"output constr: "<<constr<<std::endl;
            printf("time cost: %.4f seconds\n", time_seconds);
            std::cout<<"output result: "<<result<<std::endl;
            if (result == "undesignable"){
                std::cout<<"output number of rival structures: "<<y_rivals.size()<<std::endl;
                for(auto y: y_rivals)
                    std::cout<<"output  rival: "<<y<<std::endl;
            }else if(result == "designable"){
                std::cout<<"output x_umfe: "<<x_umfe<<std::endl;
            }
        }
    }else if (alg == "n1" || alg == "n2" || alg == "n3"){ /* edges evaluation  */
        std::unordered_map<std::string, std::string> struct2seq = loadlib_eterna("data/eterna_umfe_unsolved.csv");
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
    }
    return 0;
}
