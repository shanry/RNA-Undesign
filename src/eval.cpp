/*
 *eval.cpp*
 Evaluate the energy of a given RNA structure. Recording each energy unit.

 author: Wei Yu Tang, Tianshuo Zhou (based on code by He Zhang)
 edited by: 09/2023
*/

#include <fstream>
#include <iostream>
#include <sys/time.h>
#include <stack>
#include <tuple>
#include <cassert>
#include <unordered_map>
#include <algorithm>
#include <string>
#include <map>
#include <set>

#include "LinearFold.h"

#include "Utils/utility_v.h"
#include "Utils/utility_v_max.h"

#include <cstdio>
#include <cassert>
#include <stack>
#include <string>
#include <vector>
#include <map>
#include <tuple>
#include <unordered_map>
#include <set>
#include <cmath>

#define BASE 1000
#define SPECIAL_HP
#define SPECIAL_HP_3
using namespace std;

enum loops {
    hairpin,
    stacking,
    bulge,
    interior,
    multi_inside,
    multi_outside,
    external
};

struct hash_tuple {
  
    template <class T1, class T2, class T3>
  
    size_t operator()(
        const tuple<T1, T2, T3>& x)
        const
    {
        return get<0>(x)*pow(BASE,0)
               + get<1>(x)*pow(BASE,1)
               + get<2>(x)*pow(BASE,2);
    }
};

struct hash_tuple2 {
  
    template <class T1, class T2, class T3, class T4>
  
    size_t operator()(
        const tuple<T1, T2, T3, T4>& x)
        const
    {
        return get<0>(x)*pow(BASE,0)
               + get<1>(x)*pow(BASE,1)
               + get<2>(x)*pow(BASE,2)
               + get<3> (x)*pow(BASE,3);
    }
};

vector<vector<int>> find_critical(string ref1, string ref2, bool is_verbose) {
    assert(ref1.size() == ref2.size());
    int n = ref1.length();

    // map[{i, j, type of loops}] = {y or y', indices}
    unordered_map<tuple<int, int, loops>, pair<int, vector<int>>, hash_tuple> critical_loops;
    unordered_map<tuple<int, int, int, int>, pair<int, vector<int>>, hash_tuple2> critical_bulge;
    unordered_map<tuple<int, int, int, int>, pair<int, vector<int>>, hash_tuple2> critical_internal;
    set<int> critical_positions;

    unordered_map<int, vector<pair<int, int>>> inside;

    stack<pair<int, int>> stk; // tuple of (index, page)
    tuple<int, int> inner_loop;

    for (int t = 0; t < 2; t++) {
        string& ref = (t == 0? ref2 : ref1);
        inside.clear();

        for (int j = 0; j < n; j++) {
            if (ref[j] == '(') {
                if (!stk.empty()) { // +1 for outer loop page
                    stk.top().second ++;
                }
                stk.push(make_pair(j, 0)); // init page=0
            }

            else if (ref[j] == ')') {
                assert(!stk.empty());
                tuple<int, int> top = stk.top();
                int i = get<0>(top), page = get<1>(top);
                stk.pop();

                if (page == 0) { // hairpin
                    tuple<int, int, loops> loop = make_tuple(i, j, hairpin);
                    if (critical_loops.find(loop) != critical_loops.end()) {
                        critical_loops.erase(loop);
                    } else {
                        bool set_hairpin = false;
#ifdef SPECIAL_HP_3
                        if (j - i - 1 == 3) {
                            critical_loops[loop] = make_pair(t, vector<int> {i, j, i+1, i+2, j-1});
                            set_hairpin = true;
                        }
#endif

#ifdef SPECIAL_HP_4
                        if (j - i - 1 == 4) {
                            critical_loops[loop] = make_pair(t, vector<int> {i, j, i+1, i+2, i+3, j-1});
                            set_hairpin = true;
                        }
#endif

#ifdef SPECIAL_HP_6
                        if (j - i - 1 == 6) {
                            critical_loops[loop] = make_pair(t, vector<int> {i, j, i+1, i+2, i+3, i+4, i+5, j-1});
                            set_hairpin = true;
                        }
#endif

                        if (!set_hairpin)
                            critical_loops[loop] = make_pair(t, vector<int> {i, j, i+1, j-1});
                    }
                }

                else if (page == 1) { //single
                    int p = get<0>(inner_loop), q = get<1>(inner_loop);


                    if (p == i+1 && q == j-1) {
                        // stacking
                        tuple<int, int, loops> loop = make_tuple(i, j, stacking);

                        if (critical_loops.find(loop) != critical_loops.end()) {
                            critical_loops.erase(loop);
                        } else {
                            critical_loops[loop] = make_pair(t, vector<int> {i, j, p, q});
                        }
                    } else if (p == i+1 || q == j-1) {
                        // bugle
                        tuple<int, int, int, int> loop = make_tuple(i, p, q, j);

                        if (critical_bulge.find(loop) != critical_bulge.end()) {
                            critical_bulge.erase(loop);
                        } else {
                            critical_bulge[loop] = make_pair(t, vector<int> {i, j, p, q});
                        }
                    } else {
                        // internal
                        tuple<int, int, int, int> loop = make_tuple(i, p, q, j);

                        if (critical_internal.find(loop) != critical_internal.end()) {
                            critical_internal.erase(loop);
                        } else {
                            critical_internal[loop] = make_pair(t, vector<int> {i, j, p, q, i+1, j-1, p-1, q+1});
                        }
                    }
                }

                else { // multi
                    tuple<int, int, loops> loop = make_tuple(i, j, multi_outside);

                    if (critical_loops.find(loop) != critical_loops.end()) {
                        critical_loops.erase(loop);
                    } else {
                        critical_loops[loop] = make_pair(t, vector<int> {i, j, i+1, j-1});
                    }

                    for (auto& item: inside[i]) {
                        int p = item.first, q = item.second;
                        tuple<int, int, loops> loop = make_tuple(p, q, multi_inside);

                        if (critical_loops.find(loop) != critical_loops.end()) {
                            critical_loops.erase(loop);
                        } else {
                            critical_loops[loop] = make_pair(t, vector<int> {p, q, p-1, q+1});
                        }
                    }
                }

                //update inner_loop
                inner_loop = make_tuple(i, j);

                // possible M
                if (!stk.empty()){
                    inside[stk.top().first].push_back(make_pair(i, j));
                }

                // external loop
                if (stk.empty()) {
                    tuple<int, int, loops> loop = make_tuple(i, j, external);

                    if (critical_loops.find(loop) != critical_loops.end()) {
                        critical_loops.erase(loop);
                    } else {
                        critical_loops[loop] = make_pair(t, vector<int> {i, j, i-1, j+1});
                    }
                }
            }
        }
    }
    vector<string> loop_names = {"Hairpin", "Stacking", "Bulge", "Internal", "Multi_inside", "Multi_outside", "External"};
    for (auto &item: critical_loops) {
        for (int &x: item.second.second) {
            if (x >= 0 && x < n) // external loops indices may be out of bound
                critical_positions.insert(x);
        }
    }

    for (auto &item: critical_bulge) {
        for (int &x: item.second.second) {
            critical_positions.insert(x);
        }
    }

    for (auto &item: critical_internal) {
        for (int &x: item.second.second) {
            critical_positions.insert(x);
        }
    }
    vector<vector<int>> cr_loops;
    if (is_verbose){
        printf("critical loops start\n");
        printf("%lu\n", critical_loops.size() + critical_internal.size());
    }
    for (auto &item: critical_loops) {
        // print: is ref1, loop type, indices...
        if (is_verbose)
            printf("%d %d ", item.second.first, get<2>(item.first));
        vector<int> indexed_loop = {item.second.first, get<2>(item.first)};
        for (int &x: item.second.second) {
            if (is_verbose)
                printf("%d ", x);
            indexed_loop.push_back(x);
        }
        if (is_verbose)
            printf("\n");
        cr_loops.push_back(indexed_loop);
    }

    for (auto &item: critical_bulge) {
        // print: is ref1, loop type, indices...
        loops type = bulge;
        vector<int> indexed_loop = {item.second.first, type};
        for (int &x: item.second.second) {
            indexed_loop.push_back(x);
        }
        cr_loops.push_back(indexed_loop);
    }

    for (auto &item: critical_internal) {
        // print: is ref1, loop type, indices...
        loops type = interior;
        vector<int> indexed_loop = {item.second.first, type};
        for (int &x: item.second.second) {
            indexed_loop.push_back(x);
        }
        cr_loops.push_back(indexed_loop);
    }

    if (is_verbose) {
        for (auto &item: critical_loops) {
            printf("%s (%d, %d) ref%d : ", loop_names[get<2>(item.first)].c_str(), get<0>(item.first), get<1>(item.first), -item.second.first+2);
            for (int &x: item.second.second) {
                printf("%d ", x);
            }
            printf("\n");
        }

        for (auto &item: critical_internal) {
            printf("Internal (%d, %d), (%d, %d) ref%d : ", get<0>(item.first), get<3>(item.first), get<1>(item.first), get<2>(item.first), -item.second.first+2);
            for (int &x: item.second.second) {
                printf("%d ", x);
            }
            printf("\n");
        }

        printf("critical positions: ");
        for (int x: critical_positions) {
            printf("%d, ", x);
        }
        printf("\n");
    }

    return cr_loops;
}

vector<vector<int>> find_critical_plus(string ref1, string ref2, set<int>& critical_positions, bool is_verbose) {
    assert(ref1.size() == ref2.size());
    int n = ref1.length();

    // map[{i, j, type of loops}] = {y or y', indices}
    unordered_map<tuple<int, int, loops>, pair<int, vector<int>>, hash_tuple> critical_loops;
    unordered_map<tuple<int, int, int, int>, pair<int, vector<int>>, hash_tuple2> critical_bulge;
    unordered_map<tuple<int, int, int, int>, pair<int, vector<int>>, hash_tuple2> critical_internal;

    unordered_map<int, vector<pair<int, int>>> inside;

    stack<pair<int, int>> stk; // tuple of (index, page)
    tuple<int, int> inner_loop;

    for (int t = 0; t < 2; t++) {
        string& ref = (t == 0? ref2 : ref1);
        inside.clear();

        for (int j = 0; j < n; j++) {
            if (ref[j] == '(') {
                if (!stk.empty()) { // +1 for outer loop page
                    stk.top().second ++;
                }
                stk.push(make_pair(j, 0)); // init page=0
            }

            else if (ref[j] == ')') {
                assert(!stk.empty());
                tuple<int, int> top = stk.top();
                int i = get<0>(top), page = get<1>(top);
                stk.pop();

                if (page == 0) { // hairpin
                    tuple<int, int, loops> loop = make_tuple(i, j, hairpin);
                    if (critical_loops.find(loop) != critical_loops.end()) {
                        critical_loops.erase(loop);
                    } else {
                        bool set_hairpin = false;
#ifdef SPECIAL_HP_3
                        if (j - i - 1 == 3) {
                            critical_loops[loop] = make_pair(t, vector<int> {i, j, i+1, i+2, j-1});
                            set_hairpin = true;
                        }
#endif

#ifdef SPECIAL_HP_4
                        if (j - i - 1 == 4) {
                            critical_loops[loop] = make_pair(t, vector<int> {i, j, i+1, i+2, i+3, j-1});
                            set_hairpin = true;
                        }
#endif

#ifdef SPECIAL_HP_6
                        if (j - i - 1 == 6) {
                            critical_loops[loop] = make_pair(t, vector<int> {i, j, i+1, i+2, i+3, i+4, i+5, j-1});
                            set_hairpin = true;
                        }
#endif

                        if (!set_hairpin)
                            critical_loops[loop] = make_pair(t, vector<int> {i, j, i+1, j-1});
                    }
                }

                else if (page == 1) { //single
                    int p = get<0>(inner_loop), q = get<1>(inner_loop);


                    if (p == i+1 && q == j-1) {
                        // stacking
                        tuple<int, int, loops> loop = make_tuple(i, j, stacking);

                        if (critical_loops.find(loop) != critical_loops.end()) {
                            critical_loops.erase(loop);
                        } else {
                            critical_loops[loop] = make_pair(t, vector<int> {i, j, p, q});
                        }
                    } else if (p == i+1 || q == j-1) {
                        // bugle
                        tuple<int, int, int, int> loop = make_tuple(i, p, q, j);

                        if (critical_bulge.find(loop) != critical_bulge.end()) {
                            critical_bulge.erase(loop);
                        } else {
                            critical_bulge[loop] = make_pair(t, vector<int> {i, j, p, q});
                        }
                    } else {
                        // internal
                        tuple<int, int, int, int> loop = make_tuple(i, p, q, j);

                        if (critical_internal.find(loop) != critical_internal.end()) {
                            critical_internal.erase(loop);
                        } else {
                            critical_internal[loop] = make_pair(t, vector<int> {i, j, p, q, i+1, j-1, p-1, q+1});
                        }
                    }
                }

                else { // multi
                    tuple<int, int, loops> loop = make_tuple(i, j, multi_outside);

                    if (critical_loops.find(loop) != critical_loops.end()) {
                        critical_loops.erase(loop);
                    } else {
                        critical_loops[loop] = make_pair(t, vector<int> {i, j, i+1, j-1});
                    }

                    for (auto& item: inside[i]) {
                        int p = item.first, q = item.second;
                        tuple<int, int, loops> loop = make_tuple(p, q, multi_inside);

                        if (critical_loops.find(loop) != critical_loops.end()) {
                            critical_loops.erase(loop);
                        } else {
                            critical_loops[loop] = make_pair(t, vector<int> {p, q, p-1, q+1});
                        }
                    }
                }

                //update inner_loop
                inner_loop = make_tuple(i, j);

                // possible M
                if (!stk.empty()){
                    inside[stk.top().first].push_back(make_pair(i, j));
                }

                // external loop
                if (stk.empty()) {
                    tuple<int, int, loops> loop = make_tuple(i, j, external);

                    if (critical_loops.find(loop) != critical_loops.end()) {
                        critical_loops.erase(loop);
                    } else {
                        critical_loops[loop] = make_pair(t, vector<int> {i, j, i-1, j+1});
                    }
                }
            }
        }
    }
    vector<string> loop_names = {"Hairpin", "Stacking", "Bulge", "Internal", "Multi_inside", "Multi_outside", "External"};
    for (auto &item: critical_loops) {
        for (int &x: item.second.second) {
            if (x >= 0 && x < n) // external loops indices may be out of bound
                critical_positions.insert(x);
        }
    }

    for (auto &item: critical_bulge) {
        for (int &x: item.second.second) {
            critical_positions.insert(x);
        }
    }

    for (auto &item: critical_internal) {
        for (int &x: item.second.second) {
            critical_positions.insert(x);
        }
    }
    vector<vector<int>> cr_loops;
    if(is_verbose){
        printf("critical loops start\n");
        printf("%lu\n", critical_loops.size() + critical_bulge.size() + critical_internal.size());
    }
    for (auto &item: critical_loops) {
        // print: is ref1, loop type, indices...
        if(is_verbose)
            printf("%d %d: ", item.second.first, get<2>(item.first));
        vector<int> indexed_loop = {item.second.first, get<2>(item.first)};
        for (int &x: item.second.second) {
            if(is_verbose)
                printf("%d ", x);
            indexed_loop.push_back(x);
        }
        if(is_verbose)
            printf("\n");
        cr_loops.push_back(indexed_loop);
    }

    for (auto &item: critical_bulge) {
        // print: is ref1, loop type, indices...
        loops type = bulge;
        if(is_verbose)
            printf("%d %d: ", item.second.first, type);
        vector<int> indexed_loop = {item.second.first, type};
        for (int &x: item.second.second) {
            if(is_verbose)
                printf("%d ", x);
            indexed_loop.push_back(x);
        }
        if(is_verbose)
            printf("\n");
        cr_loops.push_back(indexed_loop);
    }

    for (auto &item: critical_internal) {
        // print: is ref1, loop type, indices...
        loops type = interior;
        if(is_verbose)
            printf("%d %d: ", item.second.first, type);
        vector<int> indexed_loop = {item.second.first, type};
        for (int &x: item.second.second) {
            if(is_verbose)
                printf("%d ", x);
            indexed_loop.push_back(x);
        }
        if(is_verbose)
            printf("\n");
        cr_loops.push_back(indexed_loop);
    }

    if (is_verbose) {
        for (auto &item: critical_loops) {
            printf("%s (%d, %d) ref%d : ", loop_names[get<2>(item.first)].c_str(), get<0>(item.first), get<1>(item.first), -item.second.first+2);
            for (int &x: item.second.second) {
                printf("%d ", x);
            }
            printf("\n");
        }

        for (auto &item: critical_bulge) {
            printf("Bulge (%d, %d), (%d, %d) ref%d : ", get<0>(item.first), get<3>(item.first), get<1>(item.first), get<2>(item.first), -item.second.first+2);
            for (int &x: item.second.second) {
                printf("%d ", x);
            }
            printf("\n");
        }

        for (auto &item: critical_internal) {
            printf("Internal (%d, %d), (%d, %d) ref%d : ", get<0>(item.first), get<3>(item.first), get<1>(item.first), get<2>(item.first), -item.second.first+2);
            for (int &x: item.second.second) {
                printf("%d ", x);
            }
            printf("\n");
        }
    }
    if (is_verbose){
        printf("critical loops:\n");
        for(auto cr_loop: cr_loops){
            for(auto value: cr_loop)
                printf("%d\t", value);
            printf("\n");
        }
    }
    if (is_verbose){
        printf("critical positions: ");
        for (int x: critical_positions) {
            printf("%d, ", x);
        }
        printf("\n");
    }
    return cr_loops;
}

long diff_eval(string& seq, vector<vector<int>>& cr_loops, bool& is_verbose, int& dangle_model) {
    int n = seq.length();
    
    vector<int> if_tetraloops;
    vector<int> if_hexaloops;
    vector<int> if_triloops;

    v_init_tetra_hex_tri(seq, n, if_tetraloops, if_hexaloops, if_triloops); // calculate if_tetraloops, if_hexaloops, if_triloops
    
    vector<int> eval_nucs(n);
    for (int i = 0; i < n; i++) {
      eval_nucs[i] = GET_ACGU_NUM_V(seq[i]); // lhuang: explicitly use Vienna coding (not very nice)
    }

    long score = 0;
    long energy_ref1 = 0, energy_ref2 = 0;
    for (auto &item: cr_loops) {
        bool is_ref1 = item[0];
        loops type = (loops) item[1];
        
        int i = item[2], j = item[3];

        int nuci = eval_nucs[i];
        int nucj = eval_nucs[j];
        int nuci1 = (i + 1) < n ? eval_nucs[i + 1] : -1;
        int nucj_1 = (j - 1) > -1 ? eval_nucs[j - 1] : -1;
        int nuci_1 = (i-1>-1) ? eval_nucs[i-1] : -1; // only for calculating v_score_M1
        int nucj1 = (j+1) < n ? eval_nucs[j+1] : -1; // only for calculating v_score_M1

        if (type == hairpin) {
            int tetra_hex_tri = -1;

            if (j-i-1 == 4) // 6:tetra
                tetra_hex_tri = if_tetraloops[i];
            else if (j-i-1 == 6) // 8:hexa
                tetra_hex_tri = if_hexaloops[i];
            else if (j-i-1 == 3) // 5:tri
                tetra_hex_tri = if_triloops[i];

            score = - v_score_special_hairpin(i, j, nuci, nuci1, nucj_1, nucj, tetra_hex_tri);

            if (is_verbose)
                printf("Hairpin loop ( %d, %d) %c%c ref%d: %.2f\n", i, j, seq[i], seq[j], -is_ref1+2, score / -100.0);
        
        } else if (type == stacking || type == bulge || type == interior) {
            int p = item[4], q = item[5];

            int nucp_1 = eval_nucs[p-1], nucp = eval_nucs[p], nucq = eval_nucs[q], nucq1 = eval_nucs[q+1];

            score = - v_score_single(i,j,p,q, nuci, nuci1, nucj_1, nucj,
                                            nucp_1, nucp, nucq, nucq1);
            if (is_verbose)
                printf("Stacking/Bulge loop ( %d, %d) %c%c; ( %d, %d) %c%c ref%d: %.2f\n", i, j, seq[i], seq[j], p, q, seq[p],seq[q], -is_ref1+2, score / -100.0);
       
        } else if (type == multi_inside) {
            score = - v_score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, n, dangle_model);

            if (is_verbose)
                printf("Multi Inside (%d, %d) %c%c ref%d: %.2f \n", i, j, seq[i], seq[j], -is_ref1+2, score / -100.0);
        } else if (type == multi_outside) {
            score = - v_score_multi(i, j, nuci, nuci1, nucj_1, nucj, n, dangle_model);
            
            if (is_verbose)
                printf("Multi Outside (%d, %d) %c%c ref%d: %.2f \n", i, j, seq[i], seq[j], -is_ref1+2, score / -100.0);
        } else if (type == external) {
            int k = i - 1;
            int nuck = k > -1 ? eval_nucs[k] : -1;
            int nuck1 = eval_nucs[k+1];

            score = - v_score_external_paired(k+1, j, nuck, nuck1,
                                              nucj, nucj1, n, dangle_model);

            if (is_verbose)
                printf("External loop (%d, %d) %c%c ref%d: %.2f \n", i, j, seq[i], seq[j], -is_ref1+2, score / -100.0);
        }

        if (is_ref1)
            energy_ref1 += score;
        else
            energy_ref2 += score;
    }
    if (is_verbose)
        printf("ref1 energy = %.2f, ref2 energy = %.2f\n", energy_ref1 / -100.0, energy_ref2 / -100.0);
    return energy_ref1 - energy_ref2;
}

long linear_eval(string& seq, string& ref, bool& is_verbose, int& dangle_model) {

    int seq_length = seq.length();

    vector<int> if_tetraloops;
    vector<int> if_hexaloops;
    vector<int> if_triloops;

    v_init_tetra_hex_tri(seq, seq_length, if_tetraloops, if_hexaloops, if_triloops); // calculate if_tetraloops, if_hexaloops, if_triloops

    vector<int> eval_nucs;
    eval_nucs.clear();
    eval_nucs.resize(seq_length);
    for (int i = 0; i < seq_length; ++i) {
      eval_nucs[i] = GET_ACGU_NUM_V(seq[i]); // lhuang: explicitly use Vienna coding (not very nice)
    }

    long total_energy = 0;
    long external_energy = 0;
    long M1_energy[seq_length];
    long multi_number_unpaired[seq_length];
    // int external_number_unpaired = 0;

    stack<pair<int, int>> stk; // tuple of (index, page)
    tuple<int, int> inner_loop;

    for (int j=0; j<seq_length; j++) {
        M1_energy[j] = 0; // init multi of position j
        multi_number_unpaired[j] = 0;

        if (ref[j] == '.') {
            if (!stk.empty())
                multi_number_unpaired[stk.top().first] += 1;
        }

        else if (ref[j] == '(') {
            if (!stk.empty()) { // +1 for outer loop page
                stk.top().second ++;
            }
            stk.push(make_pair(j, 0)); // init page=0
        }

        else if (ref[j] == ')') {
            assert(!stk.empty());
            tuple<int, int> top = stk.top();
            int i = get<0>(top), page = get<1>(top);
            stk.pop();

            int nuci = eval_nucs[i];
            int nucj = eval_nucs[j];
            int nuci1 = (i + 1) < seq_length ? eval_nucs[i + 1] : -1;
            int nucj_1 = (j - 1) > -1 ? eval_nucs[j - 1] : -1;
            int nuci_1 = (i-1>-1) ? eval_nucs[i-1] : -1; // only for calculating v_score_M1
            int nucj1 = (j+1) < seq_length ? eval_nucs[j+1] : -1; // only for calculating v_score_M1

            if (page == 0) { // hairpin
                int tetra_hex_tri = -1;
                if (j-i-1 == 4) // 6:tetra
                    tetra_hex_tri = if_tetraloops[i];
                else if (j-i-1 == 6) // 8:hexa
                    tetra_hex_tri = if_hexaloops[i];
                else if (j-i-1 == 3) // 5:tri
                    tetra_hex_tri = if_triloops[i];
                
                int newscore = - v_score_special_hairpin(i, j, nuci, nuci1, nucj_1, nucj, tetra_hex_tri);
                if (is_verbose)
                    printf("Hairpin loop ( %d, %d) %c%c : %.2f\n", i+1, j+1, seq[i], seq[j], newscore / -100.0);
                total_energy += newscore;
            }

            else if (page == 1) { //single
                int p = get<0>(inner_loop), q = get<1>(inner_loop);

                int nucp_1 = eval_nucs[p-1], nucp = eval_nucs[p], nucq = eval_nucs[q], nucq1 = eval_nucs[q+1];

                int newscore = - v_score_single(i,j,p,q, nuci, nuci1, nucj_1, nucj,
                                                  nucp_1, nucp, nucq, nucq1);
                if (is_verbose)
                    printf("Interior loop ( %d, %d) %c%c; ( %d, %d) %c%c : %.2f\n", i+1, j+1, seq[i], seq[j], p+1, q+1, seq[p],seq[q], newscore / -100.0);
                total_energy += newscore;
            }

            else { //multi
                int multi_score = 0;
                multi_score += M1_energy[i];
                multi_score += - v_score_multi(i, j, nuci, nuci1, nucj_1, nucj, seq_length, dangle_model);
                multi_score += - v_score_multi_unpaired(i+1, i + multi_number_unpaired[i]); // current model is 0
                if (is_verbose)
                    printf("Multi loop ( %d, %d) %c%c : %.2f\n", i+1, j+1, seq[i], seq[j], multi_score / -100.0);
                total_energy += multi_score;
            }

            //update inner_loop
            inner_loop = make_tuple(i, j);

            // possible M
            if (!stk.empty())
                M1_energy[stk.top().first] += - v_score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length, dangle_model);

            // check if adding external energy
            if (stk.empty()) {
                int k = i - 1;
                int nuck = k > -1 ? eval_nucs[k] : -1;
                int nuck1 = eval_nucs[k+1];
                external_energy +=  - v_score_external_paired(k+1, j, nuck, nuck1,
                                                            nucj, nucj1, seq_length, dangle_model);
                // external_energy += 0; currently external unpaired is 0
            }
        }
    }
    if (is_verbose)
        printf("External loop : %.2f\n", external_energy / -100.0);
    total_energy += external_energy;
    return total_energy;
}

long max_diff(int n, vector<vector<int>>& cr_loops, bool is_verbose, int dangle_model) {
    long score = 0;
    long energy_ref1 = 0, energy_ref2 = 0;

    // vector<int> if_tetraloops;
    // vector<int> if_hexaloops;
    // vector<int> if_triloops;

    // v_init_tetra_hex_tri(seq, seq_length, if_tetraloops, if_hexaloops, if_triloops);

    for (auto &item: cr_loops) {
        bool is_ref1 = item[0];
        loops type = (loops) item[1];

        int i = item[2], j = item[3];

        if (type == hairpin) {
            int tetra_hex_tri = -1;

            // if (j-i-1 == 4) // 6:tetra
            //     tetra_hex_tri = if_tetraloops[i];
            // else if (j-i-1 == 6) // 8:hexa
            //     tetra_hex_tri = if_hexaloops[i];
            // else if (j-i-1 == 3) // 5:tri
            //     tetra_hex_tri = if_triloops[i];

            if (is_ref1)
                score = - v_score_hairpin_max(i, j, tetra_hex_tri);
            else
                score = - v_score_hairpin_min(i, j, tetra_hex_tri);

            if (is_verbose)
                printf("Hairpin loop ( %d, %d) ref%d: %.2f\n", i, j, -is_ref1+2, score / -100.0);
        
        } else if (type == stacking || type == bulge || type == interior) {
            int p = item[4], q = item[5];

            if (is_ref1)
                score = - v_score_single_max(i,j,p,q);
            else
                score = - v_score_single_min(i,j,p,q);

            if (is_verbose)
                printf("Stacking/Bulge loop ( %d, %d) ; ( %d, %d) ref%d: %.2f\n", i, j, p, q, -is_ref1+2, score / -100.0);
       
        } else if (type == multi_inside) {
            if (is_ref1) 
                score = - v_score_M1_max(i, j, j, n, dangle_model);
            else
                score = - v_score_M1_min(i, j, j, n, dangle_model);

            if (is_verbose)
                printf("Multi Inside (%d, %d) ref%d: %.2f \n", i, j, -is_ref1+2, score / -100.0);
        } else if (type == multi_outside) {
            if (is_ref1)
                score = - v_score_multi_max(i, j, n, dangle_model);
            else
                score = - v_score_multi_min(i, j, n, dangle_model);
            
            if (is_verbose)
                printf("Multi Outside (%d, %d) ref%d: %.2f \n", i, j, -is_ref1+2, score / -100.0);
        } else if (type == external) {
            if (is_ref1)
                score = - v_score_external_paired_max(i, j, n, dangle_model);
            else
                score = - v_score_external_paired_min(i, j, n, dangle_model);

            if (is_verbose)
                printf("External loop (%d, %d) ref%d: %.2f \n", i, j, -is_ref1+2, score / -100.0);
        }

        if (is_ref1)
            energy_ref1 += score;
        else
            energy_ref2 += score;
    }

    return energy_ref1 - energy_ref2;
}

bool test_diff(string seq, string ref1, string ref2, bool is_verbose, int dangle_model){
    long energy_ref1 = linear_eval(seq, ref1, is_verbose, dangle_model);
    long energy_ref2 = linear_eval(seq, ref2, is_verbose, dangle_model);
    printf("ref1 energy: %.2f\n", energy_ref1/-100.0);
    printf("ref2 energy: %.2f\n", energy_ref2/-100.0);

    vector<vector<int>> cr_loops = find_critical(ref1, ref2, is_verbose);
    long delta_energy = diff_eval(seq, cr_loops, is_verbose, dangle_model);
    // weiyu : changed diff_eval, need to modify test
    // long delta_energy = diff_eval(seq, ref1, ref2, is_verbose, dangle_model);
    bool equal = (energy_ref1 - energy_ref2)==delta_energy;
    printf("pass test: %s\n", equal ? "true" : "false");
    printf("e1 - e2: %.2f\n", (energy_ref1 - energy_ref2) / -100.0);
    printf("delta  : %.2f\n", delta_energy/-100.0);
    return 1;
}

// Tianshuo: LF code

#ifdef lv
    bool comparefunc(std::pair<int,State> a, std::pair<int,State> b) {
        return a.first > b.first;
    }

    void BeamCKYParser::sort_keys(std::unordered_map<int, State> &map, std::vector<std::pair<int,State>> &sorted_keys) {
        sorted_keys.clear();
        for(auto &kv : map) {
            sorted_keys.push_back(kv);
        }
        sort(sorted_keys.begin(), sorted_keys.end(), comparefunc);    
    }
#endif

void BeamCKYParser::get_parentheses(char* result, string& seq) {
    memset(result, '.', seq_length);
    result[seq_length] = 0;

    stack<tuple<int, int, State>> stk;
    stk.push(make_tuple(0, seq_length-1, bestC[seq_length-1]));

    if(is_verbose){
            printf(">verbose\n");
    }
    // verbose stuff
    vector<pair<int,int>> multi_todo;
    unordered_map<int,int> mbp; // multi bp
    double total_energy = .0;
    double external_energy = .0;

    while ( !stk.empty() ) {
        tuple<int, int, State> top = stk.top();
        int i = get<0>(top), j = get<1>(top);
        State& state = get<2>(top);
        stk.pop();

        int k, p, q;

        switch (state.manner) {
            case MANNER_H:
                // this state should not be traced
                break;
            case MANNER_HAIRPIN:
                {
                    result[i] = '(';
                    result[j] = ')';
                    if(is_verbose) {
                        int tetra_hex_tri = -1;
                        if (j-i-1 == 4) // 6:tetra
                            tetra_hex_tri = if_tetraloops[i];
                        else if (j-i-1 == 6) // 8:hexa
                            tetra_hex_tri = if_hexaloops[i];
                        else if (j-i-1 == 3) // 5:tri
                            tetra_hex_tri = if_triloops[i];
                        int nuci = nucs[i], nucj = nucs[j];
                        int nuci1 = (i + 1) < seq_length ? nucs[i + 1] : -1;
                        int nucj_1 = (j - 1) > -1 ? nucs[j - 1] : -1;

                        value_type newscore = - v_score_hairpin(i, j, nuci, nuci1, nucj_1, nucj, tetra_hex_tri);
                        printf("Hairpin loop ( %d, %d) %c%c : %.2f\n", i+1, j+1, seq[i], seq[j], newscore / -100.0);
                        total_energy += newscore;
                    }
                }
                break;
            case MANNER_SINGLE:
                {
                    result[i] = '(';
                    result[j] = ')';
                    p = i + state.trace.paddings.l1;
                    q = j - state.trace.paddings.l2;
                    stk.push(make_tuple(p, q, bestP[q][p]));
                    if(is_verbose) {
                        int nuci = nucs[i], nuci1 = nucs[i+1], nucj_1 = nucs[j-1], nucj = nucs[j];
                        int nucp_1 = nucs[p-1], nucp = nucs[p], nucq = nucs[q], nucq1 = nucs[q+1];

                        value_type newscore = -v_score_single(i,j,p,q, nuci, nuci1, nucj_1, nucj,
                                                          nucp_1, nucp, nucq, nucq1);
                        printf("Interior loop ( %d, %d) %c%c; ( %d, %d) %c%c : %.2f\n", i+1, j+1, seq[i], seq[j], p+1, q+1, seq[p],seq[q], newscore / -100.0);
                        total_energy += newscore;
                    }
                }
                break;
            case MANNER_HELIX:
                {
                    result[i] = '(';
                    result[j] = ')';
                    stk.push(make_tuple(i+1, j-1, bestP[j-1][i+1]));
                    if(is_verbose){
                        p = i + 1;
                        q = j - 1;
                        int nuci = nucs[i], nuci1 = nucs[i+1], nucj_1 = nucs[j-1], nucj = nucs[j];
                        int nucp_1 = nucs[p-1], nucp = nucs[p], nucq = nucs[q], nucq1 = nucs[q+1];

                        value_type newscore = -v_score_single(i,j,p,q, nuci, nuci1, nucj_1, nucj,
                                                          nucp_1, nucp, nucq, nucq1);
                        printf("Interior loop ( %d, %d) %c%c; ( %d, %d) %c%c : %.2f\n", i+1, j+1, seq[i], seq[j], p+1, q+1, seq[p],seq[q], newscore / -100.0);
                        total_energy += newscore;
                    }
                }
                break;
            case MANNER_MULTI: 
                p = i + state.trace.paddings.l1;
                q = j - state.trace.paddings.l2;
                stk.push(make_tuple(p, q, bestM2[q][p]));
                break;
            case MANNER_MULTI_eq_MULTI_plus_U:
                p = i + state.trace.paddings.l1;
                q = j - state.trace.paddings.l2;
                stk.push(make_tuple(p, q, bestM2[q][p]));
                break;
            case MANNER_P_eq_MULTI:
                result[i] = '(';
                result[j] = ')';
                stk.push(make_tuple(i, j, bestMulti[j][i]));
                if(is_verbose) {
                    multi_todo.push_back(make_pair(i,j));
                }
                break;
            case MANNER_M2_eq_M_plus_P:
                k = state.trace.split;
                stk.push(make_tuple(i, k, bestM[k][i]));
                stk.push(make_tuple(k+1, j, bestP[j][k+1]));
                if(is_verbose)
                    mbp[k+1] = j;
                break;
            case MANNER_M_eq_M2:
                stk.push(make_tuple(i, j, bestM2[j][i]));
                break;
            case MANNER_M_eq_M_plus_U:
                stk.push(make_tuple(i, j-1, bestM[j-1][i]));
                break;
            case MANNER_M_eq_P:
                stk.push(make_tuple(i, j, bestP[j][i]));
                if(is_verbose)
                    mbp[i] = j;
                break;
            case MANNER_C_eq_C_plus_U:
                k = j - 1;
                if (k != -1)
                    stk.push(make_tuple(0, k, bestC[k]));
                if (is_verbose) 
                    external_energy += - v_score_external_unpaired(0, 0); // zero at this moment
                break;
            case MANNER_C_eq_C_plus_P:
                {
                    k = state.trace.split;
                    if (k != -1) {
                        stk.push(make_tuple(0, k, bestC[k]));
                        stk.push(make_tuple(k+1, j, bestP[j][k+1]));
                    }
                    else {
                        stk.push(make_tuple(i, j, bestP[j][i]));
                    }
                    if (is_verbose) {
                        int nuck = k > -1 ? nucs[k] : -1;
                        int nuck1 = nucs[k+1], nucj = nucs[j];
                        int nucj1 = (j + 1) < seq_length ? nucs[j + 1] : -1;
                        external_energy +=  - v_score_external_paired(k+1, j, nuck, nuck1,
                                                                    nucj, nucj1, seq_length, dangle_model);
                    }
                }
                break;
            default:  // MANNER_NONE or other cases
                if (use_constraints){
                    printf("We can't find a valid structure for this sequence and constraint.\n");
                    printf("There are two minor restrictions in our real system:\n");
                    printf("the length of an interior loop is bounded by 30nt \n");
                    printf("(a standard limit found in most existing RNA folding software such as CONTRAfold)\n");
                    printf("so is the leftmost (50-end) unpaired segment of a multiloop (new constraint).\n");
                    exit(1);
                } 
                printf("wrong manner at %d, %d: manner %d\n", i, j, state.manner); fflush(stdout);
                assert(false);
                
        }
    }

    if(is_verbose) {
        for (auto item : multi_todo) {
            int i = item.first;
            int j = item.second;
            int nuci = nucs[i], nuci1 = nucs[i+1], nucj_1 = nucs[j-1], nucj = nucs[j];
            value_type multi_energy = - v_score_multi(i, j, nuci, nuci1, nucj_1, nucj, seq_length, dangle_model);
            int num_unpaired = 0;
            for (int k=i+1; k<j; ++k) {
                if (result[k] == '.')
                    num_unpaired += 1;
                else if (result[k] == '(') {
                    int p = k, q = mbp[k];
                    int nucp_1 = nucs[p-1], nucp = nucs[p], nucq = nucs[q], nucq1 = nucs[q+1];

                    multi_energy += - v_score_M1(p, q, q, nucp_1, nucp, nucq, nucq1, seq_length, dangle_model);
                    k = q;
                }
            }
            multi_energy += - v_score_multi_unpaired(1, num_unpaired);

            printf("Multi loop ( %d, %d) %c%c : %.2f\n", i+1, j+1, seq[i], seq[j], multi_energy / -100.0);
            total_energy += multi_energy;
        }

        printf("External loop : %.2f\n", external_energy / -100.0);
        total_energy += external_energy;

#ifndef lv
            printf("Energy(kcal/mol): %.2f\n", total_energy / -100.0);
#endif
    }

    return;
}

unsigned long quickselect_partition(vector<pair<value_type, int>>& scores, unsigned long lower, unsigned long upper) {
    value_type pivot = scores[upper].first;
    while (lower < upper) {
        while (scores[lower].first < pivot) ++lower;
        while (scores[upper].first > pivot) --upper;
        if (scores[lower].first == scores[upper].first) ++lower;
        else if (lower < upper) swap(scores[lower], scores[upper]);

    }
    return upper;
}

// in-place quick-select
value_type quickselect(vector<pair<value_type, int>>& scores, unsigned long lower, unsigned long upper, unsigned long k) {
    if ( lower == upper ) return scores[lower].first;
    unsigned long split = quickselect_partition(scores, lower, upper);
    unsigned long length = split - lower + 1;
    if (length == k) return scores[split].first;
    else if (k  < length) return quickselect(scores, lower, split-1, k);
    else return quickselect(scores, split+1, upper, k - length);
}

value_type BeamCKYParser::beam_prune(std::unordered_map<int, State> &beamstep) {
    scores.clear();
    for (auto &item : beamstep) {
        int i = item.first;
        State &cand = item.second;
        int k = i - 1;
        value_type newscore;
        // lisiz: for _V, avoid -inf-int=+inf
        if ((k >= 0) && (bestC[k].score == VALUE_MIN)) newscore = VALUE_MIN;
        else newscore = (k >= 0 ? bestC[k].score : 0) + cand.score;
        scores.push_back(make_pair(newscore, i));
    }
    if (scores.size() <= beam) return VALUE_MIN;
    value_type threshold = quickselect(scores, 0, scores.size() - 1, scores.size() - beam);
    for (auto &p : scores) {
        if (p.first < threshold) beamstep.erase(p.second);
    }

    return threshold;
}

void BeamCKYParser::sortM(value_type threshold,
                          std::unordered_map<int, State> &beamstep,
                          std::vector<std::pair<value_type, int>> &sorted_stepM) {
    sorted_stepM.clear();
    if (threshold == VALUE_MIN) {
        // no beam pruning before, so scores vector not usable
        for (auto &item : beamstep) {
            int i = item.first;
            State &cand = item.second;
            int k = i - 1;
            value_type newscore;
            // lisiz: constraints may cause all VALUE_MIN, sorting has no use
            if ((use_constraints) && (k >= 0) && (bestC[k].score == VALUE_MIN)) newscore = cand.score;
            else newscore = (k >= 0 ? bestC[k].score : 0) + cand.score;
            sorted_stepM.push_back(make_pair(newscore, i));
        }
    } else {
        for (auto &p : scores) {
            if (p.first >= threshold) sorted_stepM.push_back(p);
        }
    }

    sort(sorted_stepM.begin(), sorted_stepM.end(), std::greater<pair<value_type, int>>());
}

//subopt zuker inside, bestP_beta sorted, subopt zuker optimized version
std::string BeamCKYParser::get_parentheses_inside_real_backtrace(int i, int j, State& state, map<tuple<BestTypes, int, int>, string>& global_visited_inside, set<pair<int,int> >& window_visited) {

    Manner manner = state.manner;
    value_type score = state.score;

    if (manner == MANNER_H)
        return "";
    else if (manner == MANNER_HAIRPIN){
        if (global_visited_inside.count(make_tuple(TYPE_P, i, j))) 
            return global_visited_inside[make_tuple(TYPE_P, i, j)];
        std::string inner(j-i-1, '.');
        global_visited_inside[make_tuple(TYPE_P, i, j)] = '(' + inner + ')';
        window_fill(window_visited, i, j, seq_length, window_size);
        return global_visited_inside[make_tuple(TYPE_P, i, j)];
    }

    else if (manner == MANNER_SINGLE){
        if (global_visited_inside.count(make_tuple(TYPE_P,i,j)))
            return global_visited_inside[make_tuple(TYPE_P,i,j)];
        auto p = i + state.trace.paddings.l1;
        auto q = j - state.trace.paddings.l2;
        std::string left(state.trace.paddings.l1-1, '.');
        std::string right(state.trace.paddings.l2-1,'.');

        if (!global_visited_inside.count(make_tuple(TYPE_P,p,q)))
            global_visited_inside[make_tuple(TYPE_P,p,q)] = get_parentheses_inside_real_backtrace(p, q, bestP[q][p],global_visited_inside, window_visited);

        global_visited_inside[make_tuple(TYPE_P,i,j)] = '(' + left + global_visited_inside[make_tuple(TYPE_P,p,q)] + right + ')';
        window_fill(window_visited, i, j, seq_length, window_size);
        return global_visited_inside[make_tuple(TYPE_P,i,j)];
    }

    else if (manner == MANNER_HELIX){
        if (global_visited_inside.count(make_tuple(TYPE_P,i,j)))
            return global_visited_inside[make_tuple(TYPE_P,i,j)];
        if (!global_visited_inside.count(make_tuple(TYPE_P,i+1,j-1)))
            global_visited_inside[make_tuple(TYPE_P,i+1,j-1)] = get_parentheses_inside_real_backtrace(i+1, j-1, bestP[j-1][i+1],global_visited_inside, window_visited);
        global_visited_inside[make_tuple(TYPE_P,i,j)] = '(' + global_visited_inside[make_tuple(TYPE_P,i+1,j-1)] + ')';
        window_fill(window_visited, i, j, seq_length, window_size);
        return global_visited_inside[make_tuple(TYPE_P,i,j)];
    }

    else if (manner == MANNER_MULTI){
        // Multi = (....M2.....)  Multi: M; M2: N; M1: m
        if (global_visited_inside.count(make_tuple(TYPE_Multi,i,j)))
            return global_visited_inside[make_tuple(TYPE_Multi,i,j)];
        auto p = i + state.trace.paddings.l1;
        auto q = j - state.trace.paddings.l2;
        std::string left(state.trace.paddings.l1-1, '.');
        std::string right(state.trace.paddings.l2-1,'.');

        if (!global_visited_inside.count(make_tuple(TYPE_M2,q,p)))
            global_visited_inside[make_tuple(TYPE_M2,q,p)] = get_parentheses_inside_real_backtrace(p, q, bestM2[q][p],global_visited_inside, window_visited);
        global_visited_inside[make_tuple(TYPE_Multi,i,j)] = left + global_visited_inside[make_tuple(TYPE_M2,q,p)] + right;
        return global_visited_inside[make_tuple(TYPE_Multi,i,j)];
    }
    else if (manner == MANNER_MULTI_eq_MULTI_plus_U){
        if (global_visited_inside.count(make_tuple(TYPE_Multi,i,j)))
            return global_visited_inside[make_tuple(TYPE_Multi,i,j)];
        auto p = i + state.trace.paddings.l1;
        auto q = j - state.trace.paddings.l2;
        std::string left(state.trace.paddings.l1-1, '.');
        std::string right(state.trace.paddings.l2-1,'.');
        if (!global_visited_inside.count(make_tuple(TYPE_M2,p,q)))
            global_visited_inside[make_tuple(TYPE_M2,p,q)] = get_parentheses_inside_real_backtrace(p, q, bestM2[q][p],global_visited_inside, window_visited);
        global_visited_inside[make_tuple(TYPE_Multi,i,j)] = left + global_visited_inside[make_tuple(TYPE_M2,p,q)] + right;
        return global_visited_inside[make_tuple(TYPE_Multi,i,j)];
    }

    else if (manner == MANNER_P_eq_MULTI){
        if (global_visited_inside.count(make_tuple(TYPE_P,i,j)))
            return global_visited_inside[make_tuple(TYPE_P,i,j)];
        if (!global_visited_inside.count(make_tuple(TYPE_Multi,i,j)))
            global_visited_inside[make_tuple(TYPE_Multi,i,j)] = get_parentheses_inside_real_backtrace(i, j, bestMulti[j][i],global_visited_inside, window_visited);
        global_visited_inside[make_tuple(TYPE_P,i,j)] = '(' + global_visited_inside[make_tuple(TYPE_Multi,i,j)] + ')';
        window_fill(window_visited, i, j, seq_length, window_size);
        return global_visited_inside[make_tuple(TYPE_P,i,j)];
    }

    else if (manner == MANNER_M2_eq_M_plus_P){
        if (global_visited_inside.count(make_tuple(TYPE_M2,i,j)))
            return global_visited_inside[make_tuple(TYPE_M2,i,j)];
        auto k = state.trace.split;
        if (!global_visited_inside.count(make_tuple(TYPE_M,i,k)))
            global_visited_inside[make_tuple(TYPE_M,i,k)] = get_parentheses_inside_real_backtrace(i,k,bestM[k][i],global_visited_inside, window_visited);
        if (!global_visited_inside.count(make_tuple(TYPE_P,k+1,j)))
            global_visited_inside[make_tuple(TYPE_P,k+1,j)] = get_parentheses_inside_real_backtrace(k+1, j, bestP[j][k+1],global_visited_inside, window_visited);
        global_visited_inside[make_tuple(TYPE_M2,i,j)] = global_visited_inside[make_tuple(TYPE_M,i,k)] + global_visited_inside[make_tuple(TYPE_P,k+1,j)];
        return global_visited_inside[make_tuple(TYPE_M2,i,j)];
    }
    else if (manner == MANNER_M_eq_M2){
        if (global_visited_inside.count(make_tuple(TYPE_M,i,j)))
            return global_visited_inside[make_tuple(TYPE_M,i,j)];
        if (!global_visited_inside.count(make_tuple(TYPE_M2,i,j)))
            global_visited_inside[make_tuple(TYPE_M2,i,j)] = get_parentheses_inside_real_backtrace(i, j, bestM2[j][i],global_visited_inside, window_visited);
        global_visited_inside[make_tuple(TYPE_M,i,j)] = global_visited_inside[make_tuple(TYPE_M2,i,j)];
        return global_visited_inside[make_tuple(TYPE_M,i,j)];
    }
    else if (manner == MANNER_M_eq_M_plus_U){
        // Multi: M; M2: N; M1: m
        if (global_visited_inside.count(make_tuple(TYPE_M,i,j)))
            return global_visited_inside[make_tuple(TYPE_M,i,j)];
        if (!global_visited_inside.count(make_tuple(TYPE_M,i,j-1)))
            global_visited_inside[make_tuple(TYPE_M,i,j-1)] = get_parentheses_inside_real_backtrace(i, j-1, bestM[j-1][i],global_visited_inside, window_visited);
        global_visited_inside[make_tuple(TYPE_M,i,j)] = global_visited_inside[make_tuple(TYPE_M,i,j-1)] + '.';
        return global_visited_inside[make_tuple(TYPE_M,i,j)];
    }
    else if (manner == MANNER_M_eq_P){
        // Multi: M; M2: N; M1: m
        if (global_visited_inside.count(make_tuple(TYPE_M,i,j)))
            return global_visited_inside[make_tuple(TYPE_M,i,j)];
        if (!global_visited_inside.count(make_tuple(TYPE_P,i,j)))
            global_visited_inside[make_tuple(TYPE_P,i,j)] = get_parentheses_inside_real_backtrace(i, j, bestP[j][i],global_visited_inside, window_visited);
        global_visited_inside[make_tuple(TYPE_M,i,j)] = global_visited_inside[make_tuple(TYPE_P,i,j)];
        return global_visited_inside[make_tuple(TYPE_M,i,j)];
    }
    else if (manner == MANNER_C_eq_C_plus_U){
        assert(i == 0);
        if (global_visited_inside.count(make_tuple(TYPE_C,i,j)))
            return global_visited_inside[make_tuple(TYPE_C,i,j)];
        auto k = j - 1;
        if (!global_visited_inside.count(make_tuple(TYPE_C,i,k))) {
            if (k != -1)
                global_visited_inside[make_tuple(TYPE_C,i,k)] = get_parentheses_inside_real_backtrace(0, k, bestC[k],global_visited_inside, window_visited);
            else 
                global_visited_inside[make_tuple(TYPE_C,i,k)] = "";
        }
        global_visited_inside[make_tuple(TYPE_C,i,j)] = global_visited_inside[make_tuple(TYPE_C,i,k)] + '.';
        return global_visited_inside[make_tuple(TYPE_C,i,j)];
    }
    else if (manner == MANNER_C_eq_C_plus_P){
        assert(i == 0);
        if (global_visited_inside.count(make_tuple(TYPE_C,i,j)))
            return global_visited_inside[make_tuple(TYPE_C,i,j)];
        auto k = state.trace.split;
        if (!global_visited_inside.count(make_tuple(TYPE_C,i,k))){
            if (k != -1)
                global_visited_inside[make_tuple(TYPE_C,i,k)] = get_parentheses_inside_real_backtrace(0, k, bestC[k],global_visited_inside, window_visited);
            else
                global_visited_inside[make_tuple(TYPE_C,i,k)] = "";
        }
        if (!global_visited_inside.count(make_tuple(TYPE_P,k+1,j)))
            global_visited_inside[make_tuple(TYPE_P,k+1,j)] = get_parentheses_inside_real_backtrace(k+1, j, bestP[j][k+1],global_visited_inside, window_visited);
        global_visited_inside[make_tuple(TYPE_C,i,j)] = global_visited_inside[make_tuple(TYPE_C,i,k)] + global_visited_inside[make_tuple(TYPE_P,k+1,j)];
        return global_visited_inside[make_tuple(TYPE_C,i,j)];
    }
    else {
        printf("wrong manner inside at %d, %d: manner %d\n", i, j, state.manner); fflush(stdout);
        assert(false);
    }
}

//subopt zuker inside, sorted bestP_beta, the optimized opt zuker altorithm
pair<string, string> BeamCKYParser::get_parentheses_outside_real_backtrace(int i, int j, State& state_beta, map<tuple<BestTypes, int, int>, pair<string, string>>& global_visited_outside, map<tuple<BestTypes, int, int>, string>& global_visited_inside, set<pair<int,int>>& window_visited) {

    Manner manner = state_beta.manner;
    if (manner == MANNER_SINGLE){
        // (p, (i, j), q) from (i, j) to (p, q)
        if (global_visited_outside.count(make_tuple(TYPE_P, i, j)))
            return global_visited_outside[make_tuple(TYPE_P, i, j)];

        auto p = i - state_beta.trace.paddings.l1;  //l1 = (i-p)
        auto q = j + state_beta.trace.paddings.l2;  //l2 = (q-j)

        std::string left(state_beta.trace.paddings.l1-1, '.');
        std::string right(state_beta.trace.paddings.l2-1,'.');

        if (!global_visited_outside.count(make_tuple(TYPE_P, p, q))) 
            global_visited_outside[make_tuple(TYPE_P, p, q)] = get_parentheses_outside_real_backtrace(p, q, bestP_beta[q][p], global_visited_outside, global_visited_inside, window_visited);

        auto outsider = global_visited_outside[make_tuple(TYPE_P, p, q)];
        global_visited_outside[make_tuple(TYPE_P, i, j)] = make_pair(outsider.first + '(' + left, right + ')' + outsider.second);
        window_fill(window_visited, p, q, seq_length, window_size);
        return global_visited_outside[make_tuple(TYPE_P, i, j)];
    }

    else if (manner == MANNER_HELIX){
        if (global_visited_outside.count(make_tuple(TYPE_P, i, j)))
            return global_visited_outside[make_tuple(TYPE_P, i, j)];

        if (!global_visited_outside.count(make_tuple(TYPE_P, i-1, j+1)))
            global_visited_outside[make_tuple(TYPE_P, i-1, j+1)] = get_parentheses_outside_real_backtrace(i-1, j+1, bestP_beta[j+1][i-1], global_visited_outside, global_visited_inside, window_visited);
        auto outsider = global_visited_outside[make_tuple(TYPE_P, i-1, j+1)];

        //the two () are at position i+1, j-1, not at i, j
        global_visited_outside[make_tuple(TYPE_P, i, j)] = make_pair(outsider.first + '(', ')' + outsider.second);

        window_fill(window_visited, i-1, j+1, seq_length, window_size);
        return global_visited_outside[make_tuple(TYPE_P, i, j)];
    }

    else if (manner == MANNER_MULTI){
        // p, i, j, q
        if (global_visited_outside.count(make_tuple(TYPE_M2, i, j)))
            return global_visited_outside[make_tuple(TYPE_M2, i, j)];

        auto p = i - state_beta.trace.paddings.l1;
        auto q = j + state_beta.trace.paddings.l2;

        std::string left(state_beta.trace.paddings.l1-1, '.');
        std::string right(state_beta.trace.paddings.l2-1,'.');

        if (!global_visited_outside.count(make_tuple(TYPE_Multi, p, q)))
            global_visited_outside[make_tuple(TYPE_Multi, p, q)] = get_parentheses_outside_real_backtrace(p, q, bestMulti_beta[q][p], global_visited_outside, global_visited_inside, window_visited);
        auto outsider = global_visited_outside[make_tuple(TYPE_Multi, p, q)];
        global_visited_outside[make_tuple(TYPE_M2,i,j)] = make_pair(outsider.first + left, right + outsider.second);
        return global_visited_outside[make_tuple(TYPE_M2,i,j)];
    }

    else if (manner == MANNER_MULTI_eq_MULTI_plus_U){
        // p , i , j , q
        if (global_visited_outside.count(make_tuple(TYPE_Multi, i, j)))
            return global_visited_outside[make_tuple(TYPE_Multi, i, j)];
        int j_next = j + state_beta.trace.split;
        std::string right(state_beta.trace.split ,'.'); //this padding is not the same as others, we need to add '.' to replace the previous ')'.

        if (!global_visited_outside.count(make_tuple(TYPE_Multi, i, j_next)))
            global_visited_outside[make_tuple(TYPE_Multi, i, j_next)] = get_parentheses_outside_real_backtrace(i, j_next, bestMulti_beta[j_next][i], global_visited_outside, global_visited_inside, window_visited);

        auto outsider = global_visited_outside[make_tuple(TYPE_Multi, i, j_next)];
        global_visited_outside[make_tuple(TYPE_Multi, i, j)] = make_pair(outsider.first, right + outsider.second);
        // can only add right to rightside. outsider.first, (?????? inside ?????????) , right + outsider.second
        return global_visited_outside[make_tuple(TYPE_Multi, i, j)];
    }

    else if (manner == MANNER_P_eq_MULTI){
        if (global_visited_outside.count(make_tuple(TYPE_Multi, i, j)))
            return global_visited_outside[make_tuple(TYPE_Multi, i, j)];

        if (!global_visited_outside.count(make_tuple(TYPE_P, i, j)))
            global_visited_outside[make_tuple(TYPE_P, i, j)] = get_parentheses_outside_real_backtrace(i, j, bestP_beta[j][i], global_visited_outside, global_visited_inside, window_visited);

        auto outsider = global_visited_outside[make_tuple(TYPE_P, i, j)];
        global_visited_outside[make_tuple(TYPE_Multi, i, j)] = make_pair(outsider.first + '(', ')' + outsider.second);
        window_fill(window_visited, i, j, seq_length, window_size);

        return global_visited_outside[make_tuple(TYPE_Multi, i, j)];
    }

    else if (manner == MANNER_M2_eq_M_plus_P){
        // M2 = M1 + P
        // m, k=i-1,          i, j 

        int mm, kk, ii, jj;

        // M
        // M2 = M1 + P, it is from M; we need M2 outside and P inside
        if (state_beta.trace.split > j){
            mm = i;
            kk = j;
            ii = j+1;
            jj = state_beta.trace.split;

            if (global_visited_outside.count(make_tuple(TYPE_M,mm,kk)))
                return global_visited_outside[make_tuple(TYPE_M,mm,kk)];

            if (!global_visited_inside.count(make_tuple(TYPE_P,ii,jj)))
                global_visited_inside[make_tuple(TYPE_P,ii,jj)] = get_parentheses_inside_real_backtrace(ii, jj, bestP[jj][ii], global_visited_inside, window_visited);

            auto P_inside = global_visited_inside[make_tuple(TYPE_P,ii,jj)];

            if (!global_visited_outside.count(make_tuple(TYPE_M2,mm,jj)))
                global_visited_outside[make_tuple(TYPE_M2,mm,jj)] = get_parentheses_outside_real_backtrace(mm, jj, bestM2_beta[jj][mm], global_visited_outside, global_visited_inside, window_visited);

            auto outsider = global_visited_outside[make_tuple(TYPE_M2,mm,jj)];
            global_visited_outside[make_tuple(TYPE_M,mm,kk)] = make_pair(outsider.first, P_inside + outsider.second);
            return global_visited_outside[make_tuple(TYPE_M,mm,kk)];
        }

        // M2 = M + P, it is from P; we need M2 outside and M1 inside
        else{
            assert(state_beta.trace.split < i);
            mm = state_beta.trace.split;
            kk = i-1;
            ii = i;
            jj = j;

            if (global_visited_outside.count(make_tuple(TYPE_P,ii,jj)))
                return global_visited_outside[make_tuple(TYPE_P,ii,jj)];

            if (!global_visited_inside.count(make_tuple(TYPE_M,mm,kk)))
                global_visited_inside[make_tuple(TYPE_M,mm,kk)] = get_parentheses_inside_real_backtrace(mm, kk, bestM[kk][mm],global_visited_inside, window_visited);
            auto M1_inside = global_visited_inside[make_tuple(TYPE_M,mm,kk)];

            if (!global_visited_outside.count(make_tuple(TYPE_M2,mm,jj)))
                global_visited_outside[make_tuple(TYPE_M2,mm,jj)] = get_parentheses_outside_real_backtrace(mm, jj, bestM2_beta[jj][mm], global_visited_outside, global_visited_inside, window_visited);

            auto outsider = global_visited_outside[make_tuple(TYPE_M2,mm,jj)];
            global_visited_outside[make_tuple(TYPE_P,ii,jj)] = make_pair(outsider.first +  M1_inside, outsider.second);
            return global_visited_outside[make_tuple(TYPE_P,ii,jj)];
        }
    }

    else if (manner == MANNER_M_eq_M2){
        // M1 = M2
        // Multi: M; M2: N; M1: m

        if (global_visited_outside.count(make_tuple(TYPE_M2, i, j))) 
            return global_visited_outside[make_tuple(TYPE_M2, i, j)];

        if (!global_visited_outside.count(make_tuple(TYPE_M, i, j))) 
            global_visited_outside[make_tuple(TYPE_M, i, j)] = get_parentheses_outside_real_backtrace(i, j, bestM_beta[j][i], global_visited_outside, global_visited_inside, window_visited);

        global_visited_outside[make_tuple(TYPE_M2, i, j)] = global_visited_outside[make_tuple(TYPE_M, i, j)];
        return global_visited_outside[make_tuple(TYPE_M2, i, j)];
    }

    else if (manner == MANNER_M_eq_M_plus_U){
        // Multi: M; M2: N; M1(M): m

        if (global_visited_outside.count(make_tuple(TYPE_M,i,j)))
            return global_visited_outside[make_tuple(TYPE_M,i,j)];

        if (!global_visited_outside.count(make_tuple(TYPE_M,i,j+1)))
            global_visited_outside[make_tuple(TYPE_M,i,j+1)] = get_parentheses_outside_real_backtrace(i, j+1, bestM_beta[j+1][i], global_visited_outside, global_visited_inside, window_visited);

        auto outsider = global_visited_outside[make_tuple(TYPE_M,i,j+1)];
        global_visited_outside[make_tuple(TYPE_M,i,j)] = make_pair(outsider.first, '.' + outsider.second);
        return global_visited_outside[make_tuple(TYPE_M,i,j)];
    }

    else if (manner == MANNER_M_eq_P){
        if (global_visited_outside.count(make_tuple(TYPE_P,i,j)))
            return global_visited_outside[make_tuple(TYPE_P,i,j)];

        if (!global_visited_outside.count(make_tuple(TYPE_M,i,j))) 
            global_visited_outside[make_tuple(TYPE_M,i,j)] = get_parentheses_outside_real_backtrace(i, j, bestM_beta[j][i], global_visited_outside, global_visited_inside, window_visited);

        global_visited_outside[make_tuple(TYPE_P, i, j)] = global_visited_outside[make_tuple(TYPE_M, i, j)];
        return global_visited_outside[make_tuple(TYPE_P,i, j)];
    }

    else if (manner == MANNER_C_eq_C_plus_U){

        assert(i == 0);
        assert(j+1 != seq_length); //since it is even outside the base case, where base case is j + 1 == seq_length - 1

        if (global_visited_outside.count(make_tuple(TYPE_C, i, j)))
            return global_visited_outside[make_tuple(TYPE_C, i, j)];

        if (!global_visited_outside.count(make_tuple(TYPE_C, i, j+1))){
            if (j + 1 < seq_length - 1)
                global_visited_outside[make_tuple(TYPE_C, i, j+1)] = get_parentheses_outside_real_backtrace(0, j+1, bestC_beta[j+1], global_visited_outside, global_visited_inside, window_visited);
            else
                global_visited_outside[make_tuple(TYPE_C, i, j+1)] = make_pair("","");
        }
        auto outsider = global_visited_outside[make_tuple(TYPE_C, i, j+1)];
        global_visited_outside[make_tuple(TYPE_C, i, j)] = make_pair(outsider.first, '.' + outsider.second);
        return global_visited_outside[make_tuple(TYPE_C, i, j)];
    }

    else if (manner == MANNER_C_eq_C_plus_P){

        int kk, ii, jj;

        //C = c + P, it is from P, we need C outside, and c inside
        if (state_beta.trace.split == -1){
            kk = i-1;
            ii = i;
            jj = j;

            if (global_visited_outside.count(make_tuple(TYPE_P,ii,jj))) 
                return global_visited_outside[make_tuple(TYPE_P,ii,jj)];

            if (kk != -1){
                if (!global_visited_inside.count(make_tuple(TYPE_C,0,kk))) 
                    global_visited_inside[make_tuple(TYPE_C,0,kk)] = get_parentheses_inside_real_backtrace(0, kk, bestC[kk], global_visited_inside, window_visited);

                auto inside_C = global_visited_inside[make_tuple(TYPE_C,0,kk)];
                if (!global_visited_outside.count(make_tuple(TYPE_C,0,jj))) {
                    if (jj < seq_length - 1)
                        global_visited_outside[make_tuple(TYPE_C, 0, jj)] = get_parentheses_outside_real_backtrace(0, jj, bestC_beta[jj], global_visited_outside, global_visited_inside, window_visited);
                    else
                        global_visited_outside[make_tuple(TYPE_C, 0, jj)] = make_pair("","");
                }

                auto outsider = global_visited_outside[make_tuple(TYPE_C,0,jj)];
                global_visited_outside[make_tuple(TYPE_P, ii, jj)] = make_pair(outsider.first + inside_C, outsider.second);
                return global_visited_outside[make_tuple(TYPE_P, ii, jj)];
            }

            else{
                assert(ii == 0);

                if (!global_visited_outside.count(make_tuple(TYPE_C,0,jj))) {
                    if (jj < seq_length - 1)
                        global_visited_outside[make_tuple(TYPE_C,0,jj)] = get_parentheses_outside_real_backtrace(0, jj, bestC_beta[jj], global_visited_outside, global_visited_inside, window_visited);
                    else
                        global_visited_outside[make_tuple(TYPE_C,0,jj)] = make_pair("","");
                }
                global_visited_outside[make_tuple(TYPE_P, ii, jj)] = global_visited_outside[make_tuple(TYPE_C,0,jj)];
                return global_visited_outside[make_tuple(TYPE_P, ii, jj)];
            }
        }

        // C = c + P, it is from small c, we need C outside and P inside
        else {
            assert(i == 0);
            assert(state_beta.trace.split != -1);

            kk = j;
            ii = j+1;
            jj = state_beta.trace.split;

            if (global_visited_outside.count(make_tuple(TYPE_C,0,kk)))
                return global_visited_outside[make_tuple(TYPE_C,0,kk)];

            if (!global_visited_inside.count(make_tuple(TYPE_P, ii, jj))) 
                global_visited_inside[make_tuple(TYPE_P, ii, jj)] = get_parentheses_inside_real_backtrace(ii, jj, bestP[jj][ii], global_visited_inside, window_visited);

            auto inside_P = global_visited_inside[make_tuple(TYPE_P, ii, jj)];

            if (!global_visited_outside.count(make_tuple(TYPE_C,0,jj))){
                if (jj < seq_length - 1)
                    global_visited_outside[make_tuple(TYPE_C,0,jj)] = get_parentheses_outside_real_backtrace(0, jj, bestC_beta[jj], global_visited_outside, global_visited_inside, window_visited);
                else
                    global_visited_outside[make_tuple(TYPE_C,0,jj)] = make_pair("","");
            }

            auto outsider = global_visited_outside[make_tuple(TYPE_C,0,jj)];
            global_visited_outside[make_tuple(TYPE_C,0,kk)] = make_pair(outsider.first, inside_P + outsider.second);
            return global_visited_outside[make_tuple(TYPE_C,0,kk)];
        }
    }

    else {
        printf("wrong manner outside at %d, %d: manner %d\n", i, j, state_beta.manner); fflush(stdout);
        assert(false);
    }
}


void BeamCKYParser::prepare(unsigned len) {
    seq_length = len;

    bestH.clear();
    bestH.resize(seq_length);
    bestP.clear();
    bestP.resize(seq_length);
    bestM2.clear();
    bestM2.resize(seq_length);
    bestM.clear();
    bestM.resize(seq_length);
    bestC.clear();
    bestC.resize(seq_length);
    bestMulti.clear();
    bestMulti.resize(seq_length);

#ifdef is_cube_pruning
        sorted_bestM.clear();
        sorted_bestM.resize(seq_length);
#endif

    nucs.clear();
    nucs.resize(seq_length);

    scores.reserve(seq_length);

    if (use_constraints){
        allow_unpaired_position.clear();
        allow_unpaired_position.resize(seq_length);

        allow_unpaired_range.clear();
        allow_unpaired_range.resize(seq_length);
    }
}

// lisiz, constraints
bool BeamCKYParser::allow_paired(int i, int j, vector<int>* cons, char nuci, char nucj) {
    return ((*cons)[i] == -1 || (*cons)[i] == j) && ((*cons)[j] == -1 || (*cons)[j] == i) && _allowed_pairs[nuci][nucj];
}

BeamCKYParser::DecoderResult BeamCKYParser::parse(string& seq, vector<int>* cons) {

    struct timeval parse_starttime, parse_endtime;

    // number of states
    unsigned long nos_H = 0, nos_P = 0, nos_M2 = 0,
            nos_M = 0, nos_C = 0, nos_Multi = 0;
    gettimeofday(&parse_starttime, NULL);
    prepare(static_cast<unsigned>(seq.length()));

    for (int i = 0; i < seq_length; ++i)
        nucs[i] = GET_ACGU_NUM(seq[i]);

    // lisiz, constraints
    if (use_constraints) {
        for (int i=0; i<seq_length; i++){
            int cons_idx = (*cons)[i];
            allow_unpaired_position[i] = cons_idx == -1 || cons_idx == -2;
            if (cons_idx > -1){
                if (!_allowed_pairs[nucs[i]][nucs[cons_idx]]){
                    printf("Constrains on non-classical base pairs (non AU, CG, GU pairs)\n");
                    exit(1);
                }
            }
        }
        int firstpair = seq_length;
        for (int i=seq_length-1; i>-1; i--){
            allow_unpaired_range[i] = firstpair;
            if ((*cons)[i] >= 0)
                firstpair = i;
        }
    }

    vector<int> next_pair[NOTON];
    {
        if (use_constraints){
            for (int nuci = 0; nuci < NOTON; ++nuci) {
                next_pair[nuci].resize(seq_length, -1);
                int next = -1;
                for (int j = seq_length-1; j >=0; --j) {
                    next_pair[nuci][j] = next;
                    if ((*cons)[j] > -2 && _allowed_pairs[nuci][nucs[j]]) next = j;
                }
            }
        } else {
            for (int nuci = 0; nuci < NOTON; ++nuci) {
                next_pair[nuci].resize(seq_length, -1);
                int next = -1;
                for (int j = seq_length-1; j >=0; --j) {
                    next_pair[nuci][j] = next;
                    if (_allowed_pairs[nuci][nucs[j]]) next = j;
                }
            }
        }
    }

#ifdef SPECIAL_HP
#ifdef lv
    v_init_tetra_hex_tri(seq, seq_length, if_tetraloops, if_hexaloops, if_triloops);
#else
    if (is_verbose)
        v_init_tetra_hex_tri(seq, seq_length, if_tetraloops, if_hexaloops, if_triloops);
#endif
#endif

    // start CKY decoding
#ifdef lv
        if(seq_length > 0) bestC[0].set(- v_score_external_unpaired(0, 0), MANNER_C_eq_C_plus_U);
        if(seq_length > 1) bestC[1].set(- v_score_external_unpaired(0, 1), MANNER_C_eq_C_plus_U);
#else
        if(seq_length > 0) bestC[0].set(score_external_unpaired(0, 0), MANNER_C_eq_C_plus_U);
        if(seq_length > 1) bestC[1].set(score_external_unpaired(0, 1), MANNER_C_eq_C_plus_U);
#endif
    ++nos_C;

    // from left to right
    for(int j = 0; j < seq_length; ++j) {
        int nucj = nucs[j];
        int nucj1 = (j+1) < seq_length ? nucs[j+1] : -1;

        unordered_map<int, State>& beamstepH = bestH[j];
        unordered_map<int, State>& beamstepMulti = bestMulti[j];
        unordered_map<int, State>& beamstepP = bestP[j];
        unordered_map<int, State>& beamstepM2 = bestM2[j];
        unordered_map<int, State>& beamstepM = bestM[j];
        State& beamstepC = bestC[j];

        // beam of H
        {
            if (beam > 0 && beamstepH.size() > beam) beam_prune(beamstepH);

            {
                int jnext = next_pair[nucj][j];
                if (no_sharp_turn) while (jnext - j < 4 && jnext != -1) jnext = next_pair[nucj][jnext];

                // lisiz, constriants
                if (use_constraints){
                    if (!allow_unpaired_position[j]){
                        jnext = (*cons)[j] > j ? (*cons)[j] : -1; // lisiz: j must be left bracket, jump to the constrainted pair (j, j') directly
                    }
                    if (jnext != -1){
                        int nucjnext = nucs[jnext];
                        if (jnext > allow_unpaired_range[j] || !allow_paired(j, jnext, cons, nucj, nucjnext))  // lisiz: avoid cross constrainted brackets or unallowed pairs
                            jnext = -1;
                    }
                }

                if (jnext != -1) {
                    int nucjnext = nucs[jnext];
                    int nucjnext_1 = (jnext - 1) > -1 ? nucs[jnext - 1] : -1;

                    value_type newscore;

#ifdef lv
                        int tetra_hex_tri = -1;
#ifdef SPECIAL_HP
                        if (jnext-j-1 == 4) // 6:tetra
                            tetra_hex_tri = if_tetraloops[j];
                        else if (jnext-j-1 == 6) // 8:hexa
                            tetra_hex_tri = if_hexaloops[j];
                        else if (jnext-j-1 == 3) // 5:tri
                            tetra_hex_tri = if_triloops[j];
#endif
                        newscore = - v_score_hairpin(j, jnext, nucj, nucj1, nucjnext_1, nucjnext, tetra_hex_tri);
#else
                        newscore = score_hairpin(j, jnext, nucj, nucj1, nucjnext_1, nucjnext);
#endif
                    // this candidate must be the best one at [j, jnext]
                    // so no need to check the score
                    update_if_better(bestH[jnext][j], newscore, MANNER_H);
                    ++ nos_H;
                }
            }

            {
                // for every state h in H[j]
                //   1. extend h(i, j) to h(i, jnext)
                //   2. generate p(i, j)
#ifdef lv
                sort_keys(beamstepH, keys);
                for (auto &item : keys) {
#else
                for (auto &item : beamstepH) {
#endif
                    int i = item.first;
                    // printf("%d\n", i);
                    State &state = item.second;
                    int nuci = nucs[i];
                    int jnext = next_pair[nuci][j];

                    // 2. generate p(i, j)
                    // lisiz, change the order because of the constriants
                    {
                        update_if_better(beamstepP[i], state.score, MANNER_HAIRPIN);
                        ++ nos_P;
                    }

                    // lisiz, constraints
                    if (jnext != -1 && use_constraints){
                        int nucjnext = nucs[jnext];
                        if (jnext > allow_unpaired_range[i] || !allow_paired(i, jnext, cons, nuci, nucjnext))
                            continue;
                    }

                    if (jnext != -1) {
                        int nuci1 = (i + 1) < seq_length ? nucs[i + 1] : -1;
                        int nucjnext = nucs[jnext];
                        int nucjnext_1 = (jnext - 1) > -1 ? nucs[jnext - 1] : -1;

                        // 1. extend h(i, j) to h(i, jnext)
                        value_type newscore;

#ifdef lv
                            int tetra_hex_tri = -1;
#ifdef SPECIAL_HP
                            if (jnext-i-1 == 4) // 6:tetra
                                tetra_hex_tri = if_tetraloops[i];
                            else if (jnext-i-1 == 6) // 8:hexa
                                tetra_hex_tri = if_hexaloops[i];
                            else if (jnext-i-1 == 3) // 5:tri
                                tetra_hex_tri = if_triloops[i];
#endif
                            newscore = - v_score_hairpin(i, jnext, nuci, nuci1, nucjnext_1, nucjnext, tetra_hex_tri);
#else
                            newscore = score_hairpin(i, jnext, nuci, nuci1, nucjnext_1, nucjnext);
#endif
                        // this candidate must be the best one at [i, jnext]
                        // so no need to check the score
                        update_if_better(bestH[jnext][i], newscore, MANNER_H);
                        ++nos_H;
                    }
                }
            }
        }
        if (j == 0) continue;

        // beam of Multi
        {
            if (beam > 0 && beamstepMulti.size() > beam) beam_prune(beamstepMulti);

            // for every state in Multi[j]
            //   1. extend (i, j) to (i, jnext)
            //   2. generate P (i, j)
#ifdef lv
            sort_keys(beamstepMulti, keys);
            for (auto &item : keys) {
#else
            for(auto& item : beamstepMulti) {
#endif
                int i = item.first;
                State& state = item.second;
                int nuci = nucs[i];
                int nuci1 = nucs[i+1];
                int jnext = next_pair[nuci][j];

                // 2. generate P (i, j)
                // lisiz, change the order because of the constraits
                {
                    value_type newscore;
#ifdef lv
                        newscore = state.score - v_score_multi(i, j, nuci, nuci1, nucs[j-1], nucj, seq_length, dangle_model);
#else
                        newscore = state.score + score_multi(i, j, nuci, nuci1, nucs[j-1], nucj, seq_length);
#endif
                    update_if_better(beamstepP[i], newscore, MANNER_P_eq_MULTI);
                    ++ nos_P;
                }

                // lisiz cnstriants
                if (jnext != -1 && use_constraints){
                    int nucjnext = nucs[jnext];
                    if (jnext > allow_unpaired_range[j] || !allow_paired(i, jnext, cons, nuci, nucjnext))
                        continue;
                }

                // 1. extend (i, j) to (i, jnext)
                {
                    char new_l1 = state.trace.paddings.l1;
                    int new_l2 = state.trace.paddings.l2 + jnext - j;
                    // if (jnext != -1 && new_l1 + new_l2 <= SINGLE_MAX_LEN) {
                    if (jnext != -1) {
                        // 1. extend (i, j) to (i, jnext)
                        value_type newscore;
#ifdef lv
                            newscore = state.score - v_score_multi_unpaired(j, jnext - 1);
#else
                            newscore = state.score + score_multi_unpaired(j, jnext - 1);
#endif
                        // this candidate must be the best one at [i, jnext]
                        // so no need to check the score
                        update_if_better(bestMulti[jnext][i], newscore, MANNER_MULTI_eq_MULTI_plus_U,
                                         new_l1,
                                         new_l2
                        );
                        ++nos_Multi;
                    }
                }
            }
        }

        // beam of P
        {
            if (beam > 0 && beamstepP.size() > beam) beam_prune(beamstepP);

            // for every state in P[j]
            //   1. generate new helix/bulge
            //   2. M = P
            //   3. M2 = M + P
            //   4. C = C + P
#ifdef is_cube_pruning
            bool use_cube_pruning = beam > MIN_CUBE_PRUNING_SIZE
                                    && beamstepP.size() > MIN_CUBE_PRUNING_SIZE;
#else
            bool use_cube_pruning = false;
#endif               

#ifdef lv
            sort_keys(beamstepP, keys);
            for (auto &item : keys) {
#else
            for(auto& item : beamstepP) {
#endif
                int i = item.first;
                State& state = item.second;
                int nuci = nucs[i];
                int nuci_1 = (i-1>-1) ? nucs[i-1] : -1;

                // 2. M = P // lhuang: check j (this P is not the last non-closing P in multi)
                if(i > 0 && j < seq_length-6){
                    value_type newscore;
#ifdef lv
                        newscore = - v_score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length, dangle_model) + state.score;
#else
                        newscore = score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length) + state.score;
#endif
                    update_if_better(beamstepM[i], newscore, MANNER_M_eq_P);
                    ++ nos_M;
                }
                //printf(" M = P at %d\n", j); fflush(stdout);

                // 3. M2 = M + P // lhuang: check j < n-1
                if(!use_cube_pruning && j < seq_length - 1) {
                    int k = i - 1;
                    if ( k > 0 && !bestM[k].empty()) {
                        value_type M1_score;
#ifdef lv
                            M1_score = - v_score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length, dangle_model) + state.score;
#else
                            M1_score = score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length) + state.score;
#endif
                        // candidate list
                        auto bestM2_iter = beamstepM2.find(i);
#ifndef is_candidate_list
                        for (auto &m : bestM[k]) {
                                int newi = m.first;
                                // eq. to first convert P to M1, then M2/M = M + M1
                                value_type newscore = M1_score + m.second.score;
                                update_if_better(beamstepM2[newi], newscore, MANNER_M2_eq_M_plus_P, k);
                                ++nos_M2;
                            }
#else
                        if (bestM2_iter==beamstepM2.end() || M1_score > bestM2_iter->second.score) {
                            for (auto &m : bestM[k]) {
                                int newi = m.first;
                                // eq. to first convert P to M1, then M2/M = M + M1
                                value_type newscore = M1_score + m.second.score;
                                update_if_better(beamstepM2[newi], newscore, MANNER_M2_eq_M_plus_P, k);
                                ++nos_M2;
                            }
                        }
#endif
                    }
                }
                //printf(" M/M2 = M + P at %d\n", j); fflush(stdout);

                // 4. C = C + P
                {
                    int k = i - 1;
                    if (k >= 0) {
                      State& prefix_C = bestC[k];
                      if (prefix_C.manner != MANNER_NONE) {
                        int nuck = nuci_1;
                        int nuck1 = nuci;
                        value_type newscore;
#ifdef lv
                            newscore = - v_score_external_paired(k+1, j, nuck, nuck1, nucj, nucj1, seq_length, dangle_model) +
                                prefix_C.score + state.score;
#else
                            newscore = score_external_paired(k+1, j, nuck, nuck1, nucj, nucj1, seq_length) +
                                prefix_C.score + state.score;
#endif
                        update_if_better(beamstepC, newscore, MANNER_C_eq_C_plus_P, k);
                        ++ nos_C;
                      }
                    } else {
                        value_type newscore;
#ifdef lv
                            newscore = - v_score_external_paired(0, j, -1, nucs[0], nucj, nucj1, seq_length, dangle_model) + state.score;
#else
                            newscore = score_external_paired(0, j, -1, nucs[0], nucj, nucj1, seq_length) + state.score;
#endif
                        update_if_better(beamstepC, newscore, MANNER_C_eq_C_plus_P, -1);
                        ++ nos_C;
                    }
                }
                //printf(" C = C + P at %d\n", j); fflush(stdout);

                // 1. generate new helix / single_branch
                // new state is of shape p..i..j..q
                if (i >0 && j<seq_length-1) {
                    value_type precomputed;
#ifdef lv
                        precomputed = 0;
#else
                        precomputed = score_junction_B(j, i, nucj, nucj1, nuci_1, nuci);
#endif
                    for (int p = i - 1; p >= std::max(i - SINGLE_MAX_LEN, 0); --p) {
                        int nucp = nucs[p];
                        int nucp1 = nucs[p + 1]; // hzhang: move here
                        int q = next_pair[nucp][j];

                        // lisiz constraints
                        if (use_constraints){
                            if (p < i-1 && !allow_unpaired_position[p+1]) // lisiz: if p+1 must be paired, break
                                break;
                            if (!allow_unpaired_position[p]){             // lisiz: if p must be paired, p must be left bracket
                                q = (*cons)[p];
                                if (q < p) break;
                            }
                        }

                        while (q != -1 && ((i - p) + (q - j) - 2 <= SINGLE_MAX_LEN)) {
                            int nucq = nucs[q];

                            // lisiz constraints
                            if (use_constraints){
                                if (q>j+1 && q > allow_unpaired_range[j])  // lisiz: if q-1 must be paired, break
                                    break;
                                if (!allow_paired(p, q, cons, nucp, nucq)) // lisiz: if p q are )(, break
                                    break;
                            }

                            int nucq_1 = nucs[q - 1];
                            if (p == i - 1 && q == j + 1) {
                                // helix
                                value_type newscore;
#ifdef lv
                                    newscore = -v_score_single(p,q,i,j, nucp, nucp1, nucq_1, nucq, nuci_1, nuci, nucj, nucj1)
                                        + state.score;
                                    // SHAPE for Vienna only
                                    if (use_shape)
                                        newscore += -(pseudo_energy_stack[p] + pseudo_energy_stack[i] + pseudo_energy_stack[j] + pseudo_energy_stack[q]);
#else
                                    newscore = score_helix(nucp, nucp1, nucq_1, nucq) + state.score;
#endif
                                update_if_better(bestP[q][p], newscore, MANNER_HELIX);
                                ++nos_P;
                            } else {
                                // single branch
                                value_type newscore;
#ifdef lv
                                    newscore = - v_score_single(p,q,i,j, nucp, nucp1, nucq_1, nucq, nuci_1, nuci, nucj, nucj1)
                                        + state.score;
#else
                                    newscore = score_junction_B(p, q, nucp, nucp1, nucq_1, nucq) +
                                        precomputed +
                                        score_single_without_junctionB(p, q, i, j, nuci_1, nuci, nucj, nucj1) +
                                        state.score;
#endif
                                update_if_better(bestP[q][p], newscore, MANNER_SINGLE,
                                                 static_cast<char>(i - p),
                                                 q - j);
                                ++nos_P;
                            }
                            q = next_pair[nucp][q];
                        }
                    }
                }
            }

	    // lhuang: check j < n-1
            if (use_cube_pruning && j < seq_length - 1) {
                // 3. M2 = M + P with cube pruning
                vector<int> valid_Ps;
                vector<value_type> M1_scores;
#ifdef lv
            sort_keys(beamstepP, keys);
            for (auto &item : keys) {
#else
            for(auto& item : beamstepP) {
#endif
                    int i = item.first;
                    State &state = item.second;
                    int nuci = nucs[i];
                    int nuci_1 = (i - 1 > -1) ? nucs[i - 1] : -1;
                    int k = i - 1;

                    // group candidate Ps
                    if (k > 0 && !bestM[k].empty()) {
                        assert(bestM[k].size() == sorted_bestM[k].size());
                        value_type M1_score;
#ifdef lv
                            M1_score = - v_score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length, dangle_model) + state.score;
#else
                            M1_score = score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length) + state.score;
#endif
                        auto bestM2_iter = beamstepM2.find(i);
#ifndef is_candidate_list
                        valid_Ps.push_back(i);
                        M1_scores.push_back(M1_score);
#else
                        if (bestM2_iter == beamstepM2.end() || M1_score > bestM2_iter->second.score) {
                            valid_Ps.push_back(i);
                            M1_scores.push_back(M1_score);
                        }
#endif
                    }
                }

                // build max heap
                // heap is of form (heuristic score, (index of i in valid_Ps, index of M in bestM[i-1]))
                vector<pair<value_type, pair<int, int>>> heap;
                for (int p = 0; p < valid_Ps.size(); ++p) {
                    int i = valid_Ps[p];
                    int k = i - 1;
                    heap.push_back(make_pair(M1_scores[p] + sorted_bestM[k][0].first,
                                             make_pair(p, 0)
                    ));
                    push_heap(heap.begin(), heap.end());
                }

                // start cube pruning
                // stop after beam size M2 states being filled
                int filled = 0;
                // exit when filled >= beam and current score < prev score
                value_type prev_score = VALUE_MIN;
                value_type current_score = VALUE_MIN;
                while ((filled < beam || current_score == prev_score) && !heap.empty()) {
                    auto &top = heap.front();
                    prev_score = current_score;
                    current_score = top.first;
                    int index_P = top.second.first;
                    int index_M = top.second.second;
                    int i = valid_Ps[top.second.first];
                    int k = i - 1;
                    int newi = sorted_bestM[k][index_M].second;
                    value_type newscore = M1_scores[index_P] + bestM[k][newi].score;
                    pop_heap(heap.begin(), heap.end());
                    heap.pop_back();

                    if (beamstepM2[newi].manner == MANNER_NONE) {
                        ++filled;
                        update_if_better(beamstepM2[newi], newscore, MANNER_M2_eq_M_plus_P, k);
                        ++nos_M2;
                    } else {
                        assert(beamstepM2[newi].score > newscore - 1e-8);
                    }

                    ++index_M;
                    while (index_M < sorted_bestM[k].size()) {
                        // candidate_score is a heuristic score
                        value_type candidate_score = M1_scores[index_P] + sorted_bestM[k][index_M].first;
                        int candidate_newi = sorted_bestM[k][index_M].second;
                        if (beamstepM2.find(candidate_newi) == beamstepM2.end()) {
                            heap.push_back(make_pair(candidate_score,
                                                     make_pair(index_P, index_M)));
                            push_heap(heap.begin(), heap.end());
                            break;
                        } else {
                            // based on the property of cube pruning, the new score must be worse
                            // than the state already inserted
                            // so we keep iterate through the candidate list to find the next
                            // candidate
                            ++index_M;
                            assert(beamstepM2[candidate_newi].score >
                                   M1_scores[index_P] + bestM[k][candidate_newi].score - 1e-8);
                        }
                    }
                }
            }
        }

        // beam of M2
        {
            if (beam > 0 && beamstepM2.size() > beam) beam_prune(beamstepM2);

            // for every state in M2[j]
            //   1. multi-loop  (by extending M2 on the left)
            //   2. M = M2
#ifdef lv
            sort_keys(beamstepM2, keys);
            for (auto &item : keys) {
#else
            for(auto& item : beamstepM2) {
#endif
                int i = item.first;
                State& state = item.second;

                // 2. M = M2
                {
                    update_if_better(beamstepM[i], state.score, MANNER_M_eq_M2);
                    ++ nos_M;
                }

                // 1. multi-loop
                {
                    for (int p = i-1; p >= std::max(i - SINGLE_MAX_LEN, 0); --p) {
                        int nucp = nucs[p];
                        int q = next_pair[nucp][j];

                        if (use_constraints){
                            if (p < i - 1 && !allow_unpaired_position[p+1])
                                break;
                            if (!allow_unpaired_position[p]){
                                q = (*cons)[p];
                                if (q < p) break;
                            }
                            if (q > j+1 && q > allow_unpaired_range[j])
                                continue;
                            int nucq = nucs[q];
                            if (!allow_paired(p, q, cons, nucp, nucq))
                                continue;
                        }

                        if (q != -1 && ((i - p - 1) <= SINGLE_MAX_LEN)) {
                            // the current shape is p..i M2 j ..q

                            value_type newscore;
#ifdef lv
                                newscore = - v_score_multi_unpaired(p+1, i-1) -
                                    v_score_multi_unpaired(j+1, q-1) + state.score;
#else
                                newscore = score_multi_unpaired(p+1, i-1) +
                                    score_multi_unpaired(j+1, q-1) + state.score;
#endif
                            update_if_better(bestMulti[q][p], newscore, MANNER_MULTI,
                                             static_cast<char>(i - p),
                                             q - j);
                            ++ nos_Multi;
                        }
                    }
                }
            }
        }

        // beam of M
        {
            value_type threshold = VALUE_MIN;
            if (beam > 0 && beamstepM.size() > beam) threshold = beam_prune(beamstepM);

#ifdef is_cube_pruning
                sortM(threshold, beamstepM, sorted_bestM[j]);
            // }
#endif

            // for every state in M[j]
            //   1. M = M + unpaired
#ifdef lv
            sort_keys(beamstepM, keys);
            for (auto &item : keys) {
#else
            for(auto& item : beamstepM) {
#endif
                int i = item.first;
                State& state = item.second;
                if (j < seq_length-1) {
                    if (use_constraints && !allow_unpaired_position[j+1]) // if j+1 must be paired
                        continue;
                    value_type newscore;
#ifdef lv
                        newscore = - v_score_multi_unpaired(j + 1, j + 1) + state.score;
#else
                        newscore = score_multi_unpaired(j + 1, j + 1) + state.score;
#endif
                    update_if_better(bestM[j+1][i], newscore, MANNER_M_eq_M_plus_U);
                    ++ nos_M;
                }
            }
        }

        // beam of C
        {
            // C = C + U
            if (j < seq_length-1) {
                if (use_constraints && !allow_unpaired_position[j+1])
                        continue;
                value_type newscore;
#ifdef lv
                    newscore = -v_score_external_unpaired(j+1, j+1) + beamstepC.score;
#else
                    newscore = score_external_unpaired(j+1, j+1) + beamstepC.score;
#endif
                update_if_better(bestC[j+1], newscore, MANNER_C_eq_C_plus_U);
                ++ nos_C;
            }
        }

    }  // end of for-loo j

    State& viterbi = bestC[seq_length-1];

    char result[seq_length+1];
    get_parentheses(result, seq);

    gettimeofday(&parse_endtime, NULL);
    double parse_elapsed_time = parse_endtime.tv_sec - parse_starttime.tv_sec + (parse_endtime.tv_usec-parse_starttime.tv_usec)/1000000.0;

    unsigned long nos_tot = nos_H + nos_P + nos_M2 + nos_Multi + nos_M + nos_C;
    if (is_verbose) {
        printf("Parse Time: %f len: %d score %f #states %lu H %lu P %lu M2 %lu Multi %lu M %lu C %lu\n",
               parse_elapsed_time, seq_length, double(viterbi.score), nos_tot,
               nos_H, nos_P, nos_M2, nos_Multi, nos_M, nos_C);
    }

    fflush(stdout);

    if (zuker){
#ifdef lv
        double printscore = (viterbi.score / -100.0);
#else
        double printscore = viterbi.score;
#endif

        // printf("%s (%.2f)\n", string(result).c_str(), printscore);
        // printf("Zuker suboptimal structures...\n");
        if (seq_length < 500)
            window_size = 2;
        else if (seq_length < 1200)
            window_size = 5;
        else if (seq_length < 2000)
            window_size = 7;
        else if (seq_length < 5000)
            window_size = 9;
        else
            window_size = 12;
        window_size = 0; // Tianshuo: disable window_size
        bestP_beta.clear();
        bestP_beta.resize(seq_length);
        bestM2_beta.clear();
        bestM2_beta.resize(seq_length);
        bestM_beta.clear();
        bestM_beta.resize(seq_length);
        bestC_beta.clear();
        bestC_beta.resize(seq_length);
        bestMulti_beta.clear();
        bestMulti_beta.resize(seq_length);

        outside(next_pair);      
        
        bool sorted_P = true;

        map<tuple<char, int, int>, string> visited; // this is for recording the structure to bottom for a specific P, C,... at i,j

        clock_t startTime;
        startTime = clock();

        // state P is sorted
        vector<tuple<value_type, int, int> > sorted_bestP_beta;
        map<tuple<BestTypes, int, int>, pair<string, string> > global_visited_outside;
        map<tuple<BestTypes, int, int>, string> global_visited_inside;
        set<pair<int,int> > window_visited;

        for(int j = 0; j < seq_length; ++j) {
            for (auto & i_elem : bestP_beta[j]){
                auto i = i_elem.first;
                if (bestP_beta[j][i].manner == 0) continue;
                value_type alpha_inside = bestP[j][i].score;
                value_type beta_outside = bestP_beta[j][i].score;
#ifdef lv
                if (abs(alpha_inside + beta_outside - float(viterbi.score))/100. > zuker_energy_delta)
#else
                if (abs(alpha_inside + beta_outside - float(viterbi.score)) > zuker_energy_delta)
#endif
                {
                    continue;
                }
                sorted_bestP_beta.push_back(make_tuple(-(alpha_inside + beta_outside), i, j));
            }
        }

        sort(sorted_bestP_beta.begin(), sorted_bestP_beta.end(), cmp);

        int count = 0;

        int num_outputs = 0;

        for (auto & item : sorted_bestP_beta){
    
            auto i = get<1>(item);
            auto j = get<2>(item);

            auto best_P_ij_outside = bestP_beta[j][i];

            auto best_P_ij_inside = bestP[j][i];

            if (window_visited.find(make_pair(i,j)) != window_visited.end()) {
                count += 1;
                continue;
            }

            visited.clear();

            window_fill(window_visited, i, j, seq_length, window_size);

            auto outsider =  get_parentheses_outside_real_backtrace(i, j, best_P_ij_outside, global_visited_outside, global_visited_inside, window_visited);

            global_visited_outside[make_tuple(TYPE_P, i, j)] = outsider;

            auto insider =  get_parentheses_inside_real_backtrace(i, j, best_P_ij_inside, global_visited_inside, window_visited);

            global_visited_inside[make_tuple(TYPE_P, i, j)] = insider;

            auto second_string = outsider.first + insider + outsider.second;

            num_outputs+=1;
            if (num_outputs > 10000 && seq_length > 20000)
                break;

// #ifdef lv
//                 printf("%s (%.2f)\n", second_string.c_str(), (bestP[j][i].score + bestP_beta[j][i].score)/(-100.));
// #else
//                 printf("%s (%.2f)\n", second_string.c_str(), bestP[j][i].score + bestP_beta[j][i].score);
// #endif
        }
    }
    return {string(result), viterbi.score, nos_tot, parse_elapsed_time};
}

BeamCKYParser::DecoderResult BeamCKYParser::parse(string& seq, vector<int>* cons, vector<string>& subopts) {

    struct timeval parse_starttime, parse_endtime;

    // number of states
    unsigned long nos_H = 0, nos_P = 0, nos_M2 = 0,
            nos_M = 0, nos_C = 0, nos_Multi = 0;
    gettimeofday(&parse_starttime, NULL);
    prepare(static_cast<unsigned>(seq.length()));

    for (int i = 0; i < seq_length; ++i)
        nucs[i] = GET_ACGU_NUM(seq[i]);

    // lisiz, constraints
    if (use_constraints) {
        for (int i=0; i<seq_length; i++){
            int cons_idx = (*cons)[i];
            allow_unpaired_position[i] = cons_idx == -1 || cons_idx == -2;
            if (cons_idx > -1){
                if (!_allowed_pairs[nucs[i]][nucs[cons_idx]]){
                    printf("Constrains on non-classical base pairs (non AU, CG, GU pairs)\n");
                    exit(1);
                }
            }
        }
        int firstpair = seq_length;
        for (int i=seq_length-1; i>-1; i--){
            allow_unpaired_range[i] = firstpair;
            if ((*cons)[i] >= 0)
                firstpair = i;
        }
    }

    vector<int> next_pair[NOTON];
    {
        if (use_constraints){
            for (int nuci = 0; nuci < NOTON; ++nuci) {
                next_pair[nuci].resize(seq_length, -1);
                int next = -1;
                for (int j = seq_length-1; j >=0; --j) {
                    next_pair[nuci][j] = next;
                    if ((*cons)[j] > -2 && _allowed_pairs[nuci][nucs[j]]) next = j;
                }
            }
        } else {
            for (int nuci = 0; nuci < NOTON; ++nuci) {
                next_pair[nuci].resize(seq_length, -1);
                int next = -1;
                for (int j = seq_length-1; j >=0; --j) {
                    next_pair[nuci][j] = next;
                    if (_allowed_pairs[nuci][nucs[j]]) next = j;
                }
            }
        }
    }

#ifdef SPECIAL_HP
#ifdef lv
    v_init_tetra_hex_tri(seq, seq_length, if_tetraloops, if_hexaloops, if_triloops);
#else
    if (is_verbose)
        v_init_tetra_hex_tri(seq, seq_length, if_tetraloops, if_hexaloops, if_triloops);
#endif
#endif

    // start CKY decoding
#ifdef lv
        if(seq_length > 0) bestC[0].set(- v_score_external_unpaired(0, 0), MANNER_C_eq_C_plus_U);
        if(seq_length > 1) bestC[1].set(- v_score_external_unpaired(0, 1), MANNER_C_eq_C_plus_U);
#else
        if(seq_length > 0) bestC[0].set(score_external_unpaired(0, 0), MANNER_C_eq_C_plus_U);
        if(seq_length > 1) bestC[1].set(score_external_unpaired(0, 1), MANNER_C_eq_C_plus_U);
#endif
    ++nos_C;

    // from left to right
    for(int j = 0; j < seq_length; ++j) {
        int nucj = nucs[j];
        int nucj1 = (j+1) < seq_length ? nucs[j+1] : -1;

        unordered_map<int, State>& beamstepH = bestH[j];
        unordered_map<int, State>& beamstepMulti = bestMulti[j];
        unordered_map<int, State>& beamstepP = bestP[j];
        unordered_map<int, State>& beamstepM2 = bestM2[j];
        unordered_map<int, State>& beamstepM = bestM[j];
        State& beamstepC = bestC[j];

        // beam of H
        {
            if (beam > 0 && beamstepH.size() > beam) beam_prune(beamstepH);

            {
                int jnext = next_pair[nucj][j];
                if (no_sharp_turn) while (jnext - j < 4 && jnext != -1) jnext = next_pair[nucj][jnext];

                // lisiz, constriants
                if (use_constraints){
                    if (!allow_unpaired_position[j]){
                        jnext = (*cons)[j] > j ? (*cons)[j] : -1; // lisiz: j must be left bracket, jump to the constrainted pair (j, j') directly
                    }
                    if (jnext != -1){
                        int nucjnext = nucs[jnext];
                        if (jnext > allow_unpaired_range[j] || !allow_paired(j, jnext, cons, nucj, nucjnext))  // lisiz: avoid cross constrainted brackets or unallowed pairs
                            jnext = -1;
                    }
                }

                if (jnext != -1) {
                    int nucjnext = nucs[jnext];
                    int nucjnext_1 = (jnext - 1) > -1 ? nucs[jnext - 1] : -1;

                    value_type newscore;

#ifdef lv
                        int tetra_hex_tri = -1;
#ifdef SPECIAL_HP
                        if (jnext-j-1 == 4) // 6:tetra
                            tetra_hex_tri = if_tetraloops[j];
                        else if (jnext-j-1 == 6) // 8:hexa
                            tetra_hex_tri = if_hexaloops[j];
                        else if (jnext-j-1 == 3) // 5:tri
                            tetra_hex_tri = if_triloops[j];
#endif
                        newscore = - v_score_hairpin(j, jnext, nucj, nucj1, nucjnext_1, nucjnext, tetra_hex_tri);
#else
                        newscore = score_hairpin(j, jnext, nucj, nucj1, nucjnext_1, nucjnext);
#endif
                    // this candidate must be the best one at [j, jnext]
                    // so no need to check the score
                    update_if_better(bestH[jnext][j], newscore, MANNER_H);
                    ++ nos_H;
                }
            }

            {
                // for every state h in H[j]
                //   1. extend h(i, j) to h(i, jnext)
                //   2. generate p(i, j)
#ifdef lv
                sort_keys(beamstepH, keys);
                for (auto &item : keys) {
#else
                for (auto &item : beamstepH) {
#endif
                    int i = item.first;
                    // printf("%d\n", i);
                    State &state = item.second;
                    int nuci = nucs[i];
                    int jnext = next_pair[nuci][j];

                    // 2. generate p(i, j)
                    // lisiz, change the order because of the constriants
                    {
                        update_if_better(beamstepP[i], state.score, MANNER_HAIRPIN);
                        ++ nos_P;
                    }

                    // lisiz, constraints
                    if (jnext != -1 && use_constraints){
                        int nucjnext = nucs[jnext];
                        if (jnext > allow_unpaired_range[i] || !allow_paired(i, jnext, cons, nuci, nucjnext))
                            continue;
                    }

                    if (jnext != -1) {
                        int nuci1 = (i + 1) < seq_length ? nucs[i + 1] : -1;
                        int nucjnext = nucs[jnext];
                        int nucjnext_1 = (jnext - 1) > -1 ? nucs[jnext - 1] : -1;

                        // 1. extend h(i, j) to h(i, jnext)
                        value_type newscore;

#ifdef lv
                            int tetra_hex_tri = -1;
#ifdef SPECIAL_HP
                            if (jnext-i-1 == 4) // 6:tetra
                                tetra_hex_tri = if_tetraloops[i];
                            else if (jnext-i-1 == 6) // 8:hexa
                                tetra_hex_tri = if_hexaloops[i];
                            else if (jnext-i-1 == 3) // 5:tri
                                tetra_hex_tri = if_triloops[i];
#endif
                            newscore = - v_score_hairpin(i, jnext, nuci, nuci1, nucjnext_1, nucjnext, tetra_hex_tri);
#else
                            newscore = score_hairpin(i, jnext, nuci, nuci1, nucjnext_1, nucjnext);
#endif
                        // this candidate must be the best one at [i, jnext]
                        // so no need to check the score
                        update_if_better(bestH[jnext][i], newscore, MANNER_H);
                        ++nos_H;
                    }
                }
            }
        }
        if (j == 0) continue;

        // beam of Multi
        {
            if (beam > 0 && beamstepMulti.size() > beam) beam_prune(beamstepMulti);

            // for every state in Multi[j]
            //   1. extend (i, j) to (i, jnext)
            //   2. generate P (i, j)
#ifdef lv
            sort_keys(beamstepMulti, keys);
            for (auto &item : keys) {
#else
            for(auto& item : beamstepMulti) {
#endif
                int i = item.first;
                State& state = item.second;
                int nuci = nucs[i];
                int nuci1 = nucs[i+1];
                int jnext = next_pair[nuci][j];

                // 2. generate P (i, j)
                // lisiz, change the order because of the constraits
                {
                    value_type newscore;
#ifdef lv
                        newscore = state.score - v_score_multi(i, j, nuci, nuci1, nucs[j-1], nucj, seq_length, dangle_model);
#else
                        newscore = state.score + score_multi(i, j, nuci, nuci1, nucs[j-1], nucj, seq_length);
#endif
                    update_if_better(beamstepP[i], newscore, MANNER_P_eq_MULTI);
                    ++ nos_P;
                }

                // lisiz cnstriants
                if (jnext != -1 && use_constraints){
                    int nucjnext = nucs[jnext];
                    if (jnext > allow_unpaired_range[j] || !allow_paired(i, jnext, cons, nuci, nucjnext))
                        continue;
                }

                // 1. extend (i, j) to (i, jnext)
                {
                    char new_l1 = state.trace.paddings.l1;
                    int new_l2 = state.trace.paddings.l2 + jnext - j;
                    // if (jnext != -1 && new_l1 + new_l2 <= SINGLE_MAX_LEN) {
                    if (jnext != -1) {
                        // 1. extend (i, j) to (i, jnext)
                        value_type newscore;
#ifdef lv
                            newscore = state.score - v_score_multi_unpaired(j, jnext - 1);
#else
                            newscore = state.score + score_multi_unpaired(j, jnext - 1);
#endif
                        // this candidate must be the best one at [i, jnext]
                        // so no need to check the score
                        update_if_better(bestMulti[jnext][i], newscore, MANNER_MULTI_eq_MULTI_plus_U,
                                         new_l1,
                                         new_l2
                        );
                        ++nos_Multi;
                    }
                }
            }
        }

        // beam of P
        {
            if (beam > 0 && beamstepP.size() > beam) beam_prune(beamstepP);

            // for every state in P[j]
            //   1. generate new helix/bulge
            //   2. M = P
            //   3. M2 = M + P
            //   4. C = C + P
#ifdef is_cube_pruning
            bool use_cube_pruning = beam > MIN_CUBE_PRUNING_SIZE
                                    && beamstepP.size() > MIN_CUBE_PRUNING_SIZE;
#else
            bool use_cube_pruning = false;
#endif               

#ifdef lv
            sort_keys(beamstepP, keys);
            for (auto &item : keys) {
#else
            for(auto& item : beamstepP) {
#endif
                int i = item.first;
                State& state = item.second;
                int nuci = nucs[i];
                int nuci_1 = (i-1>-1) ? nucs[i-1] : -1;

                // 2. M = P // lhuang: check j (this P is not the last non-closing P in multi)
                if(i > 0 && j < seq_length-6){
                    value_type newscore;
#ifdef lv
                        newscore = - v_score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length, dangle_model) + state.score;
#else
                        newscore = score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length) + state.score;
#endif
                    update_if_better(beamstepM[i], newscore, MANNER_M_eq_P);
                    ++ nos_M;
                }
                //printf(" M = P at %d\n", j); fflush(stdout);

                // 3. M2 = M + P // lhuang: check j < n-1
                if(!use_cube_pruning && j < seq_length - 1) {
                    int k = i - 1;
                    if ( k > 0 && !bestM[k].empty()) {
                        value_type M1_score;
#ifdef lv
                            M1_score = - v_score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length, dangle_model) + state.score;
#else
                            M1_score = score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length) + state.score;
#endif
                        // candidate list
                        auto bestM2_iter = beamstepM2.find(i);
#ifndef is_candidate_list
                        for (auto &m : bestM[k]) {
                                int newi = m.first;
                                // eq. to first convert P to M1, then M2/M = M + M1
                                value_type newscore = M1_score + m.second.score;
                                update_if_better(beamstepM2[newi], newscore, MANNER_M2_eq_M_plus_P, k);
                                ++nos_M2;
                            }
#else
                        if (bestM2_iter==beamstepM2.end() || M1_score > bestM2_iter->second.score) {
                            for (auto &m : bestM[k]) {
                                int newi = m.first;
                                // eq. to first convert P to M1, then M2/M = M + M1
                                value_type newscore = M1_score + m.second.score;
                                update_if_better(beamstepM2[newi], newscore, MANNER_M2_eq_M_plus_P, k);
                                ++nos_M2;
                            }
                        }
#endif
                    }
                }
                //printf(" M/M2 = M + P at %d\n", j); fflush(stdout);

                // 4. C = C + P
                {
                    int k = i - 1;
                    if (k >= 0) {
                      State& prefix_C = bestC[k];
                      if (prefix_C.manner != MANNER_NONE) {
                        int nuck = nuci_1;
                        int nuck1 = nuci;
                        value_type newscore;
#ifdef lv
                            newscore = - v_score_external_paired(k+1, j, nuck, nuck1, nucj, nucj1, seq_length, dangle_model) +
                                prefix_C.score + state.score;
#else
                            newscore = score_external_paired(k+1, j, nuck, nuck1, nucj, nucj1, seq_length) +
                                prefix_C.score + state.score;
#endif
                        update_if_better(beamstepC, newscore, MANNER_C_eq_C_plus_P, k);
                        ++ nos_C;
                      }
                    } else {
                        value_type newscore;
#ifdef lv
                            newscore = - v_score_external_paired(0, j, -1, nucs[0], nucj, nucj1, seq_length, dangle_model) + state.score;
#else
                            newscore = score_external_paired(0, j, -1, nucs[0], nucj, nucj1, seq_length) + state.score;
#endif
                        update_if_better(beamstepC, newscore, MANNER_C_eq_C_plus_P, -1);
                        ++ nos_C;
                    }
                }
                //printf(" C = C + P at %d\n", j); fflush(stdout);

                // 1. generate new helix / single_branch
                // new state is of shape p..i..j..q
                if (i >0 && j<seq_length-1) {
                    value_type precomputed;
#ifdef lv
                        precomputed = 0;
#else
                        precomputed = score_junction_B(j, i, nucj, nucj1, nuci_1, nuci);
#endif
                    for (int p = i - 1; p >= std::max(i - SINGLE_MAX_LEN, 0); --p) {
                        int nucp = nucs[p];
                        int nucp1 = nucs[p + 1]; // hzhang: move here
                        int q = next_pair[nucp][j];

                        // lisiz constraints
                        if (use_constraints){
                            if (p < i-1 && !allow_unpaired_position[p+1]) // lisiz: if p+1 must be paired, break
                                break;
                            if (!allow_unpaired_position[p]){             // lisiz: if p must be paired, p must be left bracket
                                q = (*cons)[p];
                                if (q < p) break;
                            }
                        }

                        while (q != -1 && ((i - p) + (q - j) - 2 <= SINGLE_MAX_LEN)) {
                            int nucq = nucs[q];

                            // lisiz constraints
                            if (use_constraints){
                                if (q>j+1 && q > allow_unpaired_range[j])  // lisiz: if q-1 must be paired, break
                                    break;
                                if (!allow_paired(p, q, cons, nucp, nucq)) // lisiz: if p q are )(, break
                                    break;
                            }

                            int nucq_1 = nucs[q - 1];
                            if (p == i - 1 && q == j + 1) {
                                // helix
                                value_type newscore;
#ifdef lv
                                    newscore = -v_score_single(p,q,i,j, nucp, nucp1, nucq_1, nucq, nuci_1, nuci, nucj, nucj1)
                                        + state.score;
                                    // SHAPE for Vienna only
                                    if (use_shape)
                                        newscore += -(pseudo_energy_stack[p] + pseudo_energy_stack[i] + pseudo_energy_stack[j] + pseudo_energy_stack[q]);
#else
                                    newscore = score_helix(nucp, nucp1, nucq_1, nucq) + state.score;
#endif
                                update_if_better(bestP[q][p], newscore, MANNER_HELIX);
                                ++nos_P;
                            } else {
                                // single branch
                                value_type newscore;
#ifdef lv
                                    newscore = - v_score_single(p,q,i,j, nucp, nucp1, nucq_1, nucq, nuci_1, nuci, nucj, nucj1)
                                        + state.score;
#else
                                    newscore = score_junction_B(p, q, nucp, nucp1, nucq_1, nucq) +
                                        precomputed +
                                        score_single_without_junctionB(p, q, i, j, nuci_1, nuci, nucj, nucj1) +
                                        state.score;
#endif
                                update_if_better(bestP[q][p], newscore, MANNER_SINGLE,
                                                 static_cast<char>(i - p),
                                                 q - j);
                                ++nos_P;
                            }
                            q = next_pair[nucp][q];
                        }
                    }
                }
            }

	    // lhuang: check j < n-1
            if (use_cube_pruning && j < seq_length - 1) {
                // 3. M2 = M + P with cube pruning
                vector<int> valid_Ps;
                vector<value_type> M1_scores;
#ifdef lv
            sort_keys(beamstepP, keys);
            for (auto &item : keys) {
#else
            for(auto& item : beamstepP) {
#endif
                    int i = item.first;
                    State &state = item.second;
                    int nuci = nucs[i];
                    int nuci_1 = (i - 1 > -1) ? nucs[i - 1] : -1;
                    int k = i - 1;

                    // group candidate Ps
                    if (k > 0 && !bestM[k].empty()) {
                        assert(bestM[k].size() == sorted_bestM[k].size());
                        value_type M1_score;
#ifdef lv
                            M1_score = - v_score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length, dangle_model) + state.score;
#else
                            M1_score = score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length) + state.score;
#endif
                        auto bestM2_iter = beamstepM2.find(i);
#ifndef is_candidate_list
                        valid_Ps.push_back(i);
                        M1_scores.push_back(M1_score);
#else
                        if (bestM2_iter == beamstepM2.end() || M1_score > bestM2_iter->second.score) {
                            valid_Ps.push_back(i);
                            M1_scores.push_back(M1_score);
                        }
#endif
                    }
                }

                // build max heap
                // heap is of form (heuristic score, (index of i in valid_Ps, index of M in bestM[i-1]))
                vector<pair<value_type, pair<int, int>>> heap;
                for (int p = 0; p < valid_Ps.size(); ++p) {
                    int i = valid_Ps[p];
                    int k = i - 1;
                    heap.push_back(make_pair(M1_scores[p] + sorted_bestM[k][0].first,
                                             make_pair(p, 0)
                    ));
                    push_heap(heap.begin(), heap.end());
                }

                // start cube pruning
                // stop after beam size M2 states being filled
                int filled = 0;
                // exit when filled >= beam and current score < prev score
                value_type prev_score = VALUE_MIN;
                value_type current_score = VALUE_MIN;
                while ((filled < beam || current_score == prev_score) && !heap.empty()) {
                    auto &top = heap.front();
                    prev_score = current_score;
                    current_score = top.first;
                    int index_P = top.second.first;
                    int index_M = top.second.second;
                    int i = valid_Ps[top.second.first];
                    int k = i - 1;
                    int newi = sorted_bestM[k][index_M].second;
                    value_type newscore = M1_scores[index_P] + bestM[k][newi].score;
                    pop_heap(heap.begin(), heap.end());
                    heap.pop_back();

                    if (beamstepM2[newi].manner == MANNER_NONE) {
                        ++filled;
                        update_if_better(beamstepM2[newi], newscore, MANNER_M2_eq_M_plus_P, k);
                        ++nos_M2;
                    } else {
                        assert(beamstepM2[newi].score > newscore - 1e-8);
                    }

                    ++index_M;
                    while (index_M < sorted_bestM[k].size()) {
                        // candidate_score is a heuristic score
                        value_type candidate_score = M1_scores[index_P] + sorted_bestM[k][index_M].first;
                        int candidate_newi = sorted_bestM[k][index_M].second;
                        if (beamstepM2.find(candidate_newi) == beamstepM2.end()) {
                            heap.push_back(make_pair(candidate_score,
                                                     make_pair(index_P, index_M)));
                            push_heap(heap.begin(), heap.end());
                            break;
                        } else {
                            // based on the property of cube pruning, the new score must be worse
                            // than the state already inserted
                            // so we keep iterate through the candidate list to find the next
                            // candidate
                            ++index_M;
                            assert(beamstepM2[candidate_newi].score >
                                   M1_scores[index_P] + bestM[k][candidate_newi].score - 1e-8);
                        }
                    }
                }
            }
        }

        // beam of M2
        {
            if (beam > 0 && beamstepM2.size() > beam) beam_prune(beamstepM2);

            // for every state in M2[j]
            //   1. multi-loop  (by extending M2 on the left)
            //   2. M = M2
#ifdef lv
            sort_keys(beamstepM2, keys);
            for (auto &item : keys) {
#else
            for(auto& item : beamstepM2) {
#endif
                int i = item.first;
                State& state = item.second;

                // 2. M = M2
                {
                    update_if_better(beamstepM[i], state.score, MANNER_M_eq_M2);
                    ++ nos_M;
                }

                // 1. multi-loop
                {
                    for (int p = i-1; p >= std::max(i - SINGLE_MAX_LEN, 0); --p) {
                        int nucp = nucs[p];
                        int q = next_pair[nucp][j];

                        if (use_constraints){
                            if (p < i - 1 && !allow_unpaired_position[p+1])
                                break;
                            if (!allow_unpaired_position[p]){
                                q = (*cons)[p];
                                if (q < p) break;
                            }
                            if (q > j+1 && q > allow_unpaired_range[j])
                                continue;
                            int nucq = nucs[q];
                            if (!allow_paired(p, q, cons, nucp, nucq))
                                continue;
                        }

                        if (q != -1 && ((i - p - 1) <= SINGLE_MAX_LEN)) {
                            // the current shape is p..i M2 j ..q

                            value_type newscore;
#ifdef lv
                                newscore = - v_score_multi_unpaired(p+1, i-1) -
                                    v_score_multi_unpaired(j+1, q-1) + state.score;
#else
                                newscore = score_multi_unpaired(p+1, i-1) +
                                    score_multi_unpaired(j+1, q-1) + state.score;
#endif
                            update_if_better(bestMulti[q][p], newscore, MANNER_MULTI,
                                             static_cast<char>(i - p),
                                             q - j);
                            ++ nos_Multi;
                        }
                    }
                }
            }
        }

        // beam of M
        {
            value_type threshold = VALUE_MIN;
            if (beam > 0 && beamstepM.size() > beam) threshold = beam_prune(beamstepM);

#ifdef is_cube_pruning
                sortM(threshold, beamstepM, sorted_bestM[j]);
            // }
#endif

            // for every state in M[j]
            //   1. M = M + unpaired
#ifdef lv
            sort_keys(beamstepM, keys);
            for (auto &item : keys) {
#else
            for(auto& item : beamstepM) {
#endif
                int i = item.first;
                State& state = item.second;
                if (j < seq_length-1) {
                    if (use_constraints && !allow_unpaired_position[j+1]) // if j+1 must be paired
                        continue;
                    value_type newscore;
#ifdef lv
                        newscore = - v_score_multi_unpaired(j + 1, j + 1) + state.score;
#else
                        newscore = score_multi_unpaired(j + 1, j + 1) + state.score;
#endif
                    update_if_better(bestM[j+1][i], newscore, MANNER_M_eq_M_plus_U);
                    ++ nos_M;
                }
            }
        }

        // beam of C
        {
            // C = C + U
            if (j < seq_length-1) {
                if (use_constraints && !allow_unpaired_position[j+1])
                        continue;
                value_type newscore;
#ifdef lv
                    newscore = -v_score_external_unpaired(j+1, j+1) + beamstepC.score;
#else
                    newscore = score_external_unpaired(j+1, j+1) + beamstepC.score;
#endif
                update_if_better(bestC[j+1], newscore, MANNER_C_eq_C_plus_U);
                ++ nos_C;
            }
        }

    }  // end of for-loo j

    State& viterbi = bestC[seq_length-1];

    char result[seq_length+1];
    get_parentheses(result, seq);

    gettimeofday(&parse_endtime, NULL);
    double parse_elapsed_time = parse_endtime.tv_sec - parse_starttime.tv_sec + (parse_endtime.tv_usec-parse_starttime.tv_usec)/1000000.0;

    unsigned long nos_tot = nos_H + nos_P + nos_M2 + nos_Multi + nos_M + nos_C;
    if (is_verbose) {
        printf("Parse Time: %f len: %d score %f #states %lu H %lu P %lu M2 %lu Multi %lu M %lu C %lu\n",
               parse_elapsed_time, seq_length, double(viterbi.score), nos_tot,
               nos_H, nos_P, nos_M2, nos_Multi, nos_M, nos_C);
    }

    fflush(stdout);

    if (zuker){
#ifdef lv
        double printscore = (viterbi.score / -100.0);
#else
        double printscore = viterbi.score;
#endif

        // printf("%s (%.2f)\n", string(result).c_str(), printscore);
        // subopts.push_back(string(result));
        // printf("Zuker suboptimal structures...\n");
        if (seq_length < 500)
            window_size = 2;
        else if (seq_length < 1200)
            window_size = 5;
        else if (seq_length < 2000)
            window_size = 7;
        else if (seq_length < 5000)
            window_size = 9;
        else
            window_size = 12;

        bestP_beta.clear();
        bestP_beta.resize(seq_length);
        bestM2_beta.clear();
        bestM2_beta.resize(seq_length);
        bestM_beta.clear();
        bestM_beta.resize(seq_length);
        bestC_beta.clear();
        bestC_beta.resize(seq_length);
        bestMulti_beta.clear();
        bestMulti_beta.resize(seq_length);

        outside(next_pair);      
        
        bool sorted_P = true;

        map<tuple<char, int, int>, string> visited; // this is for recording the structure to bottom for a specific P, C,... at i,j

        clock_t startTime;
        startTime = clock();

        // state P is sorted
        vector<tuple<value_type, int, int> > sorted_bestP_beta;
        map<tuple<BestTypes, int, int>, pair<string, string> > global_visited_outside;
        map<tuple<BestTypes, int, int>, string> global_visited_inside;
        set<pair<int,int> > window_visited;

        for(int j = 0; j < seq_length; ++j) {
            for (auto & i_elem : bestP_beta[j]){
                auto i = i_elem.first;
                if (bestP_beta[j][i].manner == 0) continue;
                value_type alpha_inside = bestP[j][i].score;
                value_type beta_outside = bestP_beta[j][i].score;
#ifdef lv
                if (abs(alpha_inside + beta_outside - float(viterbi.score))/100. > zuker_energy_delta)
#else
                if (abs(alpha_inside + beta_outside - float(viterbi.score)) > zuker_energy_delta)
#endif
                {
                    continue;
                }
                sorted_bestP_beta.push_back(make_tuple(-(alpha_inside + beta_outside), i, j));
            }
        }

        sort(sorted_bestP_beta.begin(), sorted_bestP_beta.end(), cmp);

        int count = 0;

        int num_outputs = 0;

        for (auto & item : sorted_bestP_beta){
    
            auto i = get<1>(item);
            auto j = get<2>(item);

            auto best_P_ij_outside = bestP_beta[j][i];

            auto best_P_ij_inside = bestP[j][i];

            if (window_visited.find(make_pair(i,j)) != window_visited.end()) {
                count += 1;
                continue;
            }

            visited.clear();

            window_fill(window_visited, i, j, seq_length, window_size);

            auto outsider =  get_parentheses_outside_real_backtrace(i, j, best_P_ij_outside, global_visited_outside, global_visited_inside, window_visited);

            global_visited_outside[make_tuple(TYPE_P, i, j)] = outsider;

            auto insider =  get_parentheses_inside_real_backtrace(i, j, best_P_ij_inside, global_visited_inside, window_visited);

            global_visited_inside[make_tuple(TYPE_P, i, j)] = insider;

            auto second_string = outsider.first + insider + outsider.second;

            num_outputs+=1;
            if (num_outputs > 10000 && seq_length > 20000)
                break;
            subopts.push_back(second_string);

// #ifdef lv
//                 printf("%s (%.2f)\n", second_string.c_str(), (bestP[j][i].score + bestP_beta[j][i].score)/(-100.));
// #else
//                 printf("%s (%.2f)\n", second_string.c_str(), bestP[j][i].score + bestP_beta[j][i].score);
// #endif
        }
    }
    return {string(result), viterbi.score, nos_tot, parse_elapsed_time};
}


void BeamCKYParser::outside(vector<int> next_pair[]){
      
    struct timeval parse_starttime, parse_endtime;

    gettimeofday(&parse_starttime, NULL);

    bestC_beta[seq_length-1].score = 0.0;
    bestC_beta[seq_length-1].manner = MANNER_NONE;
    // from right to left
    value_type newscore;
    for(int j = seq_length-1; j >= 0; --j) {
        int nucj = nucs[j];
        int nucj1 = (j+1) < seq_length ? nucs[j+1] : -1;

        // important otherwise cannot iterate    
        unordered_map<int, State>& beamstepH = bestH[j];
        unordered_map<int, State>& beamstepMulti = bestMulti[j];
        unordered_map<int, State>& beamstepP = bestP[j];
        unordered_map<int, State>& beamstepM2 = bestM2[j];
        unordered_map<int, State>& beamstepM = bestM[j];
        State& beamstepC = bestC[j];

        unordered_map<int, State>& beamstepH_beta = bestH_beta[j];
        unordered_map<int, State>& beamstepMulti_beta = bestMulti_beta[j];
        unordered_map<int, State>& beamstepP_beta = bestP_beta[j];
        unordered_map<int, State>& beamstepM2_beta = bestM2_beta[j];
        unordered_map<int, State>& beamstepM_beta = bestM_beta[j];
        State& beamstepC_beta = bestC_beta[j];

        // beam of C
        {
            // C = C + U
            if (j < seq_length-1) {

#ifdef lv
                update_if_better(beamstepC_beta, bestC_beta[j+1].score, MANNER_C_eq_C_plus_U);           
#else
                newscore = score_external_unpaired(j+1, j+1);
                update_if_better(beamstepC_beta, bestC_beta[j+1].score + newscore, MANNER_C_eq_C_plus_U);

#endif
            }
        }

        if (j == 0) break; //for updating bestC_beta[0], then break need: for(int j = seq_length-1; |||||j >= 0||||; --j) {
    
        // beam of M
        {
            for(auto& item : beamstepM) {

                int i = item.first;
                State& state = beamstepM_beta[i];


                if (j < seq_length-1) {
#ifdef lv            
                    if (bestM_beta[j+1][i].manner != 0){
                        update_if_better(state, bestM_beta[j+1][i].score, MANNER_M_eq_M_plus_U);
                    }
#else
                    if (bestM_beta[j+1][i].manner != 0){
                        newscore = score_multi_unpaired(j + 1, j + 1);
                        update_if_better(state, bestM_beta[j+1][i].score + newscore, MANNER_M_eq_M_plus_U);
                    }
#endif
                }
            }
        }

        // beam of M2
        {
            for(auto& item : beamstepM2) {
                int i = item.first;
                State& state = beamstepM2_beta[i];

                // 1. multi-loop
                {
                    for (int p = i-1; p >= std::max(i - SINGLE_MAX_LEN, 0); --p) {
                        int nucp = nucs[p];
                        int q = next_pair[nucp][j];
                        if (q != -1 && ((i - p - 1) <= SINGLE_MAX_LEN)) {
#ifdef lv
                            if (bestMulti_beta[q][p].manner != 0){
                                update_if_better(state, bestMulti_beta[q][p].score, MANNER_MULTI, static_cast<char>(i - p), q - j);
                            }
#else
     
                            if (bestMulti_beta[q][p].manner != 0){
                                newscore = score_multi_unpaired(p+1, i-1) + score_multi_unpaired(j+1, q-1);
                                update_if_better(state, bestMulti_beta[q][p].score + newscore, MANNER_MULTI, static_cast<char>(i - p), q - j);
                            }
#endif
                        }
                    }
                }

                // 2. M = M2
                if (beamstepM_beta[i].manner != 0){
                    update_if_better(state, beamstepM_beta[i].score, MANNER_M_eq_M2);
                }
            }
        }
        // beam of P
        {  
            for(auto& item : beamstepP) {
                int i = item.first;
                State& state = beamstepP_beta[i];
                auto state_alpha = item.second.score; //for subopt zuker


                int nuci = nucs[i];
                int nuci_1 = (i-1>-1) ? nucs[i-1] : -1;

                if (i >0 && j<seq_length-1) {
#ifndef lv
                    value_type precomputed = score_junction_B(j, i, nucj, nucj1, nuci_1, nuci);
#endif
                    for (int p = i - 1; p >= std::max(i - SINGLE_MAX_LEN, 0); --p) {
                        int nucp = nucs[p];
                        int nucp1 = nucs[p + 1]; 
                        int q = next_pair[nucp][j];

                        while (q != -1 && ((i - p) + (q - j) - 2 <= SINGLE_MAX_LEN)) {
                            int nucq = nucs[q];
                            int nucq_1 = nucs[q - 1];

                            if (p == i - 1 && q == j + 1) {
                                // helix
#ifdef lv
                                
                                if (bestP_beta[q][p].manner != 0){
                                    int score_single = -v_score_single(p,q,i,j, nucp, nucp1, nucq_1, nucq, nuci_1, nuci, nucj, nucj1);
                                    
                                    // SHAPE for Vienna only
                                    if (use_shape)
                                    {
                                        score_single += -(pseudo_energy_stack[p] + pseudo_energy_stack[i] + pseudo_energy_stack[j] + pseudo_energy_stack[q]);
                                    }

                                    update_if_better(state, (bestP_beta[q][p].score + score_single), MANNER_HELIX, static_cast<char>(i - p), q - j);
                                }
#else
                                if (bestP_beta[q][p].manner != 0){    
                                    newscore = score_helix(nucp, nucp1, nucq_1, nucq);
                                    update_if_better(state, bestP_beta[q][p].score + newscore, MANNER_HELIX, static_cast<char>(i - p), q - j);
                                }

#endif
                            } else {
                                // single branch
#ifdef lv
                                if  (bestP_beta[q][p].manner != 0){   
                                    int score_single = - v_score_single(p,q,i,j, nucp, nucp1, nucq_1, nucq, nuci_1, nuci, nucj, nucj1);
                                    update_if_better(state, (bestP_beta[q][p].score + score_single), MANNER_SINGLE, static_cast<char>(i - p), q - j);
                                }
#else
                                if (bestP_beta[q][p].manner != 0){    
                                    newscore = score_junction_B(p, q, nucp, nucp1, nucq_1, nucq) + precomputed + score_single_without_junctionB(p, q, i, j, nuci_1, nuci, nucj, nucj1);
                                    update_if_better(state, bestP_beta[q][p].score + newscore, MANNER_SINGLE, static_cast<char>(i - p), q - j);
                                }

#endif
                            }
                            q = next_pair[nucp][q]; //beam search, q is not a key in bestP_beta[q], does not mean new q is not a key
                        }
                    }
                }

                // assert(check());

                // 2. M = P
                if(i > 0 && j < seq_length-1){

#ifdef lv
                        if (beamstepM_beta[i].manner != 0){
                            int score_M1 = - v_score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length, dangle_model);
                            update_if_better(state, (beamstepM_beta[i].score + score_M1), MANNER_M_eq_P);
                        }
#else
                        if (beamstepM_beta[i].manner != 0){    
                            newscore = score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length);
                            update_if_better(state, beamstepM_beta[i].score + newscore, MANNER_M_eq_P);
                        }
#endif
                }
                // 3. M2 = M + P
                int k = i - 1;
                if ( k > 0 && !bestM[k].empty()) {
#ifdef lv
                    int M1_score = - v_score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length, dangle_model);
                    float m1_alpha = M1_score;
#else
                    newscore = score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length);
                    float m1_alpha = newscore;
#endif
                    
                    float m1_plus_P_alpha = state_alpha + m1_alpha;
                    for (auto &m : bestM[k]) {
                        int newi = m.first;
                        State& m_state = bestM_beta[k][newi];

                        auto m_state_alpha = m.second.score;
                        if (beamstepM2_beta[newi].manner != 0){
                            update_if_better(m_state, (beamstepM2_beta[newi].score + m1_plus_P_alpha), MANNER_M2_eq_M_plus_P, j); //newi, k = i-1, i, j
                            update_if_better(state, (beamstepM2_beta[newi].score + m_state_alpha + m1_alpha), MANNER_M2_eq_M_plus_P, newi); // newi, k = i-1, i, j
                        }
                    }
                }
                // 4. C = C + P
                {   
                    int k = i - 1;
                    if (k >= 0) {
                        int nuck = nuci_1;
                        int nuck1 = nuci;

#ifdef lv
                        int score_external_paired = - v_score_external_paired(k+1, j, nuck, nuck1,
                                                                 nucj, nucj1, seq_length, dangle_model);

                        float external_paired_alpha_plus_beamstepC_beta = beamstepC_beta.score + score_external_paired;

#else
                        newscore = score_external_paired(k+1, j, nuck, nuck1, nucj, nucj1, seq_length);
                        float external_paired_alpha_plus_beamstepC_beta = beamstepC_beta.score + newscore;
#endif

                        update_if_better(bestC_beta[k], state_alpha + external_paired_alpha_plus_beamstepC_beta, MANNER_C_eq_C_plus_P, j); //0, k, i, j
                        update_if_better(state, bestC[k].score + external_paired_alpha_plus_beamstepC_beta, MANNER_C_eq_C_plus_P, -1); // 0, k, i, j

                    } else {

#ifdef lv
                        int score_external_paired = - v_score_external_paired(0, j, -1, nucs[0],
                                                                 nucj, nucj1, seq_length, dangle_model);
                        update_if_better(state, (beamstepC_beta.score + score_external_paired), MANNER_C_eq_C_plus_P, -1); //C does not exist
#else
                        newscore = score_external_paired(0, j, -1, nucs[0],
                                                             nucj, nucj1, seq_length);
                        update_if_better(state, beamstepC_beta.score + newscore, MANNER_C_eq_C_plus_P, -1); //C does not exist

#endif
                    }
                }
            }
        }

        // beam of Multi
        {
            for(auto& item : beamstepMulti) {
                int i = item.first;
                State& state = beamstepMulti_beta[i];

                int nuci = nucs[i];
                int nuci1 = nucs[i+1];
                int jnext = next_pair[nuci][j];

                // 1. extend (i, j) to (i, jnext)
                {
                    if (jnext != -1) {

#ifdef lv
                        if (bestMulti_beta[jnext][i].manner != 0){
                            update_if_better(state, (bestMulti_beta[jnext][i].score), MANNER_MULTI_eq_MULTI_plus_U, jnext - j);
                        }
#else
                        if (bestMulti_beta[jnext][i].manner != 0){
                            newscore = score_multi_unpaired(j, jnext - 1);
                            update_if_better(state, bestMulti_beta[jnext][i].score + newscore, MANNER_MULTI_eq_MULTI_plus_U, jnext - j);
                        }
#endif
                    }
                }

                // 2. generate P (i, j)
                {
#ifdef lv
                    if (beamstepP_beta[i].manner != 0){
                        int score_multi = - v_score_multi(i, j, nuci, nuci1, nucs[j-1], nucj, seq_length, dangle_model);
                        update_if_better(state, (beamstepP_beta[i].score + score_multi), MANNER_P_eq_MULTI);
                    }
#else    
                    if (beamstepP_beta[i].manner != 0){
                        newscore = score_multi(i, j, nuci, nuci1, nucs[j-1], nucj, seq_length);
                        update_if_better(state, beamstepP_beta[i].score + newscore, MANNER_P_eq_MULTI);
                    }

#endif
                }
            }
        }
        
    }  // end of for-loo j

    return;
}

BeamCKYParser::BeamCKYParser(int beam_size,
                             bool nosharpturn,
                             bool verbose,
                             bool constraints,
                             bool zuker_subopt,
                             float energy_delta,
                             string shape_file_path,
                             bool fasta,
                             int dangles)
    : beam(beam_size), 
      no_sharp_turn(nosharpturn), 
      is_verbose(verbose),
      use_constraints(constraints),
      zuker(zuker_subopt),
      zuker_energy_delta(energy_delta),
      is_fasta(fasta),
      dangle_model(dangles){
#ifdef lv
        initialize();
#else
        initialize();
        initialize_cachesingle();
#endif

    // SHAPE
    if (shape_file_path != ""){

        use_shape = true;
        int position;
        string data;

        double temp_after_mb_shape;

        ifstream in(shape_file_path);

        if (!in.good()){
            printf("Reading SHAPE file error!\n");
            assert(false);
        }


        while (!(in >> position >> data).fail()) {
            assert(int(position) == SHAPE_data.size() + 1);

            if (isdigit(int(data[0])) == 0){
                SHAPE_data.push_back(double((-1.000000)));
            }

            else {
                SHAPE_data.push_back(stod(data));
            }
        }


        for (int i = 0; i<SHAPE_data.size(); i++){
            temp_after_mb_shape = SHAPE_data[i] < 0 ? 0. : (m * log(SHAPE_data[i] + 1) + b);

            pseudo_energy_stack.push_back((int)roundf(temp_after_mb_shape * 100.));

            assert(pseudo_energy_stack.size() == i + 1 );

        }

    }

}

// -------------------------------------------------------------

// trim from end (in place)
static inline void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
        return !std::isspace(ch);
    }).base(), s.end());
}

