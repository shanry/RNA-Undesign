/*
 *eval.cpp*
 Evaluate the energy of a given RNA structure. Recording each energy unit.

 author: Wei Yu Tang, Tianshuo Zhou (based on code by He Zhang)
 edited by: 09/2023
*/

#include <iostream>
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

#include "Utils/utility_v.h"
#include "Utils/utility_v_max.h"
#include "eval.h"

#define BASE 1000

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
    printf("critical loops start\n");
    printf("%lu\n", critical_loops.size() + critical_internal.size());
    for (auto &item: critical_loops) {
        // print: is ref1, loop type, indices...
        printf("%d %d ", item.second.first, get<2>(item.first));
        vector<int> indexed_loop = {item.second.first, get<2>(item.first)};
        for (int &x: item.second.second) {
            printf("%d ", x);
            indexed_loop.push_back(x);
        }
        printf("\n");
        cr_loops.push_back(indexed_loop);
    }

    for (auto &item: critical_bulge) {
        // print: is ref1, loop type, indices...
        loops type = bulge;
        printf("%d %d ", item.second.first, type);
        for (int &x: item.second.second) {
            printf("%d ", x);
        }
        printf("\n");
    }

    for (auto &item: critical_internal) {
        // print: is ref1, loop type, indices...
        loops type = interior;
        printf("%d %d ", item.second.first, type);
        vector<int> indexed_loop = {item.second.first, type};
        for (int &x: item.second.second) {
            printf("%d ", x);
            indexed_loop.push_back(x);
        }
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

        for (auto &item: critical_internal) {
            printf("Internal (%d, %d), (%d, %d) ref%d : ", get<0>(item.first), get<3>(item.first), get<1>(item.first), get<2>(item.first), -item.second.first+2);
            for (int &x: item.second.second) {
                printf("%d ", x);
            }
            printf("\n");
        }
    }

    printf("critical positions: ");
    for (int x: critical_positions) {
        printf("%d, ", x);
    }
    printf("\n");
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
    printf("critical loops start\n");
    printf("%lu\n", critical_loops.size() + critical_internal.size());
    for (auto &item: critical_loops) {
        // print: is ref1, loop type, indices...
        printf("%d %d ", item.second.first, get<2>(item.first));
        vector<int> indexed_loop = {item.second.first, get<2>(item.first)};
        for (int &x: item.second.second) {
            printf("%d ", x);
            indexed_loop.push_back(x);
        }
        printf("\n");
        cr_loops.push_back(indexed_loop);
    }

    for (auto &item: critical_bulge) {
        // print: is ref1, loop type, indices...
        loops type = bulge;
        printf("%d %d ", item.second.first, type);
        for (int &x: item.second.second) {
            printf("%d ", x);
        }
        printf("\n");
    }

    for (auto &item: critical_internal) {
        // print: is ref1, loop type, indices...
        loops type = interior;
        printf("%d %d ", item.second.first, type);
        vector<int> indexed_loop = {item.second.first, type};
        for (int &x: item.second.second) {
            printf("%d ", x);
            indexed_loop.push_back(x);
        }
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

    printf("critical positions: ");
    for (int x: critical_positions) {
        printf("%d, ", x);
    }
    printf("\n");
    return cr_loops;
}

long diff_eval(string seq, vector<vector<int>>& cr_loops, bool is_verbose, int dangle_model) {
    int n = seq.length();
    
    // weiyu: Special Hairpin is currently off
    // vector<int> if_tetraloops;
    // vector<int> if_hexaloops;
    // vector<int> if_triloops;

    // v_init_tetra_hex_tri(seq, n, if_tetraloops, if_hexaloops, if_triloops); // calculate if_tetraloops, if_hexaloops, if_triloops
    
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

            // if (j-i-1 == 4) // 6:tetra
            //     tetra_hex_tri = if_tetraloops[i];
            // else if (j-i-1 == 6) // 8:hexa
            //     tetra_hex_tri = if_hexaloops[i];
            // else if (j-i-1 == 3) // 5:tri
            //     tetra_hex_tri = if_triloops[i];

            score = - v_score_hairpin(i, j, nuci, nuci1, nucj_1, nucj, tetra_hex_tri);

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

    printf("ref1 energy = %.2f, ref2 energy = %.2f\n", energy_ref1 / -100.0, energy_ref2 / -100.0);
    return energy_ref1 - energy_ref2;
}

long linear_eval(string seq, string ref, bool is_verbose, int dangle_model) {

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
                
                int newscore = - v_score_hairpin(i, j, nuci, nuci1, nucj_1, nucj, tetra_hex_tri);
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

int main2(int argc, char* argv[]){
    // Print the program name (argv[0])
    cout << "Program name: " << argv[0] << endl;

    // Check the number of arguments
    cout << "Number of arguments: " << argc << endl;

    // Print all command-line arguments
    for (int i = 1; i < argc; ++i) {
        cout << "Argument " << i << ": " << argv[i] << endl;
    }

    // string seq = "AAAAAAGGAAAAAAAAGCCCGAAAAGGGUGAAAAAGAACAGAAAAAAAAAAAAAAAAAAAA";
    // string ref1 = "......(.........((((.....)))).........)......................";
    // string ref2 = "................((((.....))))................................"
    string seq;
    string ref1;
    string ref2;
    if (argc >= 2){
        if (strcmp(argv[1], "test") == 0){
            while(cin >> seq && cin >> ref1 && cin >> ref2){
                cout << seq << endl;
                cout << ref1 << endl;
                cout << ref2 << endl;
                bool verbose = false;
                int dangle = 2;
                test_diff(seq, ref1, ref2, verbose, dangle);
            }
        }else if (strcmp(argv[1], "eval") == 0){
            while(cin >> seq && cin >> ref1){
                cout << seq << endl;
                cout << ref1 << endl;
                bool verbose = false;
                int dangle = 2;
                long energy = linear_eval(seq, ref1, verbose, dangle);
                printf("total energy: %.2f\n", energy/-100.0);
            }
        }else if (strcmp(argv[1], "critical") == 0){
            while(cin >> ref1 && cin >> ref2){
                cout << ref1 << endl;
                cout << ref2 << endl;
                bool verbose = false;
                find_critical(ref1, ref2, verbose);
            }
        }else if (strcmp(argv[1], "diff") == 0){
            while(cin >> seq && cin >> ref1 && cin >> ref2){
                int n;
                cin >> n;

                // cr_loops[i] = {is ref1, loop type, indices...}
                vector<vector<int>> cr_loops;
                for (int i = 0; i < n; i++) {
                    vector<int> loop (2, 0);
                    cin >> loop[0] >> loop[1];
                    int num_indices = (loop[1] == interior) ? 8 : 4;
                    int t;
                    for (int j = 0; j < num_indices; j++) {
                        cin >> t;
                        loop.push_back(t);
                    }
                    cr_loops.push_back(loop);
                }

                cout << seq << endl;
                cout << ref1 << endl;
                cout << ref2 << endl;

                bool verbose = false;
                int dangle = 2;
                long energy = diff_eval(seq, cr_loops, verbose, dangle);
                printf("energy diff: %.2f\n", energy/-100.0);
            }
        }
    }
}
