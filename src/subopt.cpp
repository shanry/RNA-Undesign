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
// #include "comps.h"
#include "eval.h"
#include <thread>

// global variables for command line arguments
double energy_delta = 5.0;
int beamsize = 0; // 0 means no pruning
bool sharpturn = false;
bool verbose = false;
int dangle = 2;

const int max_pair_diff = 3;

struct RNAStructure {
    std::string sequence;
    std::string structure;
    double energy;
};

// function to determin whether the base pairs of a structure y1 is a subset of another structure y2
bool is_subset_by_pairs(const std::string& y1, const std::string& y2) {
    std::unordered_map<int, int> subPairs = pairs_match(y1);
    std::unordered_map<int, int> fullPairs = pairs_match(y2);
    for (const auto& p : subPairs) {
        if (fullPairs.find(p.first) == fullPairs.end() || fullPairs[p.first] != p.second) {
            return false;
        }
    }
    return true;
}

// function to remove the base pairs of y1 not in y2 to make new y1 a subset of y2
std::string make_subset_by_pairs(const std::string& y1, const std::string& y2) {
    std::unordered_map<int, int> subPairs = pairs_match(y1);
    std::unordered_map<int, int> fullPairs = pairs_match(y2);
    std::string new_structure = y1;
    // remove pairs in y1 that are not in y2
    for (const auto& p : subPairs) {
        if (fullPairs.find(p.first) == fullPairs.end() || fullPairs[p.first] != p.second) {
            new_structure[p.first] = '.';
            new_structure[p.second] = '.';
        }
    }
    return new_structure;
}

void read_csv_and_fold(std::string path){
    auto df = read_csv(path.c_str());
    std::vector<std::string> targets;
    std::vector<std::string> rivals;
    for (int i = 1; i < df.size(); i++){
        auto row = df[i];
        std::string target = row[1];
        std::string sequence = row[2];
        std::string constr = std::string(target.length(), '?'); // unconstrained
        // double energy = fold_sequence(sequence, target);
        std::cout << "Processing " << i << "/" << df.size() - 1 << std::endl;
        // std::cout << "Sequence: " << sequence << std::endl;
        // std::cout << "Target: " << target << std::endl;
        // std::cout << "Constraint: " << constr << std::endl;
        // call subopt_fold
        auto structures_suboptimal = subopt_fold(sequence, constr, energy_delta, beamsize, sharpturn, verbose, dangle);
        std::cout << "Found " << structures_suboptimal.size() << " suboptimal structures within " << energy_delta << " kcal/mol:" << std::endl;
        std::string y_mfe = structures_suboptimal[0];
        bool is_subset = false;
        if (structures_suboptimal.size() == 1) {
            // calculate structural distance
            // std::cout << "MFE structure: " << y_mfe << std::endl;
            int dist = struct_dist(y_mfe, target);
            std::cout << "Structural distance: " << dist << std::endl;
            is_subset = is_subset_by_pairs(y_mfe, target);
        }else{ // find the minimum distance structure
            int min_dist = target.length();
            for (std::string y: structures_suboptimal){
                int dist = struct_dist(y, target);
                if (dist != 0 && dist < min_dist){
                    min_dist = dist;
                    y_mfe = y;
                    is_subset = is_subset_by_pairs(y, target);
                }
            }
            std::cout << "Structural distance: " << min_dist << std::endl;
        }
        // whether y_mfe is subset of target
        if (is_subset){
            std::cout << "subset? YES" << std::endl;
        }else{
            std::cout << "subset? NO" << std::endl;
            // remove the pairs in y_mfe from target to form a new structure such that the new structure is a subset of target
            std::string new_structure = make_subset_by_pairs(y_mfe, target);
            int dist = struct_dist(new_structure, target);
            std::cout << "After removing  pairs, new distance: " << dist << std::endl;
            y_mfe = new_structure;
        }
        assert (is_subset_by_pairs(y_mfe, target));
        int dist = struct_dist(y_mfe, target);
        if (dist > max_pair_diff * 2){
            // add more pairs from target to y_mfe to reduce distance to at most max_pair_diff * 2
            std::unordered_map<int, int> targetPairs = pairs_match(target);
            std::unordered_map<int, int> mfePairs = pairs_match(y_mfe);
            auto pairs_diff = std::vector<std::tuple<int, int>>(); // pairs in target but not in y_mfe
            for (const auto& p : targetPairs) {
                if (mfePairs.find(p.first) == mfePairs.end() || mfePairs[p.first] != p.second) {
                    if (p.first < p.second) {
                        pairs_diff.push_back(std::make_tuple(p.first, p.second));   
                    }
                }
            }
            // add pairs_diff.size() - max_pair_diff pairs to make distance <= max_pair_diff * 2
            int pairs_to_add = pairs_diff.size() - max_pair_diff;
            for (int j = 0; j < pairs_to_add; j++){
                int i1 = std::get<0>(pairs_diff[j]);;
                int i2 = std::get<1>(pairs_diff[j]);;
                y_mfe[i1] = '(';
                y_mfe[i2] = ')';
            }
        }
        dist = struct_dist(y_mfe, target);
        assert (dist <= max_pair_diff * 2);
        std::cout << "Final structural distance: " << dist << std::endl;
        targets.push_back(target);
        rivals.push_back(y_mfe);
        // for (auto& struct : structures_suboptimal){
        //     std::cout << struct << std::endl;
        // }
        std::cout << "----------------------------------------" << std::endl;
    }
    // write targets and rivals to txt file, separated by comma
    std::ofstream outfile("target_rival.txt");
    for (int i = 0; i < targets.size(); i++){
        outfile << targets[i] << "," << rivals[i] << std::endl;
    }
    outfile.close();
    std::cout << "Wrote target and rival structures to target_rival.txt" << std::endl;
}

int main(int argc, char** argv) {
    cxxopts::Options options("RNA Suboptimal Structure Generator", "Generates suboptimal RNA secondary structures within a specified energy range.");

    options.add_options()
        ("p,path", "Input file containing target structures and designed RNA sequences in csv format", cxxopts::value<std::string>())
        ("e,energy_delta", "Energy range for suboptimal structures", cxxopts::value<double>()->default_value("0.0"))
        ("b,beamsize", "Beam size for beam search algorithm", cxxopts::value<int>()->default_value("0"))
        ("s,sharpturn", "Enable sharp turn heuristic", cxxopts::value<bool>()->default_value("false"))
        ("d,dangle", "Dangle model to use (0, 1, 2)", cxxopts::value<int>()->default_value("2"))
        ("v,verbose", "Enable verbose output", cxxopts::value<bool>()->default_value("false"))
        ("h,help", "Print usage");

    auto result = options.parse(argc, argv);

    // print path
    if (result.count("help") || !result.count("path")) {
        std::cout << options.help() << std::endl;
        return 0;
    }

    std::string path = result["path"].as<std::string>();
    std::cout << "Input file path: " << path << std::endl;

    energy_delta = result["energy_delta"].as<double>();
    beamsize = result["beamsize"].as<int>();
    sharpturn = result["sharpturn"].as<bool>();
    verbose = result["verbose"].as<bool>();
    dangle = result["dangle"].as<int>();
    std::cout << "Energy delta: " << energy_delta << std::endl;
    std::cout << "Beam size: " << beamsize << std::endl;
    std::cout << "Sharp turn: " << sharpturn << std::endl;
    std::cout << "Verbose: " << verbose << std::endl;
    std::cout << "Dangle: " << dangle << std::endl;

    read_csv_and_fold(path);

    return 0;

}