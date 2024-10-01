/**
 * @file uniq.cpp
 * @brief This script processes RNA tree structures from JSON data, deduplicates them, and prints unique tree structures.
 *
 * This script reads RNA tree structures from a file, constructs tree representations, generates all possible rotations
 * of these trees, and identifies unique tree structures. The unique trees and their counts are then printed to the console.
 *
 * The script uses the nlohmann::json library for JSON parsing and a custom Node class to represent the tree structure.
 *
 * Dependencies:
 * - nlohmann::json library
 * - utils.h (for readLinesFromFile function)
 *
 * Classes:
 * - Node: Represents a node in the RNA tree structure.
 *
 * Functions:
 * - dedup: Reads RNA tree structures from a file, deduplicates them, and prints unique trees.
 * - main: Entry point of the script. Reads the file path from command line arguments and calls the dedup function.
 *
 * Usage:
 * Compile the script and run it with the path to the input file as an argument.
 * Example: ./uniq path/to/input/file
 */
#include <iostream>
#include <fstream>

#include "json.hpp"

#include "utils.h"

using json = nlohmann::json;
using namespace std;



// Node class definition to represent the tree structure
class Node {
public:
    string type;
    vector<int> unpaired_bases;
    vector<Node*> children;
    Node* parent;
    int child_id;

    // Constructor to initialize the Node from JSON data
    Node(json jsontree = {}, Node* parent = nullptr, int child_id = -1) {
        this->parent = parent;
        this->child_id = child_id;
        if (jsontree.is_null()) {
            type = "p";  // Default type for null
        } else if (!jsontree.empty()) {
            if (jsontree.contains("root")) {  // Root node
                type = "p";
                children.push_back(new Node(jsontree["root"], this, 0));
            } else {  // Non-root node
                type = jsontree["type"];
                unpaired_bases = jsontree.contains("loops") ? jsontree["loops"].get<vector<int>>() : vector<int>();
                if (jsontree.contains("children")) {
                    for (size_t i = 0; i < jsontree["children"].size(); ++i) {
                        Node* child = new Node(jsontree["children"][i], this, i);
                        children.push_back(child);
                    }
                }else{
                    // if self.type in ["H", "M", "I", "S", "B", "E"]:
                    if(type == "H"){
                        // pass
                    }else if(type == "M"){
                        for (size_t i = 0; i < jsontree["loops"].size()-1; ++i) {
                            Node* child_null = new Node({}, this, i);
                            children.push_back(child_null);
                        }
                    }else if(type == "I" || type == "S" || type == "B"){
                        Node* child_null = new Node({}, this, 0);
                        children.push_back(child_null);
                    }else if(type == "E"){
                        assert(false);
                    }
                }
            }
        }
    }


    Node* makeTree(int child_id){
        json json_empty = {};
        Node* tree = new Node(json_empty, nullptr, -1);
        tree->type = this->type;
        tree->parent = this;
        tree->unpaired_bases = this->unpaired_bases;

        if(this->parent != nullptr){
            Node* parent_as_child = this->parent->makeTree(this->child_id);
            // tree.children = this->children[child_id+1] + [parent_as_child] + this->children[:child_id];
            // Slicing: self.children[child_id+1:] and self.children[:child_id]
            std::vector<Node*> children_behind;
            if(this->children.size() > child_id+1){
                children_behind = std::vector<Node*>(this->children.begin() + child_id + 1, this->children.end());
            }
            std::vector<Node*> children_ahead;
            if(child_id > 0 && this->children.size() >= child_id){
                children_ahead = std::vector<Node*>(this->children.begin(), this->children.begin() + child_id);
            }
            // Concatenating the three parts
            tree->children = children_behind;  // First part from child_id+1 to the end
            tree->children.push_back(parent_as_child);  // Insert the parent_as_child in between
            tree->children.insert(tree->children.end(), children_ahead.begin(), children_ahead.end());  // Append the second slice
            // tree.unpaired_bases = this->unpaired_bases[child_id+1:] + this->unpaired_bases[:child_id+1];
            // Slicing: self.unpaired_bases[child_id+1:] and self.unpaired_bases[:child_id+1]
            vector<int> part1;
            if(this->unpaired_bases.size() > child_id+1){
                part1 = vector<int>(this->unpaired_bases.begin() + child_id + 1, this->unpaired_bases.end());
            }
            vector<int> part2;
            if(child_id >= 0 && this->unpaired_bases.size() >= child_id+1){
                part2 = vector<int>(this->unpaired_bases.begin(), this->unpaired_bases.begin() + child_id + 1);
            }
            // Concatenating the two parts
            tree->unpaired_bases = part1;  // First part of the slice
            tree->unpaired_bases.insert(tree->unpaired_bases.end(), part2.begin(), part2.end());  // Append the second part
        }
        return tree;
    }


    vector<Node*> rotated(int dep=0){
        vector<Node*> rotated_trees;
        if(dep == 0 && this->type == "E"){
            return rotated_trees;
        }
        if(dep > 0 && this->type == "p"){
            Node* rote = this->makeTree(this->child_id);
            rotated_trees.push_back(rote);
        }
        for(auto& child : this->children){
            vector<Node*> child_rotated = child->rotated(dep+1);
            rotated_trees.insert(rotated_trees.end(), child_rotated.begin(), child_rotated.end());
        }
        return rotated_trees;
    }


    // Function to print the tree structure 
    string toString() {
        string result = type + " [";
        for (size_t i = 0; i < unpaired_bases.size(); ++i) {
            if (i > 0) {
                result += ", ";
            }
            result += to_string(unpaired_bases[i]);
        }
        result += "] (";
        for (size_t i = 0; i < children.size(); ++i) {
            if (i > 0) {
                result += ", ";
            }
            result += children[i]->toString();
        }
        result += ")";
        return result;
    }
};


void dedup(string path){
    vector<std::string> lines = readLinesFromFile(path);

    map<string, vector<json>> uniqs;
    set<string> ids;

    for (auto& line : lines) {
        json js = json::parse(line);
        ids.insert(js["id"]); 
        json motif = js["motif"];
        Node* tree = new Node(motif, nullptr, -1);
        // all_rotations = set([str(tree)])
        vector<Node*> all_rotations = tree->rotated(0);
        all_rotations.push_back(tree);
        bool flag = false;
        for (auto& rotated_tree : all_rotations) {
            string tree_str = rotated_tree->toString();
            // all_rotations.add(tree_str);
            if (uniqs.find(tree_str) != uniqs.end()) {
                uniqs[tree_str].push_back(js);
                flag = true;
            }
        }
        if (!flag) {
            uniqs[tree->toString()] = vector<json>{js};
        }
    }

    for(const auto&pair : uniqs){
        std::cout<<pair.first<<endl;
    }
    std::cout<< "uniq trees count: " << uniqs.size() << endl;
    std::cout<< "uniq   ids count: " << ids.size() << endl;
}


int main(int argc, char* argv[]) {

    if (argc > 1) {
        std::string path = argv[1];
        std::cout << "path: " << path << std::endl;
        dedup(path);
    } else {
        std::cout << "No argument provided." << std::endl;
    }

}