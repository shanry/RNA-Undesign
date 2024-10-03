#include<iostream>
#include<cstdio>
#include<tuple>
#include <map>
#include<utility>
#include <cassert>
#include <cmath>
#include <algorithm>

// #include "json.hpp"
// using json = nlohmann::json;
#include "utils.h"
#include "eval.h"
#include "comps.h"

char nuc_all[] = "ACGU";
char nuc_pair_all[][3] = {"GC", "CG", "AU", "UA", "GU", "UG"};

extern std::vector<std::pair<int, int>> pairs_outside;
extern std::vector<std::pair<int, int>> pairs_inside;

TreeNode::TreeNode(int first_val, int second_val){
        first = first_val;
        second = second_val;
    }

void TreeNode::setLoop(){
        if(first == -1){
            assert(parent == NULL);
            looptype = "E";
        }else{
            // parent is the first neighbor
            neighbors.push_back(0);
            if(children.size() == 0){
                looptype = "H";
            }
            else if (children.size() == 1){
                if(first+1 == children[0]->first && second-1 == children[0]->second)
                    looptype = "S";
                else if (first+1 == children[0]->first || second-1 == children[0]->second)
                    looptype = "B";
                else if (first+1 < children[0]->first || second-1 > children[0]->second)
                    looptype = "I";
                else{
                    std::cerr << "two pairs of internal/bulge/stack are not correct!" << std::endl;
                    assert(false);
                }
            }
            else{
                looptype = "M";
            }
        }
        // children are neighbors
        for(int ic = 0; ic < children.size(); ic++)
            neighbors.push_back(ic+1);
        // lens of (unpaired) loops
        if(children.size() == 0)
            looplens.push_back(second - first - 1);
        else{
            for(int i = 0; i < children.size(); i++){
                if(i == 0)
                    looplens.push_back(children[i]->first - first - 1);
                else
                    looplens.push_back(children[i]->first - children[i-1]->second - 1);
            }
            looplens.push_back(second - children.back()->second - 1);
        }
        assert(looplens.size() == children.size()+1);
        looplen = 0;
        for(auto& seglen: looplens)
            looplen += seglen;
        json js;
        js["type"] = looptype;
        js["neighbors"] = neighbors;
        js["loops"] = looplens;
        js["child_id"] = child_id;
        jstring = js.dump();
    }

void TreeNode::printTree(){
    printf("first: %d, second: %d\n", first, second);
    for(int i = 0; i < children.size(); i++){
        printf("child[%d]: first: %d, second: %d\n", i, children[i]->first, children[i]->second);
    }
    printf("\n");
    for(int j = 0; j < children.size(); j++){
        children[j]->printTree();
    }       
}

void TreeNode::printTree(std::string ref){
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

void TreeNode::printTree(std::string& ref, std::string& seq){
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

void TreeNode::printTree(std::string& ref, std::string& seq, std::vector<std::pair<std::string, std::string>>& subrefs){
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

bool TreeNode::isTooLong(){
    if(looptype == "I" || looptype == "B")
        return looplen > SINGLE_MAX_LEN;
    else if(looptype == "M")
        return looplens[0] > MULTIPLE_FIRST_MAX_LEN;
    return false;
}

bool TreeNode::isTooShort(){
    if(looptype == "H")
        return looplen < 3; // minimum hairpin loop length
    return false;
}

LoopComplex::LoopComplex(int count, std::string y, std::string cs, int s, int e, TreeNode* n, int l, int r, std::vector<std::pair<int, int>> p_out, std::vector<std::pair<int, int>> p_in){
    count_uk = count;
    ref = y;
    constr = cs;
    start = s;
    end = e;
    node = n;
    left = l;
    right = r;
    ps_outside = p_out;
    ps_inside = p_in;

}

void LoopComplex::set_neighbors(std::vector<int> nbs){
    neighbors = nbs;
}

bool LoopComplex::hasLongLoop(){
    // printf("inside hasLongLoop\n");
    if(node->isTooLong())
        return true;
    else{
        for(int nb: neighbors){
            if (nb == 0){
                if (node->parent->isTooLong())
                    return true;
            }else if (node->children[nb-1]->isTooLong())
                return true;
        }
    }
    return false;
}

bool LoopComplex::hasShortLoop(){
    // printf("inside hasShortLoop\n");
    if(node->isTooShort())
        return true;
    else{
        for(int nb: neighbors){
            if (nb == 0){
                if (node->parent->isTooShort())
                    return true;
            }else if (node->children[nb-1]->isTooShort())
                return true;
        }
    }
    return false;
}

void LoopComplex::printLoopLens(){
    printf("self loop: %s, %d\n", node->looptype.c_str(), node->looplen);
    for(int nb: neighbors){
        printf("nb: %d\n", nb);
        if (nb == 0)
            printf("parent loop: %s, %d\n", node->parent->looptype.c_str(), node->parent->looplen);
        else
            printf("child %d loop: %s, %d\n", nb, node->children[nb-1]->looptype.c_str(), node->children[nb-1]->looplen);
    }
}

// Constructor that accepts two arguments to initialize the members
Constraint::Constraint(std::set<int>* index_set, std::vector<std::string>* seq_list){
    indices = index_set;
    seqs = seq_list;
}

Constraint::Constraint(std::set<int>* index_set, std::vector<std::string>* seq_list, std::string& ref){
    indices = index_set;
    seqs = seq_list;
    structure = ref;
}

// Constructor to initialize the Node from JSON data
Node::Node(json jsontree, Node* parent, int child_id){
    this->parent = parent;
    this->child_id = child_id;
    if (jsontree.is_null()) {
        type = "p";  // Default type for null
    } else if (!jsontree.empty()) {
        if (jsontree.contains("root")) {  // Root node
            if (jsontree["root"]["child_id"] == -1){
                type = "53"; // 5' and 3' end
                this->child_id = -1;
            }else{
                type = "p";
            }
            children.push_back(new Node(jsontree["root"], this, 0));
        } else {  // Non-root node
            type = jsontree["type"];
            unpaired_bases = jsontree.contains("loops") ? jsontree["loops"].get<std::vector<int>>() : std::vector<int>();
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

// Function to create a roated tree from a child ID
Node* Node::makeTree(int child_id){
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
        std::vector<int> part1;
        if(this->unpaired_bases.size() > child_id+1){
            part1 = std::vector<int>(this->unpaired_bases.begin() + child_id + 1, this->unpaired_bases.end());
        }
        std::vector<int> part2;
        if(child_id >= 0 && this->unpaired_bases.size() >= child_id+1){
            part2 = std::vector<int>(this->unpaired_bases.begin(), this->unpaired_bases.begin() + child_id + 1);
        }
        // Concatenating the two parts
        tree->unpaired_bases = part1;  // First part of the slice
        tree->unpaired_bases.insert(tree->unpaired_bases.end(), part2.begin(), part2.end());  // Append the second part
    }
    return tree;
}

// Function to get all rotated trees
std::vector<Node*> Node::rotated(int dep){
    std::vector<Node*> rotated_trees;
    if(dep == 0 && this->type == "53"){ // 5' and 3' end
        return rotated_trees;
    }
    if(dep > 0 && this->type == "p"){
        Node* rote = this->makeTree(this->child_id);
        rotated_trees.push_back(rote);
    }
    for(auto& child : this->children){
        std::vector<Node*> child_rotated = child->rotated(dep+1);
        rotated_trees.insert(rotated_trees.end(), child_rotated.begin(), child_rotated.end());
    }
    return rotated_trees;
}

// Function to print the tree structure 
std::string Node::toString() const{
    std::string result = type + " [";
    for (size_t i = 0; i < unpaired_bases.size(); ++i) {
        if (i > 0) {
            result += ", ";
        }
        result += std::to_string(unpaired_bases[i]);
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

// Function to print the tree structure 
std::string Node::toDotBracket() const{
    std::string result = "";
    if (this->type == "53"){ // root: 5' and 3' end
        result += "5";
        result += this->children[0]->toDotBracket();
        result += "3";
    }else if (this->type == "p"){ // leaf or root
        if(this->child_id == -1)
            result = this->children[0]->toDotBracket();
        else
            result = "(***)";
    }else if (this->type == "H"){ // hairpin
        result = "(" + std::string(this->unpaired_bases[0], '.') + ")";
    }else{ // multi-loop, external, internal, bulge, stack
        result = ""; //  + this->children[0]->toDotBracket() + ")" + this->children[1]->toDotBracket();
        if (this->type != "E")
            result += "(";
        result += std::string(this->unpaired_bases[0], '.');
        for(int i = 0; i < this->children.size(); i++){
            result += this->children[i]->toDotBracket();
            assert (i + 1 < this->unpaired_bases.size());
            result += std::string(this->unpaired_bases[i+1], '.');
        }
        if (this->type != "E")
            result += ")";
    }
    return result;
}

void Constraint::setStructure(std::string& ref){
    structure = ref;
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

void tree2Loops(TreeNode* root, std::string& ref, std::vector<LoopComplex>& lc_list){
    if (root->first == -1){
        printf("external: first: %d, second: %d; children: %d\n", root->first, root->second, root->children.size());
    }else if (root->children.size() == 0){
        printf("haiprin: first: %d, second: %d; children: %d\n", root->first, root->second, root->children.size());
    }else if(root->children.size() == 1){
        printf("internal: first: %d, second: %d; children: %d\n", root->first, root->second, root->children.size());
    }else if(root->children.size() > 1){
        printf("multi-loop: first: %d, second: %d; children: %d\n", root->first, root->second, root->children.size());
    }
    if (root->parent != NULL){  // make sure it is not the pseudo-node
        printf("parent: first: %d, second: %d; children: %d\n", root->parent->first, root->parent->second, root->parent->children.size());
        printf("   ref: %s\n", ref.c_str());
        std::string constr = removeNodeFromTree(root, ref);
        printf("constr: %s\n", constr.c_str());
        std::string cref = constr;
        std::replace(cref.begin(), cref.end(), '?', '.');
        printf("  cref: %s\n", cref.c_str());
        int count_unknown = 0;
        for(auto ch: constr){
            if (ch == '?')
                count_unknown += 1;
        }
        printf("? count: %d\n", count_unknown);
        int start = root->parent->first;
        int end = root->parent->second;
        if (root->parent->first < 0){
            start = 0;
            end = ref.length() - 1;
        }
        LoopComplex lc(count_unknown, cref, constr, start, end, root, root->first, root->second, pairs_outside, pairs_inside);
        lc_list.push_back(lc);
    }
    printf("\n");
    for(auto child: root->children){
        tree2Loops(child, ref, lc_list);
    }
    return;
}

std::string removeNodeFromTree(TreeNode* node, std::string ref){
    std::string constr(ref.length(), '?');
    int len_p = node->parent->second - node->parent->first + 1;
    if (node->parent->first >= 0){
        printf("first: %d, second: %d\n", node->parent->first, node->parent->second);
        printf(" len_p: %d\n", len_p);
        printf("   ref: %s\n", ref.substr(node->parent->first, len_p).c_str());
        constr[node->parent->first] = '(';
        constr[node->parent->second] = ')';
        printf("constr: %s\n", constr.substr(node->parent->first, len_p).c_str());
    }else{
        printf(" len_p: %d\n", ref.length());
        printf("   ref: %s\n", ref.c_str());
        printf("constr: %s\n", constr.c_str());
    }
    // for(int i = 0; i < ref.length(); i++){
    //     if(i <= node->parent->first || i >= node->parent->second)
    //         constr[i] = ref[i];
    // }
    for(auto child: node->children){
        if (child->children.size()){
            for(auto grandchild: child->children){
                for(int i = grandchild->first; i <= grandchild->second; i++){
                    constr[i] = ref[i];
                }
            }
        }
    }
    if (node->parent->first >= 0)
        printf("constr: %s\n", constr.substr(node->parent->first, len_p).c_str());
    else
        printf("constr: %s\n", constr.c_str());
    for(auto sibling: node->parent->children){
        if (sibling != node){
            for(int i = sibling->first; i <= sibling->second; i++){
                constr[i] = ref[i];
            }
        }
    }
    if (node->parent->first >= 0)
        return constr.substr(node->parent->first, len_p);
    else
        return constr;
}

void tree2TwoNeighbor(TreeNode* root, std::string& ref, std::vector<LoopComplex>& lc_list){
     printf("inside tree2TwoNeighbor\n");
    if (root->first == -1){
        printf("external: first: %d, second: %d; children: %d\n", root->first, root->second, root->children.size());
    }else if (root->children.size() == 0){
        printf("haiprin: first: %d, second: %d; children: %d\n", root->first, root->second, root->children.size());
    }else if(root->children.size() == 1){
        printf("internal: first: %d, second: %d; children: %d\n", root->first, root->second, root->children.size());
    }else if(root->children.size() > 1){
        printf("multi-loop: first: %d, second: %d; children: %d\n", root->first, root->second, root->children.size());
    }
    printf("before if\n");
    if (root != NULL && root->children.size() > 0){
        printf("in if\n");
        printf("parent: first: %d, second: %d; children: %d\n", root->first, root->second, root->children.size());
        printf("   ref: %s\n", ref.c_str());
        int count_neighbor = root->children.size() + 1;
        int start_neighbor = 0;
        if (root->parent == NULL){
            start_neighbor++; // from the first child
        }
        for(int i = start_neighbor; i < count_neighbor; i++){
            for(int j = i + 1; j < count_neighbor; j++){
                pairs_outside.clear();
                pairs_inside.clear();
                std::string constr = removeTwoNeighbors(root, ref, i, j);
                printf("constr: %s\n", constr.c_str());
                if(i==0){
                    pairs_inside.push_back(std::make_pair(root->first, root->second));
                }else{
                    pairs_inside.push_back(std::make_pair(root->children[i-1]->first, root->children[i-1]->second));
                }
                pairs_inside.push_back(std::make_pair(root->children[j-1]->first, root->children[j-1]->second));
                std::string cref = constr;
                std::replace(cref.begin(), cref.end(), '?', '.');
                printf("  cref: %s\n", cref.c_str());
                int count_unknown = 0;
                for(auto ch: constr){
                    if (ch == '?')
                        count_unknown += 1;
                }
                printf("? count: %d\n", count_unknown);
                int start, end;
                if(i == 0){
                    start = root->parent->first;
                    end = root->parent->second;
                    if (root->parent->first < 0){
                        start = 0;
                        end = ref.length() - 1;
                    }
                }else{
                    start = root->first;
                    end = root->second;
                    if (root->first < 0){
                        start = 0;
                        end = ref.length() - 1;
                    }
                }
                LoopComplex lc(count_unknown, cref, constr, start, end, root, root->first, root->second, pairs_outside, pairs_inside);
                lc.set_neighbors({i, j});
                lc_list.push_back(lc);
            }
        }
    }
    printf("\n");
    for(auto child: root->children){
        tree2TwoNeighbor(child, ref, lc_list);
    }
    return;
}

void tree2ThreeNeighbor(TreeNode* root, std::string& ref, std::vector<LoopComplex>& lc_list){
     printf("inside tree2ThreeNeighbor\n");
    if (root->first == -1){
        printf("external: first: %d, second: %d; children: %d\n", root->first, root->second, root->children.size());
    }else if (root->children.size() == 0){
        printf("haiprin: first: %d, second: %d; children: %d\n", root->first, root->second, root->children.size());
    }else if(root->children.size() == 1){
        printf("internal: first: %d, second: %d; children: %d\n", root->first, root->second, root->children.size());
    }else if(root->children.size() > 1){
        printf("multi-loop: first: %d, second: %d; children: %d\n", root->first, root->second, root->children.size());
    }
    // printf("before if\n");
    if (root != NULL && root->children.size() > 1){
        printf("in if\n");
        printf("parent: first: %d, second: %d; children: %d\n", root->first, root->second, root->children.size());
        printf("   ref: %s\n", ref.c_str());
        int count_neighbor = root->children.size() + 1;
        int start_neighbor = 0;
        if (root->parent == NULL){
            start_neighbor++; // from the first child
            count_neighbor--;
        }
        std::vector<std::vector<int>> powset3 = PowerSet3(start_neighbor, count_neighbor);
        for(auto ps: powset3){
            // for(int j = i + 1; j < count_neighbor; j++){
                pairs_outside.clear();
                pairs_inside.clear();
                std::string constr = removeThreeNeighbors(root, ref, ps);
                printf("constr: %s\n", constr.c_str());
                for(int j: ps){
                    if(j==0)
                        pairs_inside.push_back(std::make_pair(root->first, root->second));
                    else
                        pairs_inside.push_back(std::make_pair(root->children[j-1]->first, root->children[j-1]->second));
                }
                std::string cref = constr;
                std::replace(cref.begin(), cref.end(), '?', '.');
                printf("  cref: %s\n", cref.c_str());
                int count_unknown = 0;
                for(auto ch: constr){
                    if (ch == '?')
                        count_unknown += 1;
                }
                printf("? count: %d\n", count_unknown);
                int start, end;
                if(ps[0] == 0){
                    start = root->parent->first;
                    end = root->parent->second;
                    if (root->parent->first < 0){
                        start = 0;
                        end = ref.length() - 1;
                    }
                }else{
                    start = root->first;
                    end = root->second;
                    if (root->first < 0){
                        start = 0;
                        end = ref.length() - 1;
                    }
                }
                LoopComplex lc(count_unknown, cref, constr, start, end, root, root->first, root->second, pairs_outside, pairs_inside);
                lc.set_neighbors(ps);
                lc_list.push_back(lc);
            // }
        }
    }
    printf("\n");
    for(auto child: root->children){
        tree2ThreeNeighbor(child, ref, lc_list);
    }
    return;
}

void tree2Edges(TreeNode* root, std::string& ref, std::vector<LoopComplex>& lc_list){
    if (root->first == -1){
        printf("external: first: %d, second: %d; children: %d\n", root->first, root->second, root->children.size());
    }else if (root->children.size() == 0){
        printf("haiprin: first: %d, second: %d; children: %d\n", root->first, root->second, root->children.size());
    }else if(root->children.size() == 1){
        printf("internal: first: %d, second: %d; children: %d\n", root->first, root->second, root->children.size());
    }else if(root->children.size() > 1){
        printf("multi-loop: first: %d, second: %d; children: %d\n", root->first, root->second, root->children.size());
    }
    if (root->parent != NULL){
        printf("parent: first: %d, second: %d; children: %d\n", root->parent->first, root->parent->second, root->parent->children.size());
        printf("   ref: %s\n", ref.c_str());
        pairs_inside.clear();
        pairs_outside.clear();
        std::string constr = removeEdgeFromTree(root, ref);
        printf("constr: %s\n", constr.c_str());
        pairs_inside.push_back(std::make_pair(root->first, root->second));
        std::string cref = constr;
        std::replace(cref.begin(), cref.end(), '?', '.');
        printf("  cref: %s\n", cref.c_str());
        int count_unknown = 0;
        for(auto ch: constr){
            if (ch == '?')
                count_unknown += 1;
        }
        printf("? count: %d\n", count_unknown);
        int start = root->parent->first;
        int end = root->parent->second;
        if (root->parent->first < 0){
            start = 0;
            end = ref.length() - 1;
        }
        assert(pairs_inside.size() == 1);
        LoopComplex lc(count_unknown, cref, constr, start, end, root, root->first, root->second, pairs_outside, pairs_inside);
        lc.set_neighbors({0});
        lc_list.push_back(lc);
    }
    printf("\n");
    for(auto child: root->children){
        tree2Edges(child, ref, lc_list);
    }
    return;
}

std::string removeEdgeFromTree(TreeNode* node, std::string ref){
    std::string constr(ref.length(), '?');
    int len_p = node->parent->second - node->parent->first + 1;
    if (node->parent->first >= 0){
        printf("first: %d, second: %d\n", node->parent->first, node->parent->second);
        printf(" len_p: %d\n", len_p);
        printf("   ref: %s\n", ref.substr(node->parent->first, len_p).c_str());
        constr[node->parent->first] = '(';
        constr[node->parent->second] = ')';
        printf("constr: %s\n", constr.substr(node->parent->first, len_p).c_str());
        pairs_outside.push_back(std::make_pair(node->parent->first, node->parent->second));
    }else{
        printf(" len_p: %d\n", ref.length());
        printf("   ref: %s\n", ref.c_str());
        printf("constr: %s\n", constr.c_str());
        pairs_outside.push_back(std::make_pair(-1, -1));
    }
    for(auto sibling: node->parent->children){
        if (sibling != node){
            for(int i = sibling->first; i <= sibling->second; i++){
                constr[i] = ref[i];
            }
            pairs_outside.push_back(std::make_pair(sibling->first, sibling->second));
        }
    }
    if (node->parent->first >= 0)
        printf("constr: %s\n", constr.substr(node->parent->first, len_p).c_str());
    else
        printf("constr: %s\n", constr.c_str());
    for(auto child: node->children){
        for(int i = child->first; i <= child->second; i++){
            constr[i] = ref[i];
        }
        pairs_outside.push_back(std::make_pair(child->first, child->second));
    }
    if (node->parent->first >= 0)
        return constr.substr(node->parent->first, len_p);
    else
        return constr;
}

void tree2MLoops(TreeNode* root, std::string& ref, std::vector<LoopComplex>& lc_list){
    printf("inside tree2MLoops\n");
    if (root->first == -1){
        printf("external: first: %d, second: %d; children: %d\n", root->first, root->second, root->children.size());
    }else if (root->children.size() == 0){
        printf("haiprin: first: %d, second: %d; children: %d\n", root->first, root->second, root->children.size());
    }else if(root->children.size() == 1){
        printf("internal: first: %d, second: %d; children: %d\n", root->first, root->second, root->children.size());
    }else if(root->children.size() > 1){
        printf("multi-loop: first: %d, second: %d; children: %d\n", root->first, root->second, root->children.size());
    }
    printf("before if\n");
    if (root != NULL && root->children.size() > 1){
        printf("in if\n");
        printf("parent: first: %d, second: %d; children: %d\n", root->first, root->second, root->children.size());
        printf("   ref: %s\n", ref.c_str());
        std::string constr = removeMNodeFromTree(root, ref);
        printf("constr: %s\n", constr.c_str());
        std::string cref = constr;
        std::replace(cref.begin(), cref.end(), '?', '.');
        printf("  cref: %s\n", cref.c_str());
        int count_unknown = 0;
        for(auto ch: constr){
            if (ch == '?')
                count_unknown += 1;
        }
        printf("? count: %d\n", count_unknown);
        int start = root->first;
        int end = root->second;
        if (root->first < 0){
            start = 0;
            end = ref.length() - 1;
        }
        LoopComplex lc(count_unknown, cref, constr, start, end, root, root->first, root->second, pairs_outside, pairs_inside);
        lc_list.push_back(lc);
    }
    printf("\n");
    for(auto child: root->children){
        tree2MLoops(child, ref, lc_list);
    }
    return;
}

std::string removeMNodeFromTree(TreeNode* node, std::string ref){
    std::string constr(ref.length(), '?');
    int len_p = node->second - node->first + 1;
    if (node->first >= 0){
        printf("first: %d, second: %d\n", node->first, node->second);
        printf(" len_p: %d\n", len_p);
        printf("   ref: %s\n", ref.substr(node->first, len_p).c_str());
        constr[node->first] = '(';
        constr[node->second] = ')';
        printf("constr: %s\n", constr.substr(node->parent->first, len_p).c_str());
    }else{
        printf(" len_p: %d\n", ref.length());
        printf("   ref: %s\n", ref.c_str());
        printf("constr: %s\n", constr.c_str());
    }
    // for(int i = 0; i < ref.length(); i++){
    //     if(i <= node->parent->first || i >= node->parent->second)
    //         constr[i] = ref[i];
    // }
    for (auto sibling: node->children){
        for(auto child: sibling->children){
            if (child->children.size()){
                for(auto grandchild: child->children){
                    for(int i = grandchild->first; i <= grandchild->second; i++){
                        constr[i] = ref[i];
                    }
                }
            }
        }
    }
    if (node->first >= 0)
        printf("constr: %s\n", constr.substr(node->first, len_p).c_str());
    else
        printf("constr: %s\n", constr.c_str());
    // for(auto sibling: node->parent->children){
    //     if (sibling != node){
    //         for(int i = sibling->first; i <= sibling->second; i++){
    //             constr[i] = ref[i];
    //         }
    //     }
    // }
    if (node->first >= 0)
        return constr.substr(node->first, len_p);
    else
        return constr;
}

std::string removeTwoNeighbors(TreeNode* node, std::string ref, int n1, int n2){
    std::string constr(ref.length(), '?');
    // the parent of node and a child of the current node
    printf("n1: %d, n2: %d\n", n1, n2);
    if (n1 == 0){
        int len_p = node->parent->second - node->parent->first + 1;
        if (node->parent->first >= 0){
            printf("first: %d, second: %d\n", node->parent->first, node->parent->second);
            printf(" len_p: %d\n", len_p);
            printf("   ref: %s\n", ref.substr(node->parent->first, len_p).c_str());
            constr[node->parent->first] = '(';
            constr[node->parent->second] = ')';
            printf("constr: %s\n", constr.substr(node->parent->first, len_p).c_str());
            pairs_outside.push_back(std::make_pair(node->parent->first, node->parent->second));
        }else{
            printf(" len_p: %d\n", ref.length());
            printf("   ref: %s\n", ref.c_str());
            printf("constr: %s\n", constr.c_str());
            pairs_outside.push_back(std::make_pair(-1, -1));
        }
        for(auto sibling: node->parent->children){
            if (sibling != node){
                for(int i = sibling->first; i <= sibling->second; i++){
                    constr[i] = ref[i];
                }
                pairs_outside.push_back(std::make_pair(sibling->first, sibling->second));
            }
        }
        if (node->parent->first >= 0)
            printf("constr: %s\n", constr.substr(node->parent->first, len_p).c_str());
        else
            printf("constr: %s\n", constr.c_str());

        auto child = node->children[n2-1]; // parent as neibor 0 causing decreament of children's indices by 1
        for(auto sibling_child: node->children){
            if (sibling_child != child){
                for(int i = sibling_child->first; i <= sibling_child->second; i++){
                    constr[i] = ref[i];
                }
                pairs_outside.push_back(std::make_pair(sibling_child->first, sibling_child->second));
            }
        }
        if (child->children.size()){
            for(auto grandchild: child->children){
                for(int i = grandchild->first; i <= grandchild->second; i++){
                    constr[i] = ref[i];
                }
                pairs_outside.push_back(std::make_pair(grandchild->first, grandchild->second));
            }
        }

        if (node->parent->first >= 0)
            return constr.substr(node->parent->first, len_p);
        else
            return constr;
    }
    // two childrent of the current node
    if (n1 > 0){
        int len_p = node->second - node->first + 1;
        if (node->first >= 0){
            printf("first: %d, second: %d\n", node->first, node->second);
            printf(" len_p: %d\n", len_p);
            printf("   ref: %s\n", ref.substr(node->first, len_p).c_str());
            constr[node->first] = '(';
            constr[node->second] = ')';
            printf("constr: %s\n", constr.substr(node->first, len_p).c_str());
            pairs_outside.push_back(std::make_pair(node->first, node->second));
        }else{
            printf(" len_p: %d\n", ref.length());
            printf("   ref: %s\n", ref.c_str());
            printf("constr: %s\n", constr.c_str());
            pairs_outside.push_back(std::make_pair(-1, -1));
        }
        std::vector<TreeNode*> two_children;
        two_children.push_back(node->children[n1-1]); // parent as neibor 0 causing decreament of children's indices by 1
        two_children.push_back(node->children[n2-1]);
        for(auto sibling_child: node->children){
            if ( (sibling_child != two_children[0])&& (sibling_child != two_children[1]) ){
                for(int i = sibling_child->first; i <= sibling_child->second; i++){
                    constr[i] = ref[i];
                }
                pairs_outside.push_back(std::make_pair(sibling_child->first, sibling_child->second));
            }
        }
        for (auto child: two_children){
                if (child->children.size()){
                    for(auto grandchild: child->children){
                        for(int i = grandchild->first; i <= grandchild->second; i++){
                            constr[i] = ref[i];
                        }
                    pairs_outside.push_back(std::make_pair(grandchild->first, grandchild->second));
                    }
                }
        }
        if (node->first >= 0)
            printf("constr: %s\n", constr.substr(node->first, len_p).c_str());
        else
            printf("constr: %s\n", constr.c_str());
        if (node->first >= 0)
            return constr.substr(node->first, len_p);
        else
            return constr;
    }
    assert(false);
    return constr;
}

std::string removeThreeNeighbors(TreeNode* node, std::string ref, std::vector<int>& powset){
    std::string loopcode = std::to_string(node->first) + "," + std::to_string(node->second) + ":";
    for(auto nbid: powset)
        loopcode += std::to_string(nbid) + ";";
    std::cout<<"loopcode: "<<loopcode<<std::endl;
    std::string constr(ref.length(), '?');
    // the parent of node and a child of the current node
    if (powset[0] == 0){
        int len_p = node->parent->second - node->parent->first + 1;
        if (node->parent->first >= 0){
            printf("first: %d, second: %d\n", node->parent->first, node->parent->second);
            printf(" len_p: %d\n", len_p);
            printf("   ref: %s\n", ref.substr(node->parent->first, len_p).c_str());
            constr[node->parent->first] = '(';
            constr[node->parent->second] = ')';
            printf("constr: %s\n", constr.substr(node->parent->first, len_p).c_str());
            pairs_outside.push_back(std::make_pair(node->parent->first, node->parent->second));
        }else{
            printf(" len_p: %d\n", ref.length());
            printf("   ref: %s\n", ref.c_str());
            printf("constr: %s\n", constr.c_str());
            pairs_outside.push_back(std::make_pair(-1, -1));
        }

        for(auto sibling: node->parent->children){
            if (sibling != node){
                for(int i = sibling->first; i <= sibling->second; i++){
                    constr[i] = ref[i];
                }
                pairs_outside.push_back(std::make_pair(sibling->first, sibling->second));
            }
        }

        if (node->parent->first >= 0)
            printf("constr: %s\n", constr.substr(node->parent->first, len_p).c_str());
        else
            printf("constr: %s\n", constr.c_str());

        std::set<int> setChildren(powset.begin(), powset.end());
        for(int chid = 0; chid < node->children.size(); chid++ ){
            if (setChildren.find(chid+1)==setChildren.end()){
                auto sibling_child = node->children[chid];
                for(int i = sibling_child->first; i <= sibling_child->second; i++){
                    constr[i] = ref[i];
                }
                pairs_outside.push_back(std::make_pair(sibling_child->first, sibling_child->second));
            }
        }

        for(int id_child: powset){
            if(id_child > 0){
                auto child = node->children[id_child-1]; // parent as neibor 0 causing decreament of children's indices by 1
                if (child->children.size()){
                    for(auto grandchild: child->children){
                        for(int i = grandchild->first; i <= grandchild->second; i++){
                            constr[i] = ref[i];
                        }
                        pairs_outside.push_back(std::make_pair(grandchild->first, grandchild->second));
                        std::cout<<"add grandchild to pairs_outside: "<<grandchild->first<<" "<<grandchild->second<<std::endl; 
                    }
                }
            }
        }
        if (node->parent->first >= 0)
            return constr.substr(node->parent->first, len_p);
        else
            return constr;
    }
    // two childrent of the current node
    if (powset[0] > 0){
        int len_p = node->second - node->first + 1;
        if (node->first >= 0){
            printf("first: %d, second: %d\n", node->first, node->second);
            printf(" len_p: %d\n", len_p);
            printf("   ref: %s\n", ref.substr(node->first, len_p).c_str());
            constr[node->first] = '(';
            constr[node->second] = ')';
            printf("constr: %s\n", constr.substr(node->first, len_p).c_str());
            pairs_outside.push_back(std::make_pair(node->first, node->second));
        }else{
            printf(" len_p: %d\n", ref.length());
            printf("   ref: %s\n", ref.c_str());
            printf("constr: %s\n", constr.c_str());
            pairs_outside.push_back(std::make_pair(-1, -1));
        }
        printf("mark1\n");
        std::set<int> setChildren(powset.begin(), powset.end());
        for(int chid = 0; chid < node->children.size(); chid++ ){
            if (setChildren.find(chid+1)==setChildren.end()){
                auto sibling_child = node->children[chid];
                for(int i = sibling_child->first; i <= sibling_child->second; i++){
                    constr[i] = ref[i];
                }
                pairs_outside.push_back(std::make_pair(sibling_child->first, sibling_child->second));
            }
        }
        printf("mark2\n");
        std::vector<TreeNode*> pow_children;
        std::cout<<"powset.size: "<<powset.size()<<std::endl;
        for(int id_child: powset){
            std::cout<<"id_child: "<<id_child<<std::endl;
            assert(id_child > 0);
            pow_children.push_back(node->children[id_child - 1]); // parent as neibor 0 causing decreament of children's indices by 1
        }
        std::cout<<"pow_children.size: "<<pow_children.size()<<std::endl;
        for (auto child: pow_children){
            std::cout<<"pow_children: "<<child->children.size()<<"\t"<<child->first<<"\t"<<child->second<<std::endl;
            if (child->children.size() > 0){
                for(auto grandchild: child->children){
                    for(int i = grandchild->first; i <= grandchild->second; i++){
                        constr[i] = ref[i];
                    }
                    pairs_outside.push_back(std::make_pair(grandchild->first, grandchild->second));
                    std::cout<<"add grandchild to pairs_outside: "<<grandchild->first<<" "<<grandchild->second<<std::endl; 
                }
            }
        }
        if (node->first >= 0)
            printf("constr: %s\n", constr.substr(node->first, len_p).c_str());
        else
            printf("constr: %s\n", constr.c_str());
        if (node->first >= 0)
            return constr.substr(node->first, len_p);
        else
            return constr;
    }
    assert(false);
    return constr;
}

int max_hairpin(TreeNode* root){
    int maxlen = 0;
    if (root->first!=-1&&root->children.size()==0){
        maxlen = root->second - root->first - 1;
    }
    for(TreeNode* child: root->children){
        int maxlen_child = max_single(child);
        maxlen = std::max(maxlen, maxlen_child);
    }
    return maxlen;
}

int max_single(TreeNode* root){
    int maxlen = 0;
    if (root->first!=-1&&root->children.size()==1){
        maxlen = root->looplen;
    }
    for(TreeNode* child: root->children){
        int maxlen_child = max_single(child);
        maxlen = std::max(maxlen, maxlen_child);
    }
    return maxlen;
}

int max_multi(TreeNode* root){
    int maxlen = 0;
    if (root->first!=-1&&root->children.size()>1){
        maxlen = root->looplens[0];
    }
    for(TreeNode* child: root->children){
        int maxlen_child = max_multi(child);
        maxlen = std::max(maxlen, maxlen_child);
    }
    return maxlen;
}

std::string ml_degree(TreeNode* root){
    std::string ml_str = "";
    if (root->children.size() > 1){
        ml_str += std::to_string(root->children.size()) + ":" + std::to_string(root->first) + "," + std::to_string(root->second) + ";";
         for(TreeNode* child: root->children){
            ml_str += std::to_string(child->first) + "," + std::to_string(child->second) + ";";
         }
    }
    for(TreeNode* child: root->children){
        ml_str += ml_degree(child);
    } 
    return ml_str;
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

std::string LoopComplex::jsmotif(std::string id){
    json js;
    js["id"] = id;
    if(neighbors.size() == 1){
        assert(neighbors[0]==0);
        json parent_j = json::parse(node->parent->jstring);
        json children_j;
        for(int i = 0; i < node->parent->children.size(); i++){
            if(i+1 == node->child_id)
                children_j.push_back(json::parse(node->jstring));
            else
                children_j.push_back(nullptr);
        }
        parent_j["children"] = children_j;
        js["root"] = parent_j;
    }
    else if(neighbors.size() == 2){
        std::cout<<"neighbors.size() == 2"<<std::endl;
        if(neighbors[0] == 0){
            json parent_j = json::parse(node->parent->jstring);
            json child_j = json::parse(node->jstring);
            json grandchild_j = json::parse(node->children[neighbors[1]-1]->jstring);
            json children_j;
            json grandchildren_j;
            for(int i = 0; i < node->children.size(); i++){
                if(i+1 == node->children[neighbors[1]-1]->child_id)
                    grandchildren_j.push_back(grandchild_j);
                else
                    grandchildren_j.push_back(nullptr);
            }
            child_j["children"] = grandchildren_j;
            for(int i = 0; i < node->parent->children.size(); i++){
                if(i+1 == node->child_id)
                    children_j.push_back(child_j);
                else
                    children_j.push_back(nullptr);
            }
            parent_j["children"] = children_j;
            js["root"] = parent_j;
        }else{
            // json parent_j = json::parse(node->parent->jstring);
            json child_j = json::parse(node->jstring);
            json grandchild_1j = json::parse(node->children[neighbors[0]-1]->jstring);
            json grandchild_2j = json::parse(node->children[neighbors[1]-1]->jstring);
            // json children_j;
            json grandchildren_j;
            // for(int i = 0; i < node->parent->children.size(); i++){
            //     if(i+1 == node->child_id)
            //         children_j.push_back(child_j);
            //     else
            //         children_j.push_back(nullptr);
            // }
            for(int i = 0; i < node->children.size(); i++){
                if(i+1 == node->children[neighbors[0]-1]->child_id)
                    grandchildren_j.push_back(grandchild_1j);
                else if(i+1 == node->children[neighbors[1]-1]->child_id)
                    grandchildren_j.push_back(grandchild_2j);
                else
                    grandchildren_j.push_back(nullptr);
            }
            child_j["children"] = grandchildren_j;
            // parent_j["children"] = children_j;
            js["root"] = child_j;
        }
    }
    else{
        assert(neighbors.size() > 2);
        std::cout<<"neighbors.size() > 2"<<std::endl;
        if(neighbors[0] == 0){
            json parent_j = json::parse(node->parent->jstring);
            json child_j = json::parse(node->jstring);
            std::unordered_map<int, json> gcid2json;
            for(int idx_nb = 1; idx_nb < neighbors.size(); idx_nb++){
                gcid2json[idx_nb] = json::parse(node->children[neighbors[idx_nb]-1]->jstring);
            }
            // json grandchild_j = json::parse(node->children[neighbors[1]-1]->jstring);
            json children_j;
            json grandchildren_j;
            for(int i = 0; i < node->children.size(); i++){
                if(gcid2json.find(i+1) != gcid2json.end())
                    grandchildren_j.push_back(gcid2json[i+1]);
                else
                    grandchildren_j.push_back(nullptr);
            }
            child_j["children"] = grandchildren_j;
            for(int i = 0; i < node->parent->children.size(); i++){
                if(i+1 == node->child_id)
                    children_j.push_back(child_j);
                else
                    children_j.push_back(nullptr);
            }
            parent_j["children"] = children_j;
            js["root"] = parent_j;
        }else{
            // json parent_j = json::parse(node->parent->jstring);
            json child_j = json::parse(node->jstring);
            // json grandchild_1j = json::parse(node->children[neighbors[0]-1]->jstring);
            // json grandchild_2j = json::parse(node->children[neighbors[1]-1]->jstring);
            std::unordered_map<int, json> gcid2json;
            for(int idx_nb = 1; idx_nb < neighbors.size(); idx_nb++){
                gcid2json[idx_nb] = json::parse(node->children[neighbors[idx_nb]-1]->jstring);
            }
            // json children_j;
            json grandchildren_j;
            for(int i = 0; i < node->children.size(); i++){
                // if(i+1 == node->children[neighbors[0]-1]->child_id)
                //     grandchildren_j.push_back(grandchild_1j);
                // else if(i+1 == node->children[neighbors[1]-1]->child_id)
                //     grandchildren_j.push_back(grandchild_2j);
                // else
                //     grandchildren_j.push_back(nullptr);
                if(gcid2json.find(i+1) != gcid2json.end())
                    grandchildren_j.push_back(gcid2json[i+1]);
                else
                    grandchildren_j.push_back(nullptr);
            }
            child_j["children"] = grandchildren_j;
            // parent_j["children"] = children_j;
            js["root"] = child_j;
        }
    }
    return js.dump();
}

json jsrecords(LoopComplex lc, std::string y_star, std::string y_sub, std::vector<std::string> y_rivals, std::string id){
    json js;
    js["id"] = id;
    js["y_star"] = y_star;
    js["y_sub"] = y_sub;
    js["y_rivals"] = y_rivals;
    js["start"] = lc.node->first;
    js["end"] = lc.node->second;
    js["ipairs"] = lc.ps_inside;
    js["bpairs"] = lc.ps_outside;
    js["motif"] = json::parse(lc.jsmotif(id));
    js["plotstr"] = compose_args4plot(id, y_star, lc.ps_outside, lc.ps_inside);
    return js;
}