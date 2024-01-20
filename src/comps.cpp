#include<iostream>
#include<cstdio>
#include<tuple>
#include <map>
#include<utility>
#include <cassert>
#include <cmath>
#include <algorithm>

#include "json.hpp"

#include "comps.h"
#include "utils.h"

char nuc_all[] = "ACGU";
char nuc_pair_all[][3] = {"GC", "CG", "AU", "UA", "GU", "UG"};

extern std::vector<std::pair<int, int>> pairs_outside;
extern std::vector<std::pair<int, int>> pairs_inside;

TreeNode::TreeNode(int first_val, int second_val){
        first = first_val;
        second = second_val;
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
        LoopComplex lc = {count_unknown, cref, constr, start, end, root, root->first, root->second, pairs_outside, pairs_inside};
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
                LoopComplex lc = {count_unknown, cref, constr, start, end, root, root->first, root->second, pairs_outside, pairs_inside};
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
    printf("before if\n");
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
                LoopComplex lc = {count_unknown, cref, constr, start, end, root, root->first, root->second, pairs_outside, pairs_inside};
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
        LoopComplex lc = {count_unknown, cref, constr, start, end, root, root->first, root->second, pairs_outside, pairs_inside};
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
        LoopComplex lc = {count_unknown, cref, constr, start, end, root, root->first, root->second, pairs_outside, pairs_inside};
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
            printf("constr: %s\n", constr.substr(node->parent->first, len_p).c_str());
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
            printf("constr: %s\n", constr.substr(node->parent->first, len_p).c_str());
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

int max_single(TreeNode* root){
    int maxlen = 0;
    if (root->first!=-1&&root->children.size()==1){
        maxlen = (root->children[0]->first - root->first) + (root->second - root->children[0]->second)-2;
    }
    for(TreeNode* child: root->children){
        int maxlen_child = max_single(child);
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