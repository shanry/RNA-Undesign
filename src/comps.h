#ifndef COMPS_H
#define COMPS_H

#include<vector>
#include<set>
#include<string>

struct TreeNode {
    int first; // root: -1
    int second; // 
    TreeNode* parent = NULL;
    int child_id = -1;
    std::vector<TreeNode*> children;
    std::string looptype;
    std::vector<int> neighbors;
    std::vector<int> looplens;
    std::string jstring;
    TreeNode(int first_val, int second_val);
    void setLoop();
    void printTree();
    void printTreeEnum(std::string& seq, std::string& y);
    void printTree(std::string ref);
    void printTree(std::string& ref, std::string& seq);
    void printTree(std::string& ref, std::string& seq, std::vector<std::pair<std::string, std::string>>& subrefs);
};

// Define input for decomposition algs

struct LoopComplex {
    int count_uk;
    std::string ref;
    std::string constr;
    int start;
    int end;
    TreeNode* node;
    int left;
    int right;
    std::vector<std::pair<int, int>> ps_outside;
    std::vector<std::pair<int, int>> ps_inside;
    std::vector<int> neighbors;

    LoopComplex(int count, std::string y, std::string cs, int s, int e, TreeNode* n, int l, int r, std::vector<std::pair<int, int>> p_out, std::vector<std::pair<int, int>> p_in);
    void set_neighbors(std::vector<int> nbs);

    std::string jsmotif(std::string id);
};

// Define a basic constraint structure
struct Constraint {
    std::set<int>* indices;
    std::vector<std::string>* seqs;
    std::string structure;

    // Constructor that accepts two arguments to initialize the members
    Constraint(std::set<int>* index_set, std::vector<std::string>* seq_list);

    Constraint(std::set<int>* index_set, std::vector<std::string>* seq_list, std::string& ref);

    void setStructure(std::string& ref);
};

void print(Constraint& cs);
std::string enumerate(std::vector<std::tuple<int, int>>& pairs_diff, ulong order, std::string& seq);

void tree2Loops(TreeNode* root, std::string& ref, std::vector<LoopComplex>& lc_list);
std::string removeNodeFromTree(TreeNode* node, std::string ref);
void tree2TwoNeighbor(TreeNode* root, std::string& ref, std::vector<LoopComplex>& lc_list);
void tree2ThreeNeighbor(TreeNode* root, std::string& ref, std::vector<LoopComplex>& lc_list);
void tree2Edges(TreeNode* root, std::string& ref, std::vector<LoopComplex>& lc_list);
std::string removeEdgeFromTree(TreeNode* node, std::string ref);
void tree2MLoops(TreeNode* root, std::string& ref, std::vector<LoopComplex>& lc_list);
std::string removeMNodeFromTree(TreeNode* node, std::string ref);
std::string removeTwoNeighbors(TreeNode* node, std::string ref, int n1, int n2);
std::string removeThreeNeighbors(TreeNode* node, std::string ref, std::vector<int>& powset);

int max_single(TreeNode* root);
std::string ml_degree(TreeNode* root);


#endif