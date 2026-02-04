#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <memory>
#include <numeric>
#include <optional>
#include <queue>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "uniq.h"

const bool verbose = false;

long combinations_count(int n, int r){
	if(r > n || r < 0) return 0;
	if(r == 0 || r == n) return 1;
	r = std::min(r, n - r); // C(n, r) == C(n, n - r)
	long numerator = 1;
	long denominator = 1;
	for(int i = 0; i < r; ++i){
		numerator *= (n - i);
		denominator *= (i + 1);
	}
	return numerator / denominator;
}

std::vector<std::vector<int>> power_set(int start, int set_size, int min_setsize, int max_setsize){
    std::vector<std::vector<int>> powset_all;
    // the size of a power set of n is (2**n)
    unsigned long pow_set_size = pow(2, set_size);
    unsigned long counter;
    // Run from counter 000..0 to 111..1
    for (counter = 0; counter < pow_set_size; counter++) {
        std::vector<int> powset;
        for (int j = 0; j < set_size; j++) {
            // Check if jth bit in the counter is set
            if (counter & (1 << j))
                powset.push_back(j+start);
        }
        if(powset.size() >= min_setsize && powset.size() <= max_setsize)
            powset_all.push_back(powset);
    }
    return powset_all;
}

// read motif, max_prob from file
std::unordered_map<std::string, double> read_motif_prob_max(const std::string &filename) {
    std::unordered_map<std::string, double> motif2prob;
    std::ifstream infile(filename);
    if (!infile) {
        throw std::runtime_error("Failed to open motif file: " + filename);
    }
    std::string motif;
    double max_prob;
	// take the firts two columns (4 columns in total, separated by comma ,)
    std::string line;
	while (std::getline(infile, line)) {
		if (line.empty()) continue;
		std::istringstream iss(line);
		std::string motif, prob_str;
		if (!std::getline(iss, motif, ',')) continue;
		if (!std::getline(iss, prob_str, ',')) continue;
		motif2prob[motif] = std::stod(prob_str);
		// get the rotational variants and set their probabilities too
		std::vector<std::string> rotated_variants = generate_rotated_variants(motif);
		for (const auto &variant : rotated_variants) {
			motif2prob[variant] = std::stod(prob_str);
		}
	}
    return motif2prob;
}


namespace decompose {

// Parameters for motif generation, which control the topology of motifs
const int MAX_WIDTH = 2;
const int MAX_DEPTH = 5;
const int MAX_LOOP_SIZE = 5;

struct Motif;
using Mptr = std::shared_ptr<Motif>;


struct LoopNode : std::enable_shared_from_this<LoopNode> {
	struct PairIndex {
		int open{-1};
		int close{-1};
	};

	using Ptr = std::shared_ptr<LoopNode>;

	PairIndex pair{};
	char type{'u'};  // unknown by default
	std::vector<int> unpaired_bases{};
	std::vector<Ptr> children{};
	std::optional<std::size_t> child_id{};  // index of this node in parent's children vector, if it has no value, then it is root
	std::weak_ptr<LoopNode> parent{};
	// children_within_same_motif can be used to store children that are part of the same motif during motif generation
	std::optional<std::vector<Ptr>> children_gm{}; // children in the same generated motif
	mutable bool motifs_generated{false};
	mutable std::vector<Mptr> motifs{}; // cached motifs rooted at this node
	mutable std::optional<std::size_t> decomposition_count{}; // memoized count
	mutable std::optional<double> min_product_prob{}; // memoized probability for min-product decomposition
	mutable Mptr best_motif{nullptr}; // best motif selected during min-product decomposition

	LoopNode() = default;
	LoopNode(PairIndex p, char loop_type) : pair(std::move(p)), type(loop_type) {}
	LoopNode(PairIndex p, char loop_type, std::vector<int> unpaired_bases, std::vector<Ptr> children_nodes, std::optional<std::size_t> child_id)
		: pair(std::move(p)), type(loop_type), unpaired_bases(std::move(unpaired_bases)), child_id(child_id), children(std::move(children_nodes)) {}

	static Ptr string_to_node(const std::string &structure, char root_type = 'E');
	static Ptr copy_recursive_gm(const Ptr &node);  // deep copy of a LoopNode and its children, recursively

	std::string dot_bracket() const { return to_dot_bracket(); }
	std::string to_dot_bracket() const;
	std::string to_dot_bracket_gm() const;
	std::size_t length() const;
	std::size_t length_verbose(std::ostream &os) const;
	std:: size_t length_gm() const;
	std::vector<const LoopNode *> preorder() const;
	std::vector<const LoopNode *> postordered_motifs() const;
	std::size_t loop_count() const;
	std::size_t loop_count_gm() const;

	// a motif is a contiguous set of loops; generate all motifs growing from current loop node up to given depth, width, and loop size
	std::vector<Mptr> generate_motifs() const;

	// initialize children_gm as a number of 'p' nodes corresponding to the children of this node
	void initialize_children_gm_as_pnodes();

	private:
		static void finalize_node(const Ptr &node);
		static void compute_root_unpaired(const Ptr &root, std::size_t structure_length);
		static std::size_t length_impl(const LoopNode *node, bool verbose, std::ostream *os);
};


struct Motif: std::enable_shared_from_this<Motif> {
	LoopNode::Ptr root{nullptr};  // root loop of the motif
	int depth{1};	// depth of the motif
	// std::vector<LoopNode::Ptr> loops_last_depth;  // loops in the motif

	// using Mptr = std::shared_ptr<Motif>;

	Motif() = default;
	Motif(LoopNode::Ptr r, int d) : root(std::move(r)), depth(d) {}

	// deep copy constructor, copies the root and all other loop nodes
	Motif(const Motif &other) {
		// deep copy of root
		root = LoopNode::copy_recursive_gm(other.root);
		depth = other.depth;
	}

	// bfs to get loops at the last depth
	std::vector<LoopNode::Ptr> get_loops_last_depth() const {
		// fifo queue
		std::queue<LoopNode::Ptr> q_next;
		std::queue<LoopNode::Ptr> q_last;
		q_last.push(root);
		std::vector<LoopNode::Ptr> result;

		while (!q_last.empty()) {
			LoopNode::Ptr node = q_last.front();
			q_last.pop();
			if (node->children_gm && !node->children_gm->empty()) {
				for (const LoopNode::Ptr &child : *(node->children_gm)) {
					if (child->type != 'p') { // skip 'p' nodes
						q_next.push(child);
						result.push_back(child);
					}
				}
			} 
		}
		// std::cout<<"After first level, q_next size: "<< q_next.size() << "\n";
		while (!q_next.empty()) {
			result.clear();
			q_last = q_next;
			q_next = std::queue<LoopNode::Ptr>();
			while(!q_last.empty()) {
				LoopNode::Ptr node = q_last.front();
				result.push_back(node);
				q_last.pop();
				if (node->children_gm && !node->children_gm->empty()) {
					for (const LoopNode::Ptr &child : *(node->children_gm)) {
						if (child->type != 'p') { // skip 'p' nodes
							q_next.push(child);
						}
					}
				} 
			}
		}
		return result;
	}

	// get leaf loops in the motif, exclude the root if it is a leaf
	std::vector<LoopNode::Ptr> get_leaf_loops() const {
		std::vector<LoopNode::Ptr> leaves;
		std::function<void(const LoopNode::Ptr &)> dfs = [&](const LoopNode::Ptr &node) {
			if (!node->children_gm || node->children_gm->empty()) {
				if (node != root and node->type == 'p') {
					leaves.push_back(node);
				}
				return;
			}
			for (const LoopNode::Ptr &child : *(node->children_gm)) {
				dfs(child);
			}
		};
		dfs(root);
		return leaves;
	}

	// get boundary pairs of the motif, include the pairs corresponding to the first loop (the child of the root) and the pairs of leaf loops
	std::vector<LoopNode::PairIndex> get_boundary_pairs() const {
		std::vector<LoopNode::PairIndex> boundaries;
		if (root->children_gm && !root->children_gm->empty()) {
			assert (root->children_gm->size() == 1);
			boundaries.push_back(root->children_gm->front()->pair);
		}
		std::vector<LoopNode::Ptr> leaves = get_leaf_loops();
		for (const LoopNode::Ptr &leaf : leaves) {
			// get original loop first by child_id
			auto parent_ptr = leaf->parent.lock();
			if (parent_ptr) {
				assert(leaf->child_id < parent_ptr->children.size());
				boundaries.push_back(parent_ptr->children[leaf->child_id.value()]->pair);
			}
		}
		return boundaries;
	}
};


LoopNode::Ptr LoopNode::string_to_node(const std::string &structure, char root_type) {
	if (root_type != 'E' && root_type != 'p') {
		throw std::invalid_argument("root_type must be 'E' or 'p'");
	}

	PairIndex root_pair;
	if (root_type == 'E') {
		root_pair = {-1, static_cast<int>(structure.size())};
	} else {
		root_pair = {-1, -1};
	}

	Ptr root = std::make_shared<LoopNode>(root_pair, root_type);

	std::vector<Ptr> stack;
	stack.push_back(root);

	for (int i = 0; i < static_cast<int>(structure.size()); ++i) {
		const char c = structure[i];
		if (c == '(') {
			stack.push_back(std::make_shared<LoopNode>(PairIndex{i, -1}, 'u'));
		} else if (c == ')') {
			if (stack.size() <= 1) {
				throw std::runtime_error("Unbalanced structure: too many closing parentheses");
			}
			Ptr node = stack.back();
			stack.pop_back();
			node->pair.close = i;

			Ptr parent = stack.back();
			node->parent = parent;
			parent->children.push_back(node);
			node->child_id = parent->children.size() - 1;

			if (node->pair.open < 0) {
				throw std::runtime_error("Encountered negative opening index");
			}

			finalize_node(node);
		}
	}

	if (stack.size() != 1) {
		throw std::runtime_error("Unbalanced structure: missing closing parentheses");
	}

	compute_root_unpaired(root, structure.size());
	return root;
}

LoopNode::Ptr LoopNode::copy_recursive_gm(const Ptr &node) { // parent will be set by the caller
	if (!node) {
		return nullptr;
	}
	// std::cout<<"Copying node (" << node->pair.open << ", " << node->pair.close << ") of type " << node->type << "\n";
	Ptr new_node = std::make_shared<LoopNode>(node->pair, node->type, node->unpaired_bases, node->children, node->child_id);
	// copy children_gm recursively
	if (node->children_gm) {
		assert(node->children_gm->size() == node->children.size());
		new_node->children_gm = std::vector<Ptr>{};
		new_node->children_gm->reserve(node->children_gm->size());
		for (const Ptr &child_gm : *(node->children_gm)) {
			Ptr child_gm_copy = copy_recursive_gm(child_gm);
			if (child_gm_copy) {
				child_gm_copy->parent = new_node;
				child_gm_copy->child_id = new_node->children_gm->size();
				new_node->children_gm->push_back(child_gm_copy);
			}
		}
	}
	return new_node;
}

void LoopNode::finalize_node(const Ptr &node) {
	const int span = node->pair.close - node->pair.open;
	if (span <= 2) {
		node->type = 'p';
		node->unpaired_bases.clear();
		return;
	}

	const std::size_t child_count = node->children.size();
	if (child_count == 0) {
		node->type = 'H';
		node->unpaired_bases = {node->pair.close - node->pair.open - 1};
	} else if (child_count == 1) {
		node->type = 'I';
		const Ptr &child = node->children.front();
		node->unpaired_bases = {
			child->pair.open - node->pair.open - 1,
			node->pair.close - child->pair.close - 1,
		};
	} else {
		node->type = 'M';
		node->unpaired_bases.clear();
		node->unpaired_bases.push_back(node->children.front()->pair.open - node->pair.open - 1);
		for (std::size_t j = 0; j + 1 < child_count; ++j) {
			node->unpaired_bases.push_back(
				node->children[j + 1]->pair.open - node->children[j]->pair.close - 1);
		}
		node->unpaired_bases.push_back(
			node->pair.close - node->children.back()->pair.close - 1);
	}
}

void LoopNode::compute_root_unpaired(const Ptr &root, std::size_t structure_length) {
	root->unpaired_bases.clear();
	if (root->children.empty()) {
		root->unpaired_bases.push_back(static_cast<int>(structure_length));
		return;
	}

	root->unpaired_bases.push_back(root->children.front()->pair.open - root->pair.open - 1);
	for (std::size_t j = 0; j + 1 < root->children.size(); ++j) {
		root->unpaired_bases.push_back(
			root->children[j + 1]->pair.open - root->children[j]->pair.close - 1);
	}
	root->unpaired_bases.push_back(root->pair.close - root->children.back()->pair.close - 1);
}

std::string LoopNode::to_dot_bracket() const {
	if (type == 'p') {
		if (!child_id.has_value()) {
			// assert children.size() == 1
			assert(children.size() == 1);
			return children.front()->to_dot_bracket();
		}
		return "(*)";
	}

	if (type == 'H') {
		const int loop = unpaired_bases.empty() ? 0 : unpaired_bases.front();
		return "(" + std::string(loop, '.') + ")";
	}

	std::string result;
	if (type != 'E') {
		result.push_back('(');
	}

	if (!unpaired_bases.empty()) {
		result.append(unpaired_bases.front(), '.');
	}

	for (std::size_t i = 0; i < children.size(); ++i) {
		result += children[i]->to_dot_bracket();
		if (i + 1 < unpaired_bases.size()) {
			result.append(unpaired_bases[i + 1], '.');
		} else {
			throw std::runtime_error("Mismatch between children and unpaired bases");
		}
	}

	if (type != 'E') {
		result.push_back(')');
	}

	return result;
}

std::string LoopNode::to_dot_bracket_gm() const {
	if (type == 'p') {
		if (!child_id.has_value()) {
			// assert children.size() == 1
			assert(children_gm->size() == 1);
			return children_gm->front()->to_dot_bracket_gm();
		}
		return "(*)";
	}

	if (type == 'H') {
		const int loop = unpaired_bases.empty() ? 0 : unpaired_bases.front();
		return "(" + std::string(loop, '.') + ")";
	}

	std::string result;
	if (type != 'E') {
		result.push_back('(');
	}else{
		result.push_back('5'); // to indicate external loop
	}

	if (!unpaired_bases.empty()) {
		result.append(unpaired_bases.front(), '.');
	}

	const std::vector<Ptr> &children_to_use = children_gm ? *children_gm : children;

	for (std::size_t i = 0; i < children_to_use.size(); ++i) {
		result += children_to_use[i]->to_dot_bracket_gm();
		if (i + 1 < unpaired_bases.size()) {
			result.append(unpaired_bases[i + 1], '.');
		} else {
			throw std::runtime_error("Mismatch between children and unpaired bases");
		}
	}

	if (type != 'E') {
		result.push_back(')');
	}else{
		result.push_back('3'); // to indicate external loop
	}

	return result;
}

std::size_t LoopNode::length_impl(const LoopNode *node, bool verbose, std::ostream *os) {
	std::size_t subtotal = 0;
	if (verbose && os != nullptr) {
		*os << "Computing length for node (" << node->pair.open << ", " << node->pair.close
			<< ") of type " << node->type << '\n';
	}

	for (const Ptr &child : node->children) {
		subtotal += length_impl(child.get(), verbose, os);
	}

	if (node->type == 'p') {
		if (node->child_id.has_value()) {
			subtotal += 2;
		}
	} else if (node->type == 'E') {
		subtotal += std::accumulate(node->unpaired_bases.begin(), node->unpaired_bases.end(), std::size_t{0});
	} else if (node->type == 'H' || node->type == 'I' || node->type == 'M') {
		subtotal += std::accumulate(node->unpaired_bases.begin(), node->unpaired_bases.end(), std::size_t{0}) + 2;
	} else {
		throw std::runtime_error("Unknown node type encountered during length computation");
	}

	if (verbose && os != nullptr) {
		*os << "Current length for node (" << node->pair.open << ", " << node->pair.close
			<< ") is " << subtotal << '\n';
	}

	return subtotal;
}

std::size_t LoopNode::length() const {
	return length_impl(this, false, nullptr);
}

std::size_t LoopNode::length_verbose(std::ostream &os) const {
	return length_impl(this, true, &os);
}

std:: size_t LoopNode::length_gm() const {
	std::function<std::size_t(const LoopNode *)> length_gm_impl = [&](const LoopNode *node) -> std::size_t {
		std::size_t subtotal = 0;
		const std::vector<Ptr> &children_to_use = node->children_gm ? *node->children_gm : node->children;

		for (const Ptr &child : children_to_use) {
			subtotal += length_gm_impl(child.get());
		}

		if (node->type == 'p') {
			if (node->child_id.has_value()) {
				subtotal += 2;
			}
		} else if (node->type == 'E') {
			subtotal += std::accumulate(node->unpaired_bases.begin(), node->unpaired_bases.end(), std::size_t{0});
		} else if (node->type == 'H' || node->type == 'I' || node->type == 'M') {
			subtotal += std::accumulate(node->unpaired_bases.begin(), node->unpaired_bases.end(), std::size_t{0}) + 2;
		} else {
			throw std::runtime_error("Unknown node type encountered during length computation");
		}

		return subtotal;
	};

	return length_gm_impl(this);
}

std::vector<const LoopNode *> LoopNode::preorder() const {
	std::vector<const LoopNode *> nodes;
	std::vector<const LoopNode *> stack{this};
	while (!stack.empty()) {
		const LoopNode *node = stack.back();
		stack.pop_back();
		if (node->type != 'E' && node->type != 'p') {
			nodes.push_back(node);
		}
		for (auto it = node->children.rbegin(); it != node->children.rend(); ++it) {
			stack.push_back(it->get());
		}
	}
	return nodes;
}

void LoopNode::initialize_children_gm_as_pnodes() {
	assert (children_gm == std::nullopt); // should be uninitialized
	children_gm = std::vector<Ptr>{};
	for (const Ptr &child : children) {
		Ptr child_copy = std::make_shared<LoopNode>();
		child_copy->type = 'p'; // leaf node
		child_copy->parent = shared_from_this();
		child_copy->child_id = children_gm->size();
		children_gm->push_back(Ptr(child_copy));
	}
}

std::vector<Mptr> LoopNode::generate_motifs() const {
	std::vector<Mptr> motifs;
	int count_motif_last_depth = 0;
	int count_motif_current_depth = 0;
	// Implementation of motif generation would go here
	for (int depth = 1; depth <= MAX_DEPTH; ++depth) {
		if (depth == 1) {  // single loop motif, root node of type 'p'
			// LoopNode *motif_root = new LoopNode();
			std::shared_ptr<LoopNode> motif_root = std::make_shared<LoopNode>();
			motif_root->type = 'p';
			// put a copy of this node to children_gm
			std::shared_ptr<LoopNode> current_copy = std::make_shared<LoopNode>(*this);
			// std::cout<<"number of children of current node: "<< current_copy->children.size() << "\n";
			current_copy->parent = motif_root->shared_from_this();  // set parent by the caller
			// set each child of current_copy to be a new node of type 'p'
			current_copy->initialize_children_gm_as_pnodes();
			motif_root->children_gm = std::vector<Ptr>{};
			motif_root->children_gm->push_back(current_copy); // motif_root has one child: current_copy
			std::shared_ptr<const LoopNode> self_const = shared_from_this();
			Ptr original_node = std::const_pointer_cast<LoopNode>(self_const);
			motif_root->children.push_back(original_node); // motif_root also retains pointer to original node
			current_copy->child_id = 0; // first child of motif_root
			Mptr motif = std::make_shared<Motif>(Ptr(motif_root), depth);
			motifs.push_back(motif);
			count_motif_current_depth = 1;
		}else{
			// more complex motifs can be generated here; first get the motifs of depth - 1
			count_motif_last_depth = count_motif_current_depth;
			count_motif_current_depth = 0;
			// for each motif of last depth, try to add children to increase depth
			int index_start = motifs.size() - count_motif_last_depth;
			if (verbose) {
				std::cout << "Generating motifs of depth " << depth << " from "
					  << count_motif_last_depth << " motifs of depth " << depth - 1 << "\n";
				std::cout << "Index range: " << index_start << " to " << motifs.size() - 1 << "\n";
			}
			for (int i = index_start; i < index_start + count_motif_last_depth; ++i) {
				// std::cout << "Processing motif index " << i << "\n";
				Mptr motif = motifs[i];
				int count_loop = motif->root->loop_count_gm();
				int count_more_loop = MAX_LOOP_SIZE - count_loop;
				if (verbose)
					std::cout << "Motif " << i << " has " << count_loop << " loops, can add " << count_more_loop << " more loops\n";
				if (count_more_loop <= 0) {
					continue;  // skip motifs that have reached max loop size
				}
				std::vector<LoopNode::Ptr> loops_last_depth = motif->get_loops_last_depth();
				if (verbose)
					std::cout << "Motif " << i << " has " << loops_last_depth.size() << " loops at last depth\n";
				int count_children_sum = 0;
				for (const Ptr &loop : loops_last_depth) {
					count_children_sum += loop->children.size();  // consider the children in the original structure
				}
				if (verbose)
					std::cout << "Total children available to add from last depth loops: " << count_children_sum << "\n";
				if (count_children_sum == 0) {
					continue;  // no children to add
				}
				int max_children_to_add = std::min({MAX_WIDTH, count_more_loop, count_children_sum});
				if (verbose){
					std::cout << "Generating motifs from motif " << i << " of depth " << depth - 1
						  << " with " << loops_last_depth.size() << " loops at last depth, max children to add: "
						  << max_children_to_add << "\n";
				}
				// generate all pairs of <loop_last_depth_index, child_index> to add
				std::vector<std::pair<int, int>> loop_child_pairs;
				for (std::size_t j = 0; j < loops_last_depth.size(); ++j) {
					const Ptr &loop = loops_last_depth[j];
					for (std::size_t k = 0; k < loop->children.size(); ++k) {
						loop_child_pairs.emplace_back(j, k);
					}
				}
				// power set generation to add children
				std::vector<std::vector<int>> powerset_indices = power_set(0, loop_child_pairs.size(), 1, max_children_to_add);
				std::vector<std::set<std::pair<int, int>>> powerset;
				powerset.reserve(powerset_indices.size());
				for (const auto &indices : powerset_indices) {
					std::set<std::pair<int, int>> subset;
					for (int index : indices) {
						subset.insert(loop_child_pairs[index]);
					}
					powerset.push_back(std::move(subset));
				}
				// buggy version, to be fixed
				for (const auto &subset : powerset) {
					// create a new motif by adding the selected children
					// std::cout << "Generating new motif by adding " << subset.size() << " children\n";
					Mptr new_motif = std::make_shared<Motif>(*motif);  // deep copy
					std::vector<LoopNode::Ptr> new_loops_last_depth = new_motif->get_loops_last_depth();
					// grow children from loops of last depth at the indices of subset
					for (const auto &pair : subset) {
						int loop_index = pair.first;
						int child_index = pair.second;
						Ptr loop_node = new_loops_last_depth[loop_index];
						// get the corresponding child in the original structure
						// Ptr original_loop_node = loops_last_depth[loop_index];
						if (child_index >= loop_node->children.size()) {
							throw std::runtime_error("Child index out of bounds during motif generation");
						}
						// make a copy of the child to add
						Ptr child_to_add = loop_node->children[child_index];
						// make a copy of child_to_add and set its parent and child_id
						Ptr child_to_add_copy = std::make_shared<LoopNode>(*child_to_add);
						child_to_add_copy->initialize_children_gm_as_pnodes();
						// replace the child in children_gm at child_index
						assert (loop_node->children_gm);
						if (child_index >= loop_node->children_gm->size()) {
							throw std::runtime_error("Child index out of bounds in children_gm during motif generation");
						}
						(*loop_node->children_gm)[child_index] = child_to_add_copy;
						child_to_add_copy->parent = loop_node;
						child_to_add_copy->child_id = child_index; // keep the same index as in original
					}
					new_motif->depth = depth;
					motifs.push_back(new_motif);
					count_motif_current_depth += 1;
				}
			}
		}
		if (verbose && count_motif_current_depth > 0)
			std::cout << "Generated " << count_motif_current_depth << " motifs of depth " << depth << "\n";
	}
	return motifs;
}

std::size_t LoopNode::loop_count() const {
	std::function<std::size_t(const LoopNode *)> traverse = [&](const LoopNode *node) -> std::size_t {
		std::size_t total = (node->type == 'p') ? 0 : 1;
		for (const Ptr &child : node->children) {
			total += traverse(child.get());
		}
		return total;
	};
	return traverse(this);
}

std::size_t LoopNode::loop_count_gm() const {
	std::function<std::size_t(const LoopNode *)> traverse = [&](const LoopNode *node) -> std::size_t {
		std::size_t total = (node->type == 'p') ? 0 : 1;
		const std::vector<Ptr> &children_to_use = node->children_gm ? *node->children_gm : node->children;
		for (const Ptr &child : children_to_use) {
			total += traverse(child.get());
		}
		return total;
	};
	return traverse(this);
}
} // namespace decompose


int show_decomposition(std::string structure) {
	// print max_depth, max_width, max_loop_size
	std::cout << "Motif generation parameters: MAX_DEPTH=" << decompose::MAX_DEPTH
				<< ", MAX_WIDTH=" << decompose::MAX_WIDTH
				<< ", MAX_LOOP_SIZE=" << decompose::MAX_LOOP_SIZE << '\n';
	try {
		auto root = decompose::LoopNode::string_to_node(structure);
		std::cout << "Dot-bracket: " << root->dot_bracket() << '\n';
		std::cout << "Length: " << root->length() << '\n';
		std::cout << "Loop count: " << root->loop_count() << '\n';

		// generate motifs
		std::queue<decompose::LoopNode::Ptr> q;
		q.push(root);
		// bfs to generate all motifs from each node
		while (!q.empty()) {
			decompose::LoopNode::Ptr current_node = q.front();
			q.pop();
			std::cout << "----------------------------------------\n";
			// generate motifs from current_node
			std::vector<decompose::Mptr> node_motifs = current_node->generate_motifs();
			std::cout << "From node (" << current_node->pair.open << ", " << current_node->pair.close << "): ";
			std::cout << "Generated " << node_motifs.size() << " motifs from node (" 
					  << current_node->pair.open << ", " << current_node->pair.close << ")\n";
			for (std::size_t i = 0; i < node_motifs.size(); ++i) {
				const decompose::Mptr &motif = node_motifs[i];
				std::cout << "Motif " << i + 1 << ": depth=" << motif->depth
						  << ", dot-bracket=" << motif->root->to_dot_bracket_gm()
						  << ", length=" << motif->root->length_gm() << '\n';
			}
			// enqueue children
			for (const decompose::LoopNode::Ptr &child : current_node->children) {
				q.push(child);
			}
		}

	} catch (const std::exception &ex) {
		std::cerr << "Error: " << ex.what() << std::endl;
		return 1;
	}
	return 0;
}


int random_decompose(std::string structure) {
	try {
		// fix random seed for reproducibility
		srand(time(nullptr));
		auto root = decompose::LoopNode::string_to_node(structure);
		std::cout << "Dot-bracket: " << root->dot_bracket() << '\n';
		std::cout << "Length: " << root->length() << '\n';
		std::cout << "Loop count: " << root->loop_count() << '\n';

		std::vector<decompose::Mptr> decomposed_motifs;
		// auto create_single_loop_motif = [](const decompose::LoopNode::Ptr &node) -> decompose::Mptr {
		// 	auto motif_root = std::make_shared<decompose::LoopNode>();
		// 	motif_root->type = 'p';
		// 	auto node_copy = std::make_shared<decompose::LoopNode>(*node);
		// 	node_copy->parent = motif_root;
		// 	node_copy->initialize_children_gm_as_pnodes();
		// 	motif_root->children_gm = std::vector<decompose::LoopNode::Ptr>{};
		// 	motif_root->children_gm->push_back(node_copy);
		// 	motif_root->children.push_back(node);
		// 	node_copy->child_id = 0;
		// 	return std::make_shared<decompose::Motif>(motif_root, 1);
		// };

		// generate motifs
		std::queue<decompose::LoopNode::Ptr> q;
		q.push(root);
		std::set<std::pair<int, int>> visited;
		// bfs to generate all motifs from each node
		while (!q.empty()) {
			decompose::LoopNode::Ptr current_node = q.front();
			q.pop();
			if (!current_node) {
				continue;
			}
			std::pair<int, int> current_key{current_node->pair.open, current_node->pair.close};
			if (!visited.insert(current_key).second) {
				continue; // already processed
			}
			std::cout << "----------------------------------------\n";
			// generate motifs from current_node
			std::vector<decompose::Mptr> node_motifs = current_node->generate_motifs();
			std::cout << "From node (" << current_node->pair.open << ", " << current_node->pair.close << "): ";
			std::cout << "Generated " << node_motifs.size() << " motifs from node (" 
					  << current_node->pair.open << ", " << current_node->pair.close << ")\n";
			for (std::size_t i = 0; i < node_motifs.size(); ++i) {
				const decompose::Mptr &motif = node_motifs[i];
				std::cout << "Motif " << i + 1 << ": depth=" << motif->depth
						  << ", dot-bracket=" << motif->root->to_dot_bracket_gm()
						  << ", length=" << motif->root->length_gm() << '\n';
			}
			// select a random motif to continue decomposition
			if (!node_motifs.empty()) {
				std::size_t random_index = rand() % node_motifs.size();
				decompose::Mptr selected_motif = node_motifs[random_index];
				std::cout << "Randomly selected motif " << random_index + 1 << " for further decomposition.\n";
				decomposed_motifs.push_back(selected_motif);
				// if selected motif has depth 1 and is a leaf node, skip further decomposition
				
				// enqueue leaf nodes of the selected motif
				std::vector<decompose::LoopNode::Ptr> leaf_loops = selected_motif->get_leaf_loops();
				for (const decompose::LoopNode::Ptr &leaf : leaf_loops) {
					std::cout << "Processing leaf loop (" << leaf->pair.open << ", " << leaf->pair.close << ")\n";
					if (!leaf) {
						continue;
					}
					auto parent = leaf->parent.lock();
					if (!parent || !leaf->child_id) {
						continue;
					}
					std::size_t child_index = leaf->child_id.value();
					if (child_index >= parent->children.size()) {
						continue;
					}
					auto node_in_original = parent->children[child_index];
					std::cout<< "Corresponding node in original structure: (" << node_in_original->pair.open << ", " << node_in_original->pair.close << ")\n";
					if (!node_in_original) {
						continue;
					}
					// if (node_in_original->children.empty()) {
					// 	// decomposed_motifs.push_back(create_single_loop_motif(node_in_original));
					// 	continue;
					// }
					std::pair<int, int> child_key{node_in_original->pair.open, node_in_original->pair.close};
					if (visited.count(child_key) > 0) {
						continue;
					}
					std::cout << "Enqueuing leaf loop (" << node_in_original->pair.open << ", " << node_in_original->pair.close << ")\n";
					q.push(node_in_original);
				}
			}
		}
		std::cout << "structure length: " << structure.length() << "\n";
		int sum_of_motif_lengths = 0;
		std::cout << "\nFinal decomposed motifs:\n";
		for (std::size_t i = 0; i < decomposed_motifs.size(); ++i) {
			const decompose::Mptr &motif = decomposed_motifs[i];
			std::cout << "Motif " << i + 1 << ": depth=" << motif->depth
					  << ", dot-bracket=" << motif->root->to_dot_bracket_gm()
					  << ", length=" << motif->root->length_gm() << '\n';
			sum_of_motif_lengths += motif->root->length_gm();
			// minus number of * time 2
			int count_stars = 0;
			std::string db = motif->root->to_dot_bracket_gm();
			for (char c : db) {
				if (c == '*') {
					count_stars += 1;
				}
			}
			sum_of_motif_lengths -= count_stars * 2;
		}
		std::cout << "Sum of motif lengths (adjusted): " << sum_of_motif_lengths << "\n";
		assert (sum_of_motif_lengths == structure.length());
	} catch (const std::exception &ex) {
		std::cerr << "Error: " << ex.what() << std::endl;
		// throw ex;
		std::cout << "Exception caught in random_decompose.\n";
		throw ex;
		return 1;
	}
	return 0;
}

// count number of possible decompositions via recursion, count all = sum of (count child_1 * count child_2 * ...) for each current decomposition; base case: leaf node return 1
size_t count_decompositions(const decompose::LoopNode::Ptr &node) {
	// paired-base leaves (`p`) are atomic and contribute a single decomposition
	if (node->type == 'p') {
		return 1;
	}
	if (node->decomposition_count) {
		return *node->decomposition_count;
	}
	// std::cout << "Counting decompositions for node (" << node->pair.open << ", " << node->pair.close << ") of type " << node->type << "\n";
	size_t count = 1;
	if (node->type != 'H') {
		// generate motifs from current node if they were not generated
		if (!node->motifs_generated) {
			// std::cout << "Generating motifs for node (" << node->pair.open << ", " << node->pair.close << ")\n";
			node->motifs = node->generate_motifs();
			node->motifs_generated = true;
			// std::cout << "Generated " << node->motifs.size() << " motifs for node (" << node->pair.open << ", " << node->pair.close << ")\n";
		}
		const std::vector<decompose::Mptr> &node_motifs = node->motifs;
		std::cout << "=>Counting decompositions from " << node_motifs.size() << " motifs for node (" << node->pair.open << ", " << node->pair.close << ")\n";
		for (const decompose::Mptr &motif : node_motifs) {
			size_t count_children = 1;
			// get leaf loops of motif->root
			const std::vector<decompose::LoopNode::Ptr> &children_to_use = motif->get_leaf_loops();
			// if no leaf loops, continue
			if (children_to_use.empty()) {
				std::cout << "Motif at node (" << node->pair.open << ", " << node->pair.close << ") has no leaf loops, skipping\n";
				continue;
			}
			std::cout << "Processing motif at node (" << node->pair.open << ", " << node->pair.close << ") with "
					<< children_to_use.size() << " leaf loops\n";
			// for each child, get corresponding node in original structure and recursive call
			for (const decompose::LoopNode::Ptr &child : children_to_use) {
				// get corresponding child in original structure
				auto parent = child->parent.lock();
				if (!parent || !child->child_id) {
					continue;
				}
				std::size_t child_index = child->child_id.value();
				if (child_index >= parent->children.size()) {
					continue;
				}
				auto node_in_original = parent->children[child_index];
				// recursive call
				count_children *= count_decompositions(node_in_original);
			}
			count += count_children;
		}
		std::cout << "<=Total decompositions for node (" << node->pair.open << ", " << node->pair.close << "): " << count << "\n";
	}
	node->decomposition_count = count;
	return count;
}


// utilize dynamic programming to identify the decomposition with the minimum product of motif probabilities, which is the upper bound of the probability of the entire structure
double min_product_decomposition(const decompose::LoopNode::Ptr &node, const std::unordered_map<std::string, double> &motif2prob) {
	// base case: paired-base leaves (`p`) have probability 1
	// if (node->type == 'p') {
	// 	node->min_product_prob = 1.0;
	// 	return 1.0;
	// }
	assert (node->type != 'p'); // should not be called on paired-base leaves
	// check if already computed
	if (node->min_product_prob) {
		return *node->min_product_prob;
	}
	double min_product = 1.0;
	if (node->type != 'H') {
		// generate motifs from current node if they were not generated
		if (!node->motifs_generated) {
			node->motifs = node->generate_motifs();
			node->motifs_generated = true;
		}
		const std::vector<decompose::Mptr> &node_motifs = node->motifs;
		for (const decompose::Mptr &motif : node_motifs) {
			double motif_prob = 1.0;
			std::string motif_db = motif->root->to_dot_bracket_gm();
			// check motif probability in motif2prob
			auto it = motif2prob.find(motif_db);
			if (it != motif2prob.end()) {
				motif_prob = it->second;
			} else {
				// assign a default probability for unknown motifs
				motif_prob = 1.0;
			}
			double product_children = 1.0;
			// get leaf loops of motif->root
			const std::vector<decompose::LoopNode::Ptr> &children_to_use = motif->get_leaf_loops();
			// for each child, get corresponding node in original structure and recursive call
			for (const decompose::LoopNode::Ptr &child : children_to_use) {
				// get corresponding child in original structure
				auto parent = child->parent.lock();
				if (!parent || !child->child_id) {
					continue;
				}
				std::size_t child_index = child->child_id.value();
				if (child_index >= parent->children.size()) {
					continue;
				}
				auto node_in_original = parent->children[child_index];
				// recursive call
				product_children *= min_product_decomposition(node_in_original, motif2prob);
			}
			double total_product = motif_prob * product_children;
			if (total_product <= min_product) {
				min_product = total_product;
				// set backpointer for reconstruction if needed
				node->best_motif = motif;
			}
		}
	}else{
		// hairpin loop, assign probability 1.0
		double motif_prob = 1.0;
		double product_children = 1.0; // no children
		double total_product = motif_prob * product_children;
		min_product = total_product;
		// set backpointer for reconstruction if needed
		// create a motif for hairpin loop
		auto motif_root = std::make_shared<decompose::LoopNode>();
		motif_root->type = 'p';
		auto node_copy = std::make_shared<decompose::LoopNode>(*node);
		node_copy->parent = motif_root;
		// no children to initialize
		motif_root->children_gm = std::vector<decompose::LoopNode::Ptr>{};
		motif_root->children_gm->push_back(node_copy);
		// motif_root->children.push_back(node);
		node_copy->child_id = 0;
		decompose::Mptr hairpin_motif = std::make_shared<decompose::Motif>(motif_root, 1);
		node->best_motif = hairpin_motif;
	}
	node->min_product_prob = min_product;
	return min_product;
}


double min_product_decomposition_log(const decompose::LoopNode::Ptr &node, const std::unordered_map<std::string, double> &motif2prob) {
	// base case: paired-base leaves (`p`) have probability 1
	// if (node->type == 'p') {
	// 	node->min_product_prob = 1.0;
	// 	return 1.0;
	// }
	assert (node->type != 'p'); // should not be called on paired-base leaves
	// check if already computed
	if (node->min_product_prob) {
		return *node->min_product_prob;
	}
	double min_product = 0.;   // log(1) = 0
	if (node->type != 'H') {
		// generate motifs from current node if they were not generated
		if (!node->motifs_generated) {
			node->motifs = node->generate_motifs();
			node->motifs_generated = true;
		}
		const std::vector<decompose::Mptr> &node_motifs = node->motifs;
		for (const decompose::Mptr &motif : node_motifs) {
			double motif_prob = 0.; // log(1) = 0
			std::string motif_db = motif->root->to_dot_bracket_gm();
			// check motif probability in motif2prob
			auto it = motif2prob.find(motif_db);
			if (it != motif2prob.end()) {
				double motif_prob_original = it->second;
				if (motif_prob_original <= 0.0) {
					std::cerr << "Warning: motif probability is non-positive for motif: " << motif_db << "\n";
					motif_prob_original = 1; // assign a small positive probability to avoid log(0)
				}
				motif_prob = std::log(motif_prob_original);
			} 
			double product_children = 0.; // log(1) = 0
			// get leaf loops of motif->root
			const std::vector<decompose::LoopNode::Ptr> &children_to_use = motif->get_leaf_loops();
			// for each child, get corresponding node in original structure and recursive call
			for (const decompose::LoopNode::Ptr &child : children_to_use) {
				// get corresponding child in original structure
				auto parent = child->parent.lock();
				if (!parent || !child->child_id) {
					continue;
				}
				std::size_t child_index = child->child_id.value();
				if (child_index >= parent->children.size()) {
					continue;
				}
				auto node_in_original = parent->children[child_index];
				// recursive call
				product_children += min_product_decomposition_log(node_in_original, motif2prob);
			}
			double total_product = motif_prob + product_children;
			if (total_product <= min_product) {
				min_product = total_product;
				// set backpointer for reconstruction if needed
				node->best_motif = motif;
			}
		}
	}else{
		// hairpin loop, assign probability 1.0
		double motif_prob = 0.; // log(1) = 0
		double product_children = 0.; // log(1) = 0, no children
		double total_product = motif_prob + product_children;
		min_product = total_product;
		// set backpointer for reconstruction if needed
		// create a motif for hairpin loop
		auto motif_root = std::make_shared<decompose::LoopNode>();
		motif_root->type = 'p';
		auto node_copy = std::make_shared<decompose::LoopNode>(*node);
		node_copy->parent = motif_root;
		// no children to initialize
		motif_root->children_gm = std::vector<decompose::LoopNode::Ptr>{};
		motif_root->children_gm->push_back(node_copy);
		// motif_root->children.push_back(node);
		node_copy->child_id = 0;
		decompose::Mptr hairpin_motif = std::make_shared<decompose::Motif>(motif_root, 1);
		node->best_motif = hairpin_motif;
	}
	node->min_product_prob = min_product;
	return min_product;
}


double random_product_decomposition_log(const decompose::LoopNode::Ptr &node, const std::unordered_map<std::string, double> &motif2prob) {
	assert (node->type != 'p'); // should not be called on paired-base leaves
	// check if already computed
	if (node->min_product_prob) {
		return *node->min_product_prob;
	}
	double min_product = std::log(static_cast<double>(rand()) / RAND_MAX);   // log(1) = 0
	if (node->type != 'H') {
		// generate motifs from current node if they were not generated
		if (!node->motifs_generated) {
			node->motifs = node->generate_motifs();
			node->motifs_generated = true;
		}
		const std::vector<decompose::Mptr> &node_motifs = node->motifs;
		for (const decompose::Mptr &motif : node_motifs) {
			double motif_prob = std::log(static_cast<double>(rand()) / RAND_MAX); // log(1) = 0
			double product_children = std::log(static_cast<double>(rand()) / RAND_MAX); // log(1) = 0
			// get leaf loops of motif->root
			const std::vector<decompose::LoopNode::Ptr> &children_to_use = motif->get_leaf_loops();
			// for each child, get corresponding node in original structure and recursive call
			for (const decompose::LoopNode::Ptr &child : children_to_use) {
				// get corresponding child in original structure
				auto parent = child->parent.lock();
				if (!parent || !child->child_id) {
					continue;
				}
				std::size_t child_index = child->child_id.value();
				if (child_index >= parent->children.size()) {
					continue;
				}
				auto node_in_original = parent->children[child_index];
				// recursive call
				product_children += random_product_decomposition_log(node_in_original, motif2prob);
			}
			double total_product = motif_prob + product_children;
			if (total_product <= min_product) {
				min_product = total_product;
				// set backpointer for reconstruction if needed
				node->best_motif = motif;
			}
		}
	}else{
		// hairpin loop, assign probability 1.0
		double motif_prob = std::log(static_cast<double>(rand()) / RAND_MAX); // log(1) = 0
		double product_children = std::log(static_cast<double>(rand()) / RAND_MAX); // log(1) = 0, no children
		double total_product = motif_prob + product_children;
		min_product = total_product;
		// set backpointer for reconstruction if needed
		// create a motif for hairpin loop
		auto motif_root = std::make_shared<decompose::LoopNode>();
		motif_root->type = 'p';
		auto node_copy = std::make_shared<decompose::LoopNode>(*node);
		node_copy->parent = motif_root;
		// no children to initialize
		motif_root->children_gm = std::vector<decompose::LoopNode::Ptr>{};
		motif_root->children_gm->push_back(node_copy);
		// motif_root->children.push_back(node);
		node_copy->child_id = 0;
		decompose::Mptr hairpin_motif = std::make_shared<decompose::Motif>(motif_root, 1);
		node->best_motif = hairpin_motif;
	}
	node->min_product_prob = min_product;
	return min_product;
}

// baseline: only consider the motif with the minimum probability
double min_min_decomposition_log(const decompose::LoopNode::Ptr &node, const std::unordered_map<std::string, double> &motif2prob) {
	// base case: paired-base leaves (`p`) have probability 1
	// if (node->type == 'p') {
	// 	node->min_product_prob = 1.0;
	// 	return 1.0;
	// }
	assert (node->type != 'p'); // should not be called on paired-base leaves
	// check if already computed
	if (node->min_product_prob) {
		return *node->min_product_prob;
	}
	double min_product = 0.;   // log(1) = 0
	if (node->type != 'H') {
		// generate motifs from current node if they were not generated
		if (!node->motifs_generated) {
			node->motifs = node->generate_motifs();
			node->motifs_generated = true;
		}
		const std::vector<decompose::Mptr> &node_motifs = node->motifs;
		for (const decompose::Mptr &motif : node_motifs) {
			double motif_prob = 0.; // log(1) = 0
			std::string motif_db = motif->root->to_dot_bracket_gm();
			// check motif probability in motif2prob
			auto it = motif2prob.find(motif_db);
			if (it != motif2prob.end()) {
				double motif_prob_original = it->second;
				if (motif_prob_original <= 0.0) {
					std::cerr << "Warning: motif probability is non-positive for motif: " << motif_db << "\n";
					motif_prob_original = 1; // assign a small positive probability to avoid log(0)
				}
				motif_prob = std::log(motif_prob_original);
			} 
			double product_children = 0.; // log(1) = 0
			// get leaf loops of motif->root
			const std::vector<decompose::LoopNode::Ptr> &children_to_use = motif->get_leaf_loops();
			// for each child, get corresponding node in original structure and recursive call
			for (const decompose::LoopNode::Ptr &child : children_to_use) {
				// get corresponding child in original structure
				auto parent = child->parent.lock();
				if (!parent || !child->child_id) {
					continue;
				}
				std::size_t child_index = child->child_id.value();
				if (child_index >= parent->children.size()) {
					continue;
				}
				auto node_in_original = parent->children[child_index];
				// recursive call
				product_children = std::min(product_children, min_min_decomposition_log(node_in_original, motif2prob));
			}
			double total_product = std::min(motif_prob, product_children);
			if (total_product <= min_product) {
				min_product = total_product;
				// set backpointer for reconstruction if needed
				node->best_motif = motif;
			}
		}
	}else{
		// hairpin loop, assign probability 1.0
		double motif_prob = 0.; // log(1) = 0
		double product_children = 0.; // log(1) = 0, no children
		double total_product = std::min(motif_prob, product_children);
		min_product = total_product;
		// set backpointer for reconstruction if needed
		// create a motif for hairpin loop
		auto motif_root = std::make_shared<decompose::LoopNode>();
		motif_root->type = 'p';
		auto node_copy = std::make_shared<decompose::LoopNode>(*node);
		node_copy->parent = motif_root;
		// no children to initialize
		motif_root->children_gm = std::vector<decompose::LoopNode::Ptr>{};
		motif_root->children_gm->push_back(node_copy);
		// motif_root->children.push_back(node);
		node_copy->child_id = 0;
		decompose::Mptr hairpin_motif = std::make_shared<decompose::Motif>(motif_root, 1);
		node->best_motif = hairpin_motif;
	}
	node->min_product_prob = min_product;
	return min_product;
}


//backtrack to get the motifs used in the min product decomposition
std::vector<decompose::Mptr> get_min_product_motifs(const decompose::LoopNode::Ptr &node) {
	std::vector<decompose::Mptr> motifs;
	if (!node->best_motif) {
		return motifs;
	}
	motifs.push_back(node->best_motif);
	// get leaf loops of best_motif->root
	const std::vector<decompose::LoopNode::Ptr> &children_to_use = node->best_motif->get_leaf_loops();
	// for each child, get corresponding node in original structure and recursive call
	for (const decompose::LoopNode::Ptr &child : children_to_use) {
		// get corresponding child in original structure
		auto parent = child->parent.lock();
		if (!parent || !child->child_id) {
			continue;
		}
		std::size_t child_index = child->child_id.value();
		if (child_index >= parent->children.size()) {
			continue;
		}
		auto node_in_original = parent->children[child_index];
		// recursive call
		std::vector<decompose::Mptr> child_motifs = get_min_product_motifs(node_in_original);
		motifs.insert(motifs.end(), child_motifs.begin(), child_motifs.end());
	}
	return motifs;
}

std::unordered_map<std::string, double> read_motif_prob_libs() {
	// lib files: data/short14_noex.txt.uniq.bestprob, short14_wiex.txt.bestprob, long_short.txt.bestprob, shortlong_reverse.txt.bestprob
	std::vector<std::string> lib_files = {
		"data/motif_libs/short14_noex.txt.uniq.bestprob",
		"data/motif_libs/short14_wiex.txt.bestprob",
		"data/motif_libs/long_short.txt.bestprob",
		"data/motif_libs/shortlong_reverse.txt.bestprob",
		"data/motif_libs/lib_undesignable.mprob",
		"data/motif_libs/lib_undesignable.ex.txt.mprob",
		"data/motif_libs/motif_cands.mr2.mprob",
		"data/motif_libs/eterna119_20260117192741.csv.motifs.mr.mprob",
		"data/motif_libs/eterna119_20260117192741.csv.motifs.mr2.mprob",
		"data/motif_libs/archiveii_normal400_20260117222631.csv.motifs.txt.mr2.mprob"
	};
	
	std::unordered_map<std::string, double> motif2prob;
	for (const auto &file : lib_files) {
		auto motif_probs = read_motif_prob_max(file);
		for (const auto &pair : motif_probs) {
			auto it = motif2prob.find(pair.first);
			if (it == motif2prob.end() || pair.second < it->second) {
				motif2prob[pair.first] = pair.second;
			}
		}
		std::cout << "Loaded " << motif_probs.size() << " motif probabilities from " << file << "\n";
	}
	std::cout << "Total motif probabilities loaded: " << motif2prob.size() << "\n";
	// print first 5 motif probabilities
	int n_example = 10;
	std::cout << "Example motif probabilities:\n";
	int count = 0;
	for (const auto &pair : motif2prob) {
		std::cout << "Motif: " << pair.first << ", Probability: " << pair.second << "\n";
		count += 1;
		if (count >= n_example) {
			break;
		}
	}
	return motif2prob;
}


double read_motif_time_cost(const std::string &file){
	// get the last field of each line, fields separated by comma ,
	std::ifstream infile(file);
	std::string line;
	double total_time = 0.0;
	int count = 0;
	while (std::getline(infile, line)) {
		std::stringstream ss(line);
		std::string field;
		std::vector<std::string> fields;
		while (std::getline(ss, field, ',')) {
			fields.push_back(field);
		}
		if (!fields.empty()) {
			try {
				double time = std::stod(fields.back());
				total_time += time;
				count += 1;
			} catch (const std::invalid_argument &e) {
				std::cerr << "Invalid time value in file " << file << ": " << fields.back() << "\n";
			}
		}
	}
	// print total time cost and average time cost
	if (count > 0) {
		std::cout << "Processed " << count << " time values from file " << file << "\n";
		std::cout << "Total time cost: " << total_time << ", Average time cost: " << total_time / count << "\n";
		return total_time / count;
	} else {
		std::cout << "No valid time values found in file " << file << "\n";
		return 0.0;
	}
}

std::unordered_map<std::string, double> read_motif_prob_libs_exact() {
	// lib files: data/short14_noex.txt.uniq.bestprob, short14_wiex.txt.bestprob, long_short.txt.bestprob, shortlong_reverse.txt.bestprob
	std::vector<std::string> lib_files = {
		"data/short14_noex.txt.uniq.bestprob",
		"data/short14_wiex.txt.bestprob",
		"data/long_short.txt.bestprob",
		"data/shortlong_reverse.txt.bestprob"
	};
	
	std::unordered_map<std::string, double> motif2prob;
	for (const auto &file : lib_files) {
		auto motif_probs = read_motif_prob_max(file);
		// insert if not already present or if the new probability is lower
		for (const auto &pair : motif_probs) {
			auto it = motif2prob.find(pair.first);
			if (it == motif2prob.end() || pair.second < it->second) {
				motif2prob[pair.first] = pair.second;
			}
		}
		std::cout << "Loaded " << motif_probs.size() << " motif probabilities from " << file << "\n";
	}
	std::cout << "Total motif probabilities loaded: " << motif2prob.size() << "\n";
	// print first 5 motif probabilities
	int n_example = 10;
	std::cout << "Example motif probabilities:\n";
	int count = 0;
	for (const auto &pair : motif2prob) {
		std::cout << "Motif: " << pair.first << ", Probability: " << pair.second << "\n";
		count += 1;
		if (count >= n_example) {
			break;
		}
	}
	return motif2prob;
}


std::unordered_map<std::string, double> read_motif_prob_libs_exact_simple() {
	// lib files: data/short14_noex.txt.uniq.bestprob, short14_wiex.txt.bestprob, long_short.txt.bestprob, shortlong_reverse.txt.bestprob
	std::vector<std::string> lib_files = {
		"data/short14_noex.txt.uniq.bestprob"
	};
	
	std::unordered_map<std::string, double> motif2prob;
	for (const auto &file : lib_files) {
		auto motif_probs = read_motif_prob_max(file);
		// insert if not already present or if the new probability is lower
		for (const auto &pair : motif_probs) {
			auto it = motif2prob.find(pair.first);
			if (it == motif2prob.end() || pair.second < it->second) {
				motif2prob[pair.first] = pair.second;
			}
		}
		std::cout << "Loaded " << motif_probs.size() << " motif probabilities from " << file << "\n";
	}
	std::cout << "Total motif probabilities loaded: " << motif2prob.size() << "\n";
	// print first 5 motif probabilities
	int n_example = 10;
	std::cout << "Example motif probabilities:\n";
	int count = 0;
	for (const auto &pair : motif2prob) {
		std::cout << "Motif: " << pair.first << ", Probability: " << pair.second << "\n";
		count += 1;
		if (count >= n_example) {
			break;
		}
	}
	return motif2prob;
}


void min_product(const std::string &structure, const std::unordered_map<std::string, double> &motif2prob) {
	// to be implemented
	decompose::LoopNode::Ptr root = decompose::LoopNode::string_to_node(structure);
	size_t count = count_decompositions(root);
	std::cout << "Total number of decompositions: " << count << "\n";
	// compute min product decomposition
	std::cout << "\nComputing minimum product decomposition:\n";
	double min_product_log = min_product_decomposition_log(root, motif2prob);
	double min_product_real = std::exp(min_product_log);
	std::cout << "Minimum product of motif probabilities for the entire structure: " << min_product_real << "\n";
	std::vector<decompose::Mptr> min_product_motifs = get_min_product_motifs(root);
	std::cout << "Motifs used in minimum product decomposition:\n";
	for (std::size_t i = 0; i < min_product_motifs.size(); ++i) {
		const decompose::Mptr &motif = min_product_motifs[i];
		std::cout << "Motif " << i + 1 << ", dot-bracket=" << motif->root->to_dot_bracket_gm()
				  << ", length=" << motif->root->length_gm() << '\n';
		std::vector<decompose::LoopNode::PairIndex> boundaries = motif->get_boundary_pairs();
		// print as a list of tuples in python format the same line
		std::cout << "[";
		for (std::size_t j = 0; j < boundaries.size(); ++j) {
			const decompose::LoopNode::PairIndex &p = boundaries[j];
			std::cout << "(" << p.open << ", " << p.close << ")";
			if (j < boundaries.size() - 1) {
				std::cout << ", ";
			}
		}
		std::cout << "]\n";
		// if motif in motif2prob, print its probability
		std::string motif_db = motif->root->to_dot_bracket_gm();
		auto it = motif2prob.find(motif_db);
		if (it != motif2prob.end()) {
			std::cout << "  Probability: " << it->second << "\n";
		}
	}
}


// decompose a dot-bracket string into motifs
std::vector<std::string> decompose_structure(const std::string &structure) {
	std::vector<std::string> decomposed_motifs;
	try {
		auto root = decompose::LoopNode::string_to_node(structure);
		std::queue<decompose::LoopNode::Ptr> q;
		q.push(root);
		// bfs to generate all motifs from each node
		while (!q.empty()) {
			decompose::LoopNode::Ptr current_node = q.front();
			q.pop();
			// std::cout << "----------------------------------------\n";
			// generate motifs from current_node
			std::vector<decompose::Mptr> node_motifs = current_node->generate_motifs();
			for (std::size_t i = 0; i < node_motifs.size(); ++i) {
				const decompose::Mptr &motif = node_motifs[i];
				std::string motif_db = motif->root->to_dot_bracket_gm();
				decomposed_motifs.push_back(motif_db);
			}
			// enqueue children
			for (const decompose::LoopNode::Ptr &child : current_node->children) {
				q.push(child);
			}
		}
	} catch (const std::exception &ex) {
		std::cerr << "Error in decomposing structure: " << ex.what() << std::endl;
	}
	return decomposed_motifs;
}


// compile: make lineardecompose
int main(int argc, char **argv) {
	if (argc < 2) {
		// std::cerr << "Usage: ./bin/lineardecompose minp or count"  << std::endl;
	}else{
		const std::string command = argv[1];
		std::unordered_map<std::string, double> motif2prob = read_motif_prob_libs();
		// minimum of products of motif probabilities
		if (command == "minp"){
			// read from stdin
			std::string structure;
			std::cout << "Enter dot-bracket structure: \n";
			while (std::getline(std::cin, structure)) {
				if (structure.empty()) {
					break;
				}
				// skip structures with sharp turns "()" or "(.)" or "(..)"
				if (structure.find("()") != std::string::npos || structure.find("(.)") != std::string::npos || structure.find("(..)") != std::string::npos) {
					std::cout << "Skipping structure with sharp turns: " << structure << "\n";
					std::cout << "Enter dot-bracket structure (or empty line to quit): ";
					continue;
				}
				std::cout << "Processing structure: " << structure << "\n";
				// record time
				auto start = std::chrono::high_resolution_clock::now();
				min_product(structure, motif2prob);
				auto end = std::chrono::high_resolution_clock::now();
				std::chrono::duration<double> elapsed = end - start;
				std::cout << "Time taken: " << elapsed.count() << "\n";
				std::cout << "Enter dot-bracket structure (or empty line to quit): \n";
			}
			return 0;
		}
		// count decompositions
		if (command == "count"){
			// read from stdin
			std::string structure;
			std::cout << "Enter dot-bracket structure: \n";
			while (std::getline(std::cin, structure)) {
				if (structure.empty()) {
					break;
				}
				// skip structures with sharp turns "()" or "(.)" or "(..)"
				if (structure.find("()") != std::string::npos || structure.find("(.)") != std::string::npos || structure.find("(..)") != std::string::npos) {
					std::cout << "Skipping structure with sharp turns: " << structure << "\n";
					std::cout << "Enter dot-bracket structure (or empty line to quit): ";
					continue;
				}
				std::cout << "Processing structure: " << structure << "\n";
				decompose::LoopNode::Ptr root = decompose::LoopNode::string_to_node(structure);
				size_t count = count_decompositions(root);
				std::cout << "Total number of decompositions: " << count << "\n";
				std::cout << "Enter dot-bracket structure (or empty line to quit): \n";
			}
		}
	}
	return 0;
}
