#include <iostream>
#include <memory>
#include <stack>
#include <set>

#include "uniq.h"

struct TreeNode {
    std::string type;
    std::vector<int> unpaired;
    std::vector<std::shared_ptr<TreeNode>> children;
    TreeNode *parent = nullptr;
    int child_id = -1;
    int open_idx = -1;
    int close_idx = -1;

    std::string to_dotbracket() const {
        std::string result;
        if (type == "53") {
            result += "5";
            if (!children.empty()) {
                result += children[0]->to_dotbracket();
            }
            result += "3";
        } else if (type == "p") {
            if (child_id == -1) {
                if (!children.empty()) {
                    result = children[0]->to_dotbracket();
                }
            } else {
                result = "(*)";
            }
        } else if (type == "H") {
            result.push_back('(');
            result.append(unpaired.empty() ? 0 : unpaired[0], '.');
            result.push_back(')');
        } else {
            if (type != "E") {
                result.push_back('(');
            }
            if (!unpaired.empty()) {
                result.append(unpaired[0], '.');
            }
            for (std::size_t i = 0; i < children.size(); ++i) {
                result += children[i]->to_dotbracket();
                if (i + 1 < unpaired.size()) {
                    result.append(unpaired[i + 1], '.');
                }
            }
            if (type != "E") {
                result.push_back(')');
            }
        }
        return result;
    }

    std::shared_ptr<TreeNode> clone_subtree() const {
        auto node = std::make_shared<TreeNode>();
        node->type = type;
        node->unpaired = unpaired;
        node->child_id = child_id;
        node->open_idx = open_idx;
        node->close_idx = close_idx;
        for (const auto &child : children) {
            auto copied = child->clone_subtree();
            copied->parent = node.get();
            node->children.push_back(std::move(copied));
        }
        return node;
    }

    std::shared_ptr<TreeNode> make_tree(int child_index) const {
        auto tree = std::make_shared<TreeNode>();
        tree->type = type;
        tree->unpaired = unpaired;
        tree->child_id = child_index;

        if (parent != nullptr) {
            auto parent_as_child = parent->make_tree(child_id);

            auto rotate_children = [&, this](std::size_t start, std::size_t end) {
                for (std::size_t i = start; i < end; ++i) {
                    auto copied = children[i]->clone_subtree();
                    copied->parent = tree.get();
                    tree->children.push_back(std::move(copied));
                }
            };

            std::size_t start = static_cast<std::size_t>(child_index + 1);
            if (child_index < 0) {
                start = 0;
            }
            rotate_children(start, children.size());

            parent_as_child->parent = tree.get();
            tree->children.push_back(std::move(parent_as_child));

            std::size_t end = child_index < 0 ? 0 : static_cast<std::size_t>(child_index);
            rotate_children(0, end);

            std::vector<int> rotated_unpaired;
            std::size_t up_start = static_cast<std::size_t>(child_index + 1);
            if (child_index < 0) {
                up_start = 0;
            }
            for (std::size_t i = up_start; i < unpaired.size(); ++i) {
                rotated_unpaired.push_back(unpaired[i]);
            }
            std::size_t limit = child_index < 0 ? 0 : static_cast<std::size_t>(child_index + 1);
            for (std::size_t i = 0; i < limit && i < unpaired.size(); ++i) {
                rotated_unpaired.push_back(unpaired[i]);
            }
            tree->unpaired = rotated_unpaired;
        } else {
            for (const auto &child : children) {
                auto copied = child->clone_subtree();
                copied->parent = tree.get();
                tree->children.push_back(std::move(copied));
            }
        }
        return tree;
    }

    void rotated(int depth, std::vector<std::shared_ptr<TreeNode>> &out) const {
        if (depth == 0 && type == "53") {
            return;
        } else if (depth == 0 && type == "p") {
            for (const auto &child : children) {
                child->rotated(depth + 1, out);
            }
        } else if (depth > 0 && type == "p") {
            out.push_back(make_tree(-1));
        } else {
            for (const auto &child : children) {
                child->rotated(depth + 1, out);
            }
        }
    }
};

std::shared_ptr<TreeNode> parse_dotbracket(const std::string &s) {
    auto root = std::make_shared<TreeNode>();
    root->type = "p";
    root->child_id = -1;

    struct Frame {
        std::shared_ptr<TreeNode> node;
        int open_idx;
    };

    std::vector<Frame> stack;
    stack.push_back({root, -1});

    for (std::size_t i = 0; i < s.size(); ++i) {
        char c = s[i];
        if (c == '(') {
            auto n = std::make_shared<TreeNode>();
            n->parent = stack.back().node.get();
            n->open_idx = static_cast<int>(i);
            stack.push_back({n, static_cast<int>(i)});
        } else if (c == ')') {
            if (stack.size() <= 1) {
                throw std::runtime_error("Unbalanced parentheses in dot-bracket string.");
            }
            auto frame = stack.back();
            stack.pop_back();
            int open = frame.open_idx;
            int close = static_cast<int>(i);
            frame.node->open_idx = open;
            frame.node->close_idx = close;
            auto &parent = stack.back().node;
            frame.node->child_id = static_cast<int>(parent->children.size());
            frame.node->parent = parent.get();

            if (close - open == 2) {
                frame.node->type = "p";
            } else if (frame.node->children.empty()) {
                frame.node->type = "H";
                frame.node->unpaired.push_back(close - open - 1);
            } else if (frame.node->children.size() == 1) {
                frame.node->type = "I";
                const auto &child = frame.node->children[0];
                frame.node->unpaired.push_back(child->open_idx - open - 1);
                frame.node->unpaired.push_back(close - child->close_idx - 1);
            } else {
                frame.node->type = "M";
                const auto &first_child = frame.node->children.front();
                frame.node->unpaired.push_back(first_child->open_idx - open - 1);
                for (std::size_t j = 0; j + 1 < frame.node->children.size(); ++j) {
                    const auto &left = frame.node->children[j];
                    const auto &right = frame.node->children[j + 1];
                    frame.node->unpaired.push_back(right->open_idx - left->close_idx - 1);
                }
                const auto &last_child = frame.node->children.back();
                frame.node->unpaired.push_back(close - last_child->close_idx - 1);
            }
            parent->children.push_back(std::move(frame.node));
        }
    }
    if (stack.size() != 1) {
        throw std::runtime_error("Unbalanced parentheses in dot-bracket string.");
    }
    return root;
}

// a function accept a motif in dot-bracket notation, generate all rotated variants in dot-bracket notation, and return them
std::vector<std::string> generate_rotated_variants(const std::string &motif) {
    // if motif starts with '5' and ends with '3', there is no rotation needed
    if (motif.size() >= 2 && motif.front() == '5' && motif.back() == '3') {
        return {motif};
    }
    std::vector<std::string> variants;
    auto node = parse_dotbracket(motif);
    variants.push_back(node->to_dotbracket());
    std::vector<std::shared_ptr<TreeNode>> rotated_nodes;
    node->rotated(0, rotated_nodes);
    for (const auto &rn : rotated_nodes) {
        variants.push_back(rn->to_dotbracket());
    }
    return variants;
}

// int main() {
//     std::ios::sync_with_stdio(false);
//     std::cin.tie(nullptr);

//     std::set<std::string> total_variants;

//     std::string line;
//     while (std::getline(std::cin, line)) {
//         std::string trimmed;
//         for (char c : line) {
//             if (c != '\n' && c != '\r' && c != ' ') {
//                 trimmed.push_back(c);
//             }
//         }
//         if (trimmed.empty()) {
//             continue;
//         }
//         if (trimmed == ".." || trimmed == "..." || trimmed == "....") {
//             continue;
//         }

//         try {
//             auto node = parse_dotbracket(trimmed);
//             std::vector<std::shared_ptr<TreeNode>> variants;
//             variants.push_back(node);
//             node->rotated(0, variants);

//             std::cout << "total variants: " << variants.size() << "\n";
//             for (const auto &v : variants) {
//                 std::cout << v->to_dotbracket() << "\n";
//                 total_variants.insert(v->to_dotbracket());
//             }
//         } catch (const std::exception &ex) {
//             std::cerr << "Error parsing line: " << ex.what() << "\n";
//         }
//     }
//     std::cout << "total unique variants: " << total_variants.size() << "\n";
//     return 0;
// }
