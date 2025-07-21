#!/usr/bin/env python3
"""
This script is used for deduplicating and analyzing RNA motifs.

The script reads a file containing JSON-encoded RNA motifs and performs deduplication based on the loop signature of each motif. The loop signature is determined by counting the occurrences of different loop types (e.g., hairpin, stack, bulge) in the motif. Motifs with the same loop signature are considered duplicates.

The deduplication process involves creating a tree structure representation of each motif and comparing the tree structures of different motifs to identify duplicates. The script also counts the occurrences of each motif and categorizes them by family.

The deduplicated motifs and their occurrence counts are printed to the console.

Usage:
	python uniq.py <input_file>

Args:
	input_file (str): Path to the file containing JSON-encoded RNA motifs.

Example:
	python uniq.py motifs.json
"""

import sys
import json
from collections import defaultdict

import numpy as np
import pandas as pd
import math
from matplotlib import pyplot as plt

# from parser import get_mplotstr2

oldfamily2newfamily = {'16s': '16S rRNA', 
                       '23s': '23S rRNA', 
                       '5s': '5S rRNA', 
                       'RNaseP': 'RNaseP', 
                       'grp1': 'Group I Intron ', 
                       'grp2': 'Group II Intron', 
                       'srp': 'SRP', 
                       'tRNA': 'tRNA', 
                       'telomerase': 'telomerase', 
                       'tmRNA': 'tmRNA'}
new2old = {new: old for old, new in oldfamily2newfamily.items()}

class Node:
	__slots__ = "type", "children", "unpaired_bases", "parent", "child_id", "ids"

	def __init__(self, jsontree={}, parent=None, child_id=-1):
		self.children = []
		self.unpaired_bases = []
		self.parent = parent
		self.child_id = child_id
		if jsontree is None or jsontree == {}:
			self.type = 'p'
		elif jsontree != {}:
			if 'root' in jsontree: # very root
				self.ids = [jsontree['id']]
				if jsontree["root"]["child_id"] == -1:
					self.type = "53"
					self.child_id = -1
				else:
					self.type = "p"
				self.children = [Node(jsontree['root'], parent=self, child_id=0)]
			else: # non-root
				self.type = jsontree['type']
				self.unpaired_bases = jsontree.get('loops', [])
				if 'children' in jsontree:
					for i, child in enumerate(jsontree['children']):
						subtree = Node(child, parent=self, child_id=i)
						self.children.append(subtree)
				else: # leaf node
					if self.type == "H":
						pass
					elif self.type == "M": # leaf M; (deg-1) p's
						for i in range(len(jsontree['loops'])-1):
							self.children.append(Node(None, parent=self, child_id=i)) # leaf p
					elif self.type in ["I", "S", "B"]: # single p
						self.children.append(Node(None, parent=self, child_id=0))
					elif self.type == "E":
						assert False, "E loop; TODO"
		# if self.child_id == -1:
		# 	assert len(self.children) == 1, str(jsontree)


	def pp(self, dep=0):
		print(" |" * dep, self.type, len(self.children), self.unpaired_bases, self.child_id)
		for child in self.children:
			child.pp(dep+1)


	def __str__(self):
		return self.to_dotbracket()
		# return "%s %s (%s)" % (self.type, self.unpaired_bases,
		# 					   ", ".join(map(str, self.children)))
	__repr__ = __str__

	def to_bfs(self):
		return "%s %s (%s)" % (self.type, self.unpaired_bases,
							   ", ".join(child.to_bfs() for child in self.children))

	def make_tree(self, child_id):
		tree = Node(child_id=child_id)
		tree.type = self.type
		# tree.parent = self
		if child_id != -1:
			tree.parent = self.children[child_id]
		tree.unpaired_bases = self.unpaired_bases[:]
		tree.child_id = child_id
		if self.parent is not None:
			# rotate children/parent
			parent_as_child = self.parent.make_tree(self.child_id)			
			tree.children = self.children[child_id+1:] \
			           + [parent_as_child] \
			           + self.children[:child_id]
			# rotate unpaired segments
			tree.unpaired_bases = self.unpaired_bases[child_id+1:] + self.unpaired_bases[:child_id+1]
		return tree

	def rotated(self, dep=0): # returns a lazylist of rotated trees from each leaf p node
		if dep == 0 and self.type == "53":
			return
		elif dep == 0 and self.type == "p":
			for child in self.children:
				for tree in child.rotated(dep+1):
					yield tree
		elif dep > 0 and self.type == "p":
			# make a new tree from this node
			yield self.make_tree(-1)
		else:
			for child in self.children:
				for tree in child.rotated(dep+1):
					yield tree

	# convert to dotbracket notation
	def to_dotbracket(self): 
		result = ""
		if self.type == "53":
			result += "5"
			result += self.children[0].to_dotbracket()
			result += "3"
		elif self.type == "p":
			if self.child_id == -1:
				assert len(self.children) > 0
				result = self.children[0].to_dotbracket()
			else:
				result = "(*)"
		elif self.type == "H":
			result = "(" + "." * self.unpaired_bases[0] + ")"
		else:
			result = ""
			if self.type != "E":
				result += "("
			result += "." * self.unpaired_bases[0]
			for i, child in enumerate(self.children):
				result += child.to_dotbracket()
				assert i + 1 < len(self.unpaired_bases)
				result += "." * self.unpaired_bases[i+1]
			if self.type != "E":
				result += ")"
		return result


class PairNode:

	def __init__(self, pair, children=None, type='unknown', child_id=-1, parent=None):
		self.pair = pair
		self.type = type  # no external loop type
		self.unpaired_bases = []
		if children is None:
			children = []  # Initialize to an empty list if not provided
		self.children = children
		self.child_id = child_id
		self.parent = parent
	
	@classmethod
	def string2node(cls, s): # dotbracket to Node, e.g., ((*).(...(*)).), ((*).(.((*))).)
		root = PairNode((-1, -1), type='p')
		stack = []
		stack.append(root)
		for i, c in enumerate(s):
			# print(i, c)
			if c == '(':
				stack.append(PairNode((i, -1)))
			elif c == ')':
				# print('len(stack):', len(stack))
				pair_node = stack.pop()
				pair_node.pair = (pair_node.pair[0], i)
				# print(pair_node.pair, '->', stack[-1].pair)
				# print(stack[-1].children)
				stack[-1].children.append(pair_node)
				pair_node.parent = stack[-1]
				pair_node.child_id = len(stack[-1].children) - 1
				assert pair_node.pair[0] >= 0
				if pair_node.pair[1] - pair_node.pair[0] == 2: # leaf node
					pair_node.type = 'p'
				elif len(pair_node.children) == 0: # hairpin
					pair_node.type = 'H'
					pair_node.unpaired_bases = [pair_node.pair[1] - pair_node.pair[0] - 1]
				elif len(pair_node.children) == 1: # internal loop
					pair_node.type = 'I'
					pair_node.unpaired_bases = [pair_node.children[0].pair[0] - pair_node.pair[0] - 1, pair_node.pair[1] - pair_node.children[0].pair[1] - 1]
				else: # multi-loop
					pair_node.type = 'M'
					pair_node.unpaired_bases = [pair_node.children[0].pair[0] - pair_node.pair[0] - 1]
					for j in range(len(pair_node.children) - 1):
						pair_node.unpaired_bases.append(pair_node.children[j+1].pair[0] - pair_node.children[j].pair[1] - 1)
					pair_node.unpaired_bases.append(pair_node.pair[1] - pair_node.children[-1].pair[1] - 1)
				# print(pair_node.pair, len(pair_node.children), pair_node.type, pair_node.unpaired_bases)
		assert len(stack) == 1
		return root
	
	def to_dotbracket(self):
		result = ""
		if self.type == 'p':
			if self.child_id == -1:
				assert len(self.children) == 1, str(len(self.children))
				result = self.children[0].to_dotbracket()
			else:
				result = "(*)"
		elif self.type == 'H':
			result = "(" + "." * self.unpaired_bases[0] + ")"
		else:
			result = ""
			if self.type != "E":
				result += "("
			result += "." * self.unpaired_bases[0]
			for i, child in enumerate(self.children):
				result += child.to_dotbracket()
				assert i + 1 < len(self.unpaired_bases)
				result += "." * self.unpaired_bases[i+1]
			if self.type != "E":
				result += ")"
		return result
	

	def __str__(self):
		return self.to_dotbracket()
	

	def make_tree(self, child_id):
		tree = Node(child_id=child_id)
		tree.type = self.type
		if child_id != -1:
			tree.parent = self.children[child_id]
		tree.unpaired_bases = self.unpaired_bases[:]
		tree.child_id = child_id
		if self.parent is not None:
			# rotate children/parent
			parent_as_child = self.parent.make_tree(self.child_id)			
			tree.children = self.children[child_id+1:] \
			           + [parent_as_child] \
			           + self.children[:child_id]
			# rotate unpaired segments
			tree.unpaired_bases = self.unpaired_bases[child_id+1:] + self.unpaired_bases[:child_id+1]
		return tree

	def rotated(self, dep=0): # returns a lazylist of rotated trees from each leaf p node
		if dep == 0 and self.type == "53":
			return
		elif dep == 0 and self.type == "p":
			for child in self.children:
				for tree in child.rotated(dep+1):
					yield tree
		elif dep > 0 and self.type == "p":
			# make a new tree from this node
			yield self.make_tree(-1)
		else:
			for child in self.children:
				for tree in child.rotated(dep+1):
					yield tree


def get_ref2structure():
	with open('data/ref2structure.json') as f:
		ref2structure = json.load(f)
	return ref2structure


def loop_stats(tree, cache=None): # returns {B:1, M:3, ..}
	if cache is None:
		cache = defaultdict(int)
	cache[tree['type']] += 1
	if 'children' in tree:
		for sub in tree['children']:
			if sub is not None:
				loop_stats(sub, cache)
	if 'loops' in tree:
		cache["unp"] += sum(tree['loops'])
	return cache


def get_length(js):
	bpair = js["bpairs"]
	ystar = js["y_star"]
	length = bpair[0][1] - bpair[0][0] + 1 if bpair[0][0] != -1 else len(ystar)
	for pair in bpair[1:]:
		length -= pair[1] - pair[0] - 1
	return length


def get_lonepairs(js):
	bpair = js["bpairs"]
	ipair = js["ipairs"]
	idx_pairs = set()
	for pair in bpair:
		if pair[0] == -1:
			continue
		idx_pairs.add(pair[0])
		idx_pairs.add(pair[1])
	for pair in ipair:
		if pair[0] == -1:
			continue
		idx_pairs.add(pair[0])
		idx_pairs.add(pair[1])
	
	pairs_stronglylone = []  # strongly lonely pair does not have any neighbored positions that are also in idx_pairs
	pairs_weaklylone = []  # weakly lonely pair has either left or right position that doesn't have any neighbored positions also in idx_pairs
	for pair in ipair:
		if pair[0] == -1:
			continue
		count_lone = 0
		if pair[0] - 1 not in idx_pairs and pair[0] + 1 not in idx_pairs:
			count_lone += 1
		if pair[1] - 1 not in idx_pairs and pair[1] + 1 not in idx_pairs:
			count_lone += 1
		if count_lone == 2:
			pairs_stronglylone.append(pair)
		if count_lone == 1:
			pairs_weaklylone.append(pair)

	return pairs_stronglylone, pairs_weaklylone
		

def dedup(path):
	uniqs = defaultdict(list) # loop-signature -> [motifs]
	id_uniqs = set()
	for i, line in enumerate(open(path)):
		js = json.loads(line)
		id_uniqs.add(js['id'])
		# js_full = json.loads(line)
		ids, motif = js['motif']['id'], js['motif']['root']
		signature = str(sorted(loop_stats(motif).items())) # "{B:1, M:3, ..}"
		tree = Node(js['motif'])
		all_rotations = set([str(tree)])
		for newtree in tree.rotated(0):
			all_rotations.add(str(newtree))
		for othertree, otherjss in uniqs[signature]:
			if str(othertree) in all_rotations:
				otherjss.append(js)
				break
		else: # new uniq
			uniqs[signature].append((tree, [js]))
	total = 0
	for signature in uniqs:
		for tree, jss in uniqs[signature]:
			total += len(jss)
			for js in jss:
				print(json.dumps(js['motif']))
	print("motif total:", total, "motif uniq:", sum(map(len, uniqs.values())), "struct. uniq:", len(id_uniqs))
	print(sorted(list(id_uniqs)))
	return uniqs


def dedup_lines(path):
	uniqs = defaultdict(list) # loop-signature -> [motifs]
	min_uniqs = defaultdict(list) # loop-signature -> [motifs]
	id_uniqs = set()
	lines_uniqs = []
	lengths = []
	count_motif_min = 0
	for i, line in enumerate(open(path)):
		js = json.loads(line)
		if js['ismin']:
			count_motif_min += 1
		id_uniqs.add(js['id'])
		# js_full = json.loads(line)
		ids, motif = js['motif']['id'], js['motif']['root']
		signature = str(sorted(loop_stats(motif).items())) # "{B:1, M:3, ..}"
		tree = Node(js['motif'])
		all_trees_rotated = []
		all_trees_rotated.append(tree)
		all_rotations = set([str(tree)])
		# print(str(tree))
		# print(tree.to_bfs())
		for newtree in tree.rotated(0):
			all_trees_rotated.append(newtree)
			all_rotations.add(str(newtree))
		# if js['ismin'] is False:
		# 	print(i+1, str(tree))
		# for rt in all_trees_rotated:
		# 	print(rt.to_bfs())
		# for rt in all_trees_rotated:
		# 	print(str(rt))
		for othertree, otherjss in uniqs[signature]:
			if str(othertree) in all_rotations:
				otherjss.append(js)
				break
		else: # new uniq
			uniqs[signature].append((tree, [js]))
			# print(Node.to_string(tree))
			# lines_uniqs.append(line)
			# print(tree, ids, len(lines_uniqs), signature)
			if js['ismin']:
				min_uniqs[signature].append((tree, [js]))
				lengths.append(get_length(js))
				lines_uniqs.append(line)
	total = 0
	for signature in uniqs:
		for tree, jss in uniqs[signature]:
			total += len(jss)
	# 		for js in jss:
	# 			print(json.dumps(js['motif']))
	print("lengths:", sorted(lengths))
	print("size of lengths:", len(lengths))
	print("struct. uniq:", len(id_uniqs), "\tmotif total:", total, "\tminimal motif total:", count_motif_min, "\tmotif uniq (min):", f"{sum(map(len, uniqs.values()))} ({sum(map(len, min_uniqs.values()))})")
	print("average length:", np.mean(lengths))
	print("median  length:", np.median(lengths))
	print("max    length:", np.max(lengths))
	print("min    length:", np.min(lengths))
	# print(sorted(list(id_uniqs)))
	filename = path + '.uniq'
	with open(filename, 'w') as f:
		for line in lines_uniqs:
			f.write(line)
	return uniqs


def get_motif_statistics(path):
	loop_maps = []
	loop_nums = []
	lone_nums = []
	semilone_nums = []
	for i, line in enumerate(open(path)):
		js = json.loads(line)
		# js_full = json.loads(line)
		_, motif = js['motif']['id'], js['motif']['root']
		ipairs = js['ipairs']
		loop_maps.append(loop_stats(motif))
		num_loops = 0
		for loop_type, count in loop_maps[-1].items():
			if loop_type == 'unp':
				continue
			num_loops += count
		assert num_loops == len(ipairs) + 1, str(num_loops) + ' != ' + str(len(ipairs))  + '1'
		loop_nums.append(num_loops)
		pairs_lone, pairs_semilone = get_lonepairs(js)
		lone_nums.append(len(pairs_lone))
		semilone_nums.append(len(pairs_semilone))
	
	loops_count = defaultdict(int)
	for i, loop_map in enumerate(loop_maps):
		print(loop_map, loop_nums[i])
		num_loops = 0
		for loop_type, count in loop_map.items():
			if loop_type == 'unp':
				continue
			num_loops += count
		assert num_loops == loop_nums[i], str(len(loop_map)) + ' != ' + str(loop_nums[i])
		loops_count[num_loops] += 1
	print("loops count:", loops_count)

	lone_count = defaultdict(int)
	for i, lone_num in enumerate(lone_nums):
		lone_count[lone_num] += 1
	print("strongly isolated count:", lone_count)
	lone_percent = sum([lone_count[i] for i in lone_count if i > 0]) / len(lone_nums)
	print("strongly isolated percent:", f"{lone_percent:.2%}")

	semilone_count = defaultdict(int)
	for i, semilone_num in enumerate(semilone_nums):
		semilone_count[semilone_num] += 1
	print("weakly isolated count:", semilone_count)
	semilone_percent = sum([semilone_count[i] for i in semilone_count if i > 0]) / len(semilone_nums)
	print(f"weakly isolated percent: {semilone_percent:.2%}")

	isolated_ratio = lone_percent + semilone_percent
	print(f"isolated ratio: {isolated_ratio:.2%}")


def count_occurs(uniqs):
	motif2family = defaultdict(lambda: defaultdict(int))
	motif2js = defaultdict(list)
	for signature in uniqs:
		for tree, jss in uniqs[signature]:
			for js in jss:
				family = js['id'].split('_')[0]
				motif2family[tree][family] += 1
				motif2js[tree].append(js)
	for tree in motif2family:
		if len(motif2family[tree]) > 1:
			print('tree:', tree)
			print('occurs dict:', {family: motif2family[tree][family] for family in motif2family[tree]})
			print('family count:', len(motif2family[tree]))
			print('occurs count:', len(motif2js[tree]))
			# print('plotstr:', get_mplotstr2(motif2js[tree][0]))
			print(motif2js[tree][0])
			print()


def gen_dotprnths(path):
	uniq_trees = []
	strs = []
	uniq_strs = set()
	str2id = {}
	for i, line in enumerate(open(path)):
		js = json.loads(line)
		tree = Node(js['motif'])
		if str(tree) not in uniq_strs: # if original tree is not in uniq_strs
			uniq_trees.append(tree)
			uniq_strs.add(str(tree))
			strs.append(str(tree))
			str2id[str(tree)] = i+1
			for newtree in tree.rotated(0):
				newstr = str(newtree)
				assert newstr not in uniq_strs or str2id[newstr] == i+1, str(i+1) + ':\t' + str(str2id[newstr]) + '\n' + line # then otated trees should not be in uniq_strs
				if newstr not in uniq_strs:
					uniq_strs.add(newstr)
					strs.append(newstr)
					str2id[newstr] = i+1
	print('total trees:', len(uniq_trees), 'uniq strs:', len(strs))
	with open(path + '.dotprnths', 'w') as f:
		for s in strs:
			f.write(s + '\n')
	print('output:', path + '.dotprnths')


def conditioned_count(path, func_condition):
	# ref2structure = get_ref2structure()  # for archiveII dataset
	ref2structure = dict()
	df = pd.read_csv('data/rnastralign_filtered.csv')
	for i, row in df.iterrows():
		id = row['id']
		structure = row['secondary_structure']
		ref2structure[id] = structure
	uniqs = defaultdict(list) # loop-signature -> [motifs]
	min_uniqs = defaultdict(list) # loop-signature -> [motifs]
	id_uniqs = set()
	structure_uniqs = set()
	lines_uniqs = []
	lengths = []
	count_motif_min = 0
	for i, line in enumerate(open(path)):
		js = json.loads(line)
		if not func_condition(js):
			continue
		if js['ismin']:
			count_motif_min += 1
		id_uniqs.add(js['id'])
		structure_uniqs.add(ref2structure[js['id']])
		# js_full = json.loads(line)
		ids, motif = js['motif']['id'], js['motif']['root']
		signature = str(sorted(loop_stats(motif).items())) # "{B:1, M:3, ..}"
		tree = Node(js['motif'])
		all_trees_rotated = []
		all_trees_rotated.append(tree)
		all_rotations = set([str(tree)])
		# print(str(tree))
		# print(tree.to_bfs())
		for newtree in tree.rotated(0):
			all_trees_rotated.append(newtree)
			all_rotations.add(str(newtree))
		# for rt in all_trees_rotated:
		# 	print(rt.to_bfs())
		# for rt in all_trees_rotated:
		# 	print(str(rt))
		for othertree, otherjss in uniqs[signature]:
			if str(othertree) in all_rotations:
				otherjss.append(js)
				break
		else: # new uniq
			uniqs[signature].append((tree, [js]))
			# print(Node.to_string(tree))
			lines_uniqs.append(line)
			# print(tree, ids, len(lines_uniqs), signature)
			if js['ismin']:
				min_uniqs[signature].append((tree, [js]))
			lengths.append(get_length(js))
	total = 0
	for signature in uniqs:
		for tree, jss in uniqs[signature]:
			total += len(jss)
	# 		for js in jss:
	# 			print(json.dumps(js['motif']))
	count_motif = total
	count_motif_uniq = sum(map(len, uniqs.values()))
	count_motif_uniq_min = sum(map(len, min_uniqs.values()))
	count_structure_id = len(id_uniqs)
	count_structure = len(structure_uniqs)
	print("motif total:", total, "\tminimal motif total:", count_motif_min, "\tmotif uniq (min):", f"{sum(map(len, uniqs.values()))} ({sum(map(len, min_uniqs.values()))})")
	print("id uniq:", len(id_uniqs), "\tstruct. uniq:", len(structure_uniqs))
	print("lengths:", sorted(lengths))
	# print(sorted(list(id_uniqs)))
	print('-----------------------------------')
	return count_motif_min, count_motif_uniq, count_motif_uniq_min, count_structure_id, count_structure


if __name__ == "__main__":
	uniqs = defaultdict(list) # loop-signature -> [motifs]
	path = sys.argv[1]
	print('path:', path)
	alg = sys.argv[2] if len(sys.argv) > 2 else 'dedup_lines'
	print('alg:', alg)
	if alg == 'dedup_lines':
		uniqs = dedup_lines(path)
		# uniqs = dedup(path)
		# count_occurs(uniqs)
	elif alg == 'stats':
		get_motif_statistics(path)
	elif alg == 'dotpr':
		gen_dotprnths(path)
	elif alg == 'cond':
		# family_list = ['tRNA', '5S rRNA', 'SRP', 'RNaseP', 'tmRNA', 'Group I Intron ', 'telomerase', 'Group II Intron', '16S rRNA', '23S rRNA']  # for archiveII dataset
		family_list = [ 'tRNA', 'SRP', 'RNaseP', 'tmRNA', 'group_I_intron', 'telomerase', '16S_rRNA']
		# family_list = []
		data = []
		for newfamily in family_list:
			# family = new2old[newfamily]  # for archiveII dataset 
			family = newfamily
			print('family:', family)
			count_motif, count_motif_uniq, count_motif_uniq_min, count_structure_id, count_structure = conditioned_count(path, lambda js: family in js['id'])
			data.append([newfamily, count_motif, count_motif_uniq, count_motif_uniq_min, count_structure_id, count_structure])
		df = pd.DataFrame(data, columns=['Family', 'MotifTotal', 'MotifUniq', 'MotifUniqMin', 'StructureID', 'Structure'])
		summary = ['Total', df['MotifTotal'].sum(), df['MotifUniq'].sum(), df['MotifUniqMin'].sum(), df['StructureID'].sum(), df['Structure'].sum()]
		df.loc[len(df)] = summary
		print(df)
		df.to_csv(path + '.cond.csv', index=False)
	elif alg == 'mstr': # cat data/short14_undesignable_dg0.txt | ./scripts/uniq.py xx mstr
		uniq_trees = []
		all_strs = set()
		for line in sys.stdin:
			if line.strip() == '..' or line.strip() == '...' or line.strip() == '....':
				continue
			s = line.strip()
			node = PairNode.string2node(s)
			if str(node) not in all_strs:
				uniq_trees.append(node)
				all_strs.add(str(node))
				for rt in node.rotated(0):
					if str(rt) not in all_strs:
						all_strs.add(str(rt))
		print('total trees:', len(all_strs), 'uniq strs:', len(uniq_trees))
		# for tree in uniq_trees:
		# 	print(str(tree))