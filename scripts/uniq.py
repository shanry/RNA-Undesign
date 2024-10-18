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

from parser import get_mplotstr2


above_str = "midway, above=-5pt"
below_str = "midway, below_strow=-5pt"


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
		tree.parent = self
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
	for i, line in enumerate(open(path)):
		js = json.loads(line)
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
	print("struct. uniq:", len(id_uniqs), "\tmotif total:", total, "\tmotif uniq (min):", f"{sum(map(len, uniqs.values()))} ({sum(map(len, min_uniqs.values()))})")
	print("lengths:", sorted(lengths))
	print(sorted(list(id_uniqs)))
	filename = path + '.uniq'
	with open(filename, 'w') as f:
		for line in lines_uniqs:
			f.write(line)
	return uniqs


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
			print('plotstr:', get_mplotstr2(motif2js[tree][0]))
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


if __name__ == "__main__":
	uniqs = defaultdict(list) # loop-signature -> [motifs]
	path = sys.argv[1]
	print('path:', path)
	alg = sys.argv[2] if len(sys.argv) > 2 else 'dedup_lines'
	print('alg:', alg)
	if alg == 'dedup_lines':
		uniqs = dedup_lines(path)
		# uniqs = dedup(path)
		count_occurs(uniqs)
	elif alg == 'dotpr':
		gen_dotprnths(path)