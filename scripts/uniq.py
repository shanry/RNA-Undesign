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
		if jsontree is None:
			self.type = 'p'
		elif jsontree != {}:
			if 'root' in jsontree: # very root
				self.ids = [jsontree['id']]
				self.type = "p" # top p node
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


	def pp(self, dep=0):
		print(" |" * dep, self.type, len(self.children), self.unpaired_bases, self.child_id)
		for child in self.children:
			child.pp(dep+1)


	def __str__(self):
		return "%s %s (%s)" % (self.type, self.unpaired_bases,
							   ", ".join(map(str, self.children)))
	__repr__ = __str__


	def make_tree(self, child_id):
		tree = Node()
		tree.type = self.type
		tree.parent = self
		tree.unpaired_bases = self.unpaired_bases[:]

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
		if dep ==0 and self.type == "E":
			return
		if dep > 0 and self.type == "p":
			# make a new tree from this node
			yield self.make_tree(self.child_id)
		for child in self.children:
			for tree in child.rotated(dep+1):
				yield tree


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

def check_eq(a, b):
	return True

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
	id_uniqs = set()
	lines_uniqs = []
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
			lines_uniqs.append(line)
	total = 0
	for signature in uniqs:
		for tree, jss in uniqs[signature]:
			total += len(jss)
	# 		for js in jss:
	# 			print(json.dumps(js['motif']))
	print("motif total:", total, "motif uniq:", sum(map(len, uniqs.values())), "struct. uniq:", len(id_uniqs))
	print(sorted(list(id_uniqs)))
	filename = path.replace('.txt', '.uniq.txt')
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


if __name__ == "__main__":
	uniqs = defaultdict(list) # loop-signature -> [motifs]
	path = sys.argv[1]
	print('path:', path)
	uniqs = dedup_lines(path)
	count_occurs(uniqs)