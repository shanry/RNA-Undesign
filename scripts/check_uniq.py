#!/usr/bin/env python3

import sys
import json
from collections import defaultdict

'''
cat motifs.json | ./check_uniq.py

where each line in the input is a json string like
{"id":"86","root":{"child_id":1,"children":[null,{"child_id":2,"loops":[3],"neighbors":[0],"type":"H"}],"loops":[1,1,1],"neighbors":[0,1,2],"type":"M"}}

output the list of unique motifs and counts.

e.g., check all 16s/23s:

cat raws/[12]* | grep ^{ | ./check_uniq.py

e.g., check all eterna:

cat raws/et22* | grep ^{ | ./check_uniq.py

'''

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
		if dep > 0 and self.type == "p":
			# make a new tree from this node
			yield self.make_tree(self.child_id)
		for child in self.children:
			for tree in child.rotated(dep+1):
				yield tree

	def tolatex(self, dep=0): # call with dep=1 for motif
		arity = len(self.children) + 1
		if dep == 0:
			s = "[.{\\tiny \it root} "
			s += "\\edge node[%s]{\\tiny %d};" % (above_str, self.unpaired_bases[0])
			s += "[.E "
		else:
			if arity == 1:
				s = " H " # hairpin
			elif arity == 2:
				if sum(self.unpaired_bases) == 0:
					looptype = "S" # stack
				elif min(self.unpaired_bases) == 0:
					looptype = "B" # bulge
				else:
					looptype = "I" # internal
				s = "[.%s " % looptype
			else:
				s = "[.M "
		for i, child in enumerate(self.children):
			s += "\\edge node[%s]{\\tiny %d};" % (above_str if i == 0 else below_str,
							     				  self.unpaired_bases[i+1])
			s += "[.{\\tiny \\it p} "
			s += "\\edge node[%s]{\\tiny %d};" % (above_str, child.unpaired_bases[0])

			s += child.tolatex(dep+1)
			s += "] " # p
		if arity >1: # NB?
			s += "] "
		if dep == 0:
			s += "] "
		return s		


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

if __name__ == "__main__":
	uniqs = defaultdict(list) # loop-signature -> [motifs]
	for i, line in enumerate(sys.stdin):
		js = json.loads(line)['motif']	
		ids, motif = js['id'], js['root']
		print(motif)
		signature = str(sorted(loop_stats(motif).items())) # "{B:1, M:3, ..}"
		print(signature)
		tree = Node(js)
		tree.pp()
		print(tree)
		all_rotations = set([str(tree)])
		for newtree in tree.rotated(0):
			print("rotated")
			newtree.pp()
			all_rotations.add(str(newtree))
		for othertree, otherjss in uniqs[signature]:
			if str(othertree) in all_rotations:
				print("found dup:", othertree, motif)
				otherjss.append(js)
				break
		else: # new uniq
			#motif['id'] = [motif['id']] # list of ids; for duplicates
			uniqs[signature].append((tree, [js]))
	for signature in uniqs:
		print(signature)
		for tree, jss in uniqs[signature]:
			print(tree)
			print(len(jss))
			for js in jss:
				print(" ", json.dumps(js))

		print()

	print("tot", i+1, "uniq", sum(map(len, uniqs.values())))




