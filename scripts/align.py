#!/usr/bin/env python3
"""
This script aligns unique RNA motifs with non-unique RNA motifs and generates alignment files.

The script reads two JSON files containing unique and non-unique RNA motifs, respectively. It then maps the unique motifs to the non-unique motifs, counts occurrences, and generates alignment files for both unique and non-unique motifs. Additionally, it prints out information about unique motifs that belong to more than one family.

Functions:
    align(path_uniq, path_nonuniq):
        Maps unique motifs to non-unique motifs, counts occurrences, and generates alignment files.

Usage:
    python align.py <path_to_unique_motifs_json> <path_to_non_unique_motifs_json>

Arguments:
    path_to_unique_motifs_json: Path to the JSON file containing unique RNA motifs.
    path_to_non_unique_motifs_json: Path to the JSON file containing non-unique RNA motifs.

Example:
    python align.py unique_motifs.txt non_unique_motifs.txt # each file contains one JSON object per line
"""
import sys
import json
from collections import defaultdict, Counter

from uniq import Node, get_length


# mapping between unique motifs and non-unique motifs
def align(path_uniq, path_nonuniq):
    motifs_uniq = []
    dotbracket2motif = {}
    with open(path_uniq) as f:
        for i, line in enumerate(f):
            js = json.loads(line)
            js['id_uniq'] = i+1 # 1-indexed
            js['family2count'] = {}
            js['occurrences'] = []
            js['length'] = get_length(js)
            motif_tree = Node(js['motif'])
            dotbracket2motif[str(motif_tree)] = js['id_uniq']
            for newtree in motif_tree.rotated(0):
                dotbracket2motif[str(newtree)] = i+1
            motifs_uniq.append(js)

    print('count of uniq motifs:', len(motifs_uniq))

    motifs_nonuniq = []
    with open(path_nonuniq) as f:
        counter = Counter()
        for i, line in enumerate(f):
            js = json.loads(line)
            counter[js['id']] += 1
            motif_tree = Node(js['motif'])
            motif_str = str(motif_tree)
            assert motif_str in dotbracket2motif
            id_uniq = dotbracket2motif[motif_str]
            js['id_uniq'] = id_uniq
            js['id_in_structure'] = js['id'] + '_ymotif' + str(counter[js['id']])+'_mode0'
            family = js['id'].split('_')[0]
            motifs_uniq[id_uniq-1]['family2count'][family] = motifs_uniq[id_uniq-1]['family2count'].get(family, 0) + 1
            motifs_uniq[id_uniq-1]['occurrences'].append(js['id_in_structure'])
            motifs_nonuniq.append(js)

    with open(path_uniq+'.align', 'w') as f:
        for motif in motifs_uniq:
            f.write(json.dumps(motif)+'\n')

    with open(path_nonuniq+'.align', 'w') as f:
        for motif in motifs_nonuniq:
            f.write(json.dumps(motif)+'\n')

    for motif in motifs_uniq:
        if len(motif['family2count']) > 1:
            d = {}
            d['id_uniq'] = motif['id_uniq']
            d['family2count'] = motif['family2count']
            d['occurrences'] = motif['occurrences']
            print(d)


if __name__ == "__main__":
    uniqs = defaultdict(list) # loop-signature -> [motifs]
    path_uniqs = sys.argv[1]
    path_nonuniq = sys.argv[2]
    print('path_uniqs:', path_uniqs)
    print('path_nonuniq:', path_nonuniq)
    align(path_uniqs, path_nonuniq)
