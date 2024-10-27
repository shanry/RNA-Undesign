#!/usr/bin/env python3
import os
import argparse
from collections import Counter

import json
import glob

'''
python parser.py --path log.txt

where each line in log.txt is a json string representaing a motif
{"bpairs":[[72,79]],"end":77,"id":"5s_ref2","ipairs":[[73,77]],"ismin":true,"motif":{"id":"5s_ref2","root":{"child_id":1,"children":[{"child_id":1,"loops":[3],"neighbors":[0],"type":"H"}],"loops":[0,1],"neighbors":[0,1],"type":"B"}},"plotstr":"\"5s_ref2,(.....(.((((((.....((((((.............))))..))....))))))..).((.((....(((((...).))))....)).))...).......,73 80 GREEN Fomark 73 80 0.70 0.5 colorpair 74 78 0.1667 1.0 colorpair \"","start":73,"time":0.004015999846160412,"y_rivals":["(......)"],"y_star":"(.....(.((((((.....((((((.............))))..))....))))))..).((.((....(((((...).))))....)).))...).......","y_sub":"((...).)"}

output the list of string to plot corresponding motif.

parse the log of applying undesignability detection to 5s:

python parser.py --path 5s_ref.txt.pn.log.20240221092905.txt

the output file is named as 5s_ref.txt.pn.plotstr.20240221092905.txt

'''


def extract_pairs_list(ss):
    pairs = []
    stack = []
    for i, c in enumerate(ss):
        if c != '(' and c != ')':
            pass
        elif c=="(":
            stack.append(i)
        elif c==")":
            j = stack.pop()
            pairs.append((j, i))
        else:
            raise ValueError(f"wrong structure at position {i}: {c}")
    return pairs


def shrink(i, j, bpairs, SIZE = 2):
    inew, jnew = i, j
    for pair in bpairs:
        delta = pair[1] - pair[0] - 1 - SIZE
        if pair[1] < i:
            inew -= delta
            jnew -= delta
        elif pair[0] > i and pair[1] < j:
            jnew -= delta
        elif pair[0] == i or pair[1] == j:
            jnew -= delta
    return inew, jnew


def replace_substring(string, start, end, replacement):
    return string[:start] + replacement + string[end:]


# intial version without extending the outmost pair
def get_mplotstr2(m):
    # shrink y
    # ynew = m['y_sub']
    m['bpairs'] = sorted(m['bpairs'])
    m['ipairs'] = sorted(m['ipairs'])
    if m['bpairs'][0][0] >= 0:
        ysub = m['y_star'][m['bpairs'][0][0]: m['bpairs'][0][1]+1]
        ynew = m['y_star'][m['bpairs'][0][0]: m['bpairs'][0][1]+1]
        start = m['bpairs'][0][0]
    else:
        ysub = m['y_star']
        ynew = m['y_star']
        # m['bpairs'][0] = [0, len(ysub)-1]
        start = 0
    # print(m['bpairs'])
    # print(m['ipairs'])
    # print(ynew)
    SIZE = 3
    delta = 0
    # start = m['bpairs'][0][0]
    for pair in m['bpairs'][1:]:
        # print(pair, pair[1] - pair[0] - 1 - 2)
        # print(len(ynew))
        ynew = replace_substring(ynew, pair[0]+1-start-delta, pair[1]-start-delta, "...")
        # print(len(ynew))
        delta += pair[1] - pair[0] - 1 - SIZE
        # print(f"delta: {delta}")
        # print(f"ydiff: {len(m['y_sub']) - len(ynew)}")
        assert len(ysub) - len(ynew) == delta, f"{len(ysub), len(ynew)}"
        # print()
    # shrink pairs
    bpairs = []
    ipairs = []
    for i, j in m['bpairs']:
        # print(i, j)
        if i >= 0:
            bpairs.append(shrink(i, j, m['bpairs'][1:], SIZE))
        else:
            bpairs.append((i, j))
        # print(bpairs)
        # print()
    # print(m['bpairs'])
    # print(bpairs)
    for i, j in m['ipairs']:
        ipairs.append(shrink(i, j, m['bpairs'][1:], SIZE))
    # print(m['ipairs'])
    # print(ipairs)
    # plot string 
    i0, j0 = bpairs[0]
    if i0 >= 0:
        greenplot = f"{i0+1-i0} {j0+1-i0} GREEN Fomark "
    else:
        assert len(m['y_star']) == len(m['y_sub'])
        greenplot = f"{1} {len(ynew)} GREEN Fomark "
        i0 = 0
    whiteplot = ""
    for i, j in bpairs[1:]:
        whiteplot += f"{i+1-i0} {j+1-i0} WHITE Fomark "
    pairplot = ""
    if bpairs[0][0] >= 0:
        pairplot += f"{i0+1-i0} {j0+1-i0} 0.667 0.5 colorpair "
    for i, j in bpairs[1:]:
        pairplot += f"{i+1-i0} {j+1-i0} 0.667 0.5 colorpair "
    for i, j in ipairs:
        pairplot += f"{i+1-i0} {j+1-i0} 0.1667 1.0 colorpair "
    poststr = ""
    for i, j in bpairs[1:]:
        poststr += f"{i+1-i0+1} {j+1-i0-1} 10 WHITE omark "
    plotstr = ",".join([m['id']+"_motif", ynew, greenplot + whiteplot + pairplot, poststr])
    return plotstr


# extend the outmost pair with two unpaired bases
def get_mplotstr(m):
    # shrink y
    # ynew = m['y_sub']
    m['bpairs'] = sorted(m['bpairs'])
    m['ipairs'] = sorted(m['ipairs'])
    if m['bpairs'][0][0] >= 0:
        ysub = m['y_star'][m['bpairs'][0][0]: m['bpairs'][0][1]+1]
        ynew = m['y_star'][m['bpairs'][0][0]: m['bpairs'][0][1]+1]
        start = m['bpairs'][0][0]
    else:
        ysub = m['y_star']
        ynew = m['y_star']
        # m['bpairs'][0] = [0, len(ysub)-1]
        start = 0
    # print(m['bpairs'])
    # print(m['ipairs'])
    # print(ynew)
    SIZE = 3
    delta = 0
    # start = m['bpairs'][0][0]
    for pair in m['bpairs'][1:]:
        # print(pair, pair[1] - pair[0] - 1 - 2)
        # print(len(ynew))
        ynew = replace_substring(ynew, pair[0]+1-start-delta, pair[1]-start-delta, "...")
        # print(len(ynew))
        delta += pair[1] - pair[0] - 1 - SIZE
        # print(f"delta: {delta}")
        # print(f"ydiff: {len(m['y_sub']) - len(ynew)}")
        assert len(ysub) - len(ynew) == delta, f"{len(ysub), len(ynew)}"
        # print()
    # shrink pairs
    bpairs = []
    ipairs = []
    for i, j in m['bpairs']:
        if i >= 0:
            bpairs.append(shrink(i, j, m['bpairs'][1:], SIZE))
        else:
            bpairs.append((i, j))
    for i, j in m['ipairs']:
        ipairs.append(shrink(i, j, m['bpairs'][1:], SIZE))
    i0, j0 = bpairs[0]
    i_base = i0
    flag = False
    if i0 >= 0:
        ynew = "." + ynew + "." # extend the outmost pair with two unpaired bases
        i_base -= 1
        greenplot = f"{i0+1-i_base} {j0+1-i_base} GREEN Fomark "
        flag = True
    else:
        assert len(m['y_star']) == len(m['y_sub'])
        greenplot = f"{1} {len(ynew)} GREEN Fomark "
        i_base = 0
    whiteplot = ""
    for i, j in bpairs[1:]:
        whiteplot += f"{i+1-i_base} {j+1-i_base} WHITE Fomark "
    pairplot = ""
    if bpairs[0][0] >= 0:
        pairplot += f"{i0+1-i_base} {j0+1-i_base} 0.667 0.5 colorpair "
    for i, j in bpairs[1:]:
        pairplot += f"{i+1-i_base} {j+1-i_base} 0.667 0.5 colorpair "
    for i, j in ipairs:
        pairplot += f"{i+1-i_base} {j+1-i_base} 0.1667 1.0 colorpair "
    poststr = ""
    for i, j in bpairs[1:]:
        poststr += f"{i+1-i_base+1} {j+1-i_base-1} 10 WHITE omark "
    if flag:
        poststr += f"1 1 10 WHITE omark "
        poststr += f"{len(ynew)} {len(ynew)} 10 WHITE omark "
    plotstr = ",".join([m['motif_id'], ynew, greenplot + whiteplot + pairplot, poststr])
    return plotstr

# plot motifs in original structure
def get_yplotstr(m):
    i0, j0 = m['bpairs'][0]
    if i0 >= 0:
        greenplot = f"{i0+1} {j0+1} GREEN Fomark "
    else:
        assert len(m['y_star']) == len(m['y_sub']), f"{m['id']}: {json.dumps(m)} \n {len(m['y_star'])} != {len(m['y_sub'])}"
        greenplot = f"{1} {len(m['y_star'])} GREEN Fomark "
        i0 = 0
        j0 = len(m['y_star'])-1
    whiteplot = ""
    for i, j in m['bpairs'][1:]:
        whiteplot += f"{i+1} {j+1} WHITE Fomark "
    pairplot = ""
    if m['bpairs'][0][0] >= 0:
        pairplot += f"{i0+1} {j0+1} 0.667 0.5 colorpair "
    for i, j in m['bpairs'][1:]:
        pairplot += f"{i+1} {j+1} 0.667 0.5 colorpair "
    for i, j in m['ipairs']:
        pairplot += f"{i+1} {j+1} 0.1667 1.0 colorpair "
    plotstr = ",".join([m['ymotif_id'], m['y_star'], greenplot + whiteplot + pairplot])
    return plotstr


def get_motifs(file):
    motifs_list = []
    for line in open(file):
        motifs_list.append(json.loads(line))
    return motifs_list


def get_rival_motif_plotstr(m, ynew):
    # shrink y
    # ynew = m['y_sub']
    m['bpairs'] = sorted(m['bpairs'])
    m['ipairs'] = sorted(m['ipairs'])
    if m['bpairs'][0][0] >= 0:
        ysub = m['y_star'][m['bpairs'][0][0]: m['bpairs'][0][1]+1]
        # ynew = m['y_star'][m['bpairs'][0][0]: m['bpairs'][0][1]+1]
        start = m['bpairs'][0][0]
    else:
        ysub = m['y_star']
        # ynew = m['y_star']
        # m['bpairs'][0] = [0, len(ysub)-1]
        start = 0
    # print(m['bpairs'])
    # print(m['ipairs'])
    # print(ynew)
    SIZE = 3
    delta = 0
    # start = m['bpairs'][0][0]
    for pair in m['bpairs'][1:]:
        # print(pair, pair[1] - pair[0] - 1 - 2)
        # print(len(ynew))
        ynew = replace_substring(ynew, pair[0]+1-start-delta, pair[1]-start-delta, "...")
        # print(len(ynew))
        delta += pair[1] - pair[0] - 1 - SIZE
        # print(f"delta: {delta}")
        # print(f"ydiff: {len(m['y_sub']) - len(ynew)}")
        assert len(ysub) - len(ynew) == delta, f"{len(ysub), len(ynew)}"
        # print()
    # shrink pairs
    bpairs = []
    ipairs = []
    for i, j in m['bpairs']:
        if i >= 0:
            bpairs.append(shrink(i, j, m['bpairs'][1:], SIZE))
        else:
            bpairs.append((i, j))
    for i, j in m['ipairs']:
        ipairs.append(shrink(i, j, m['bpairs'][1:], SIZE))
    i0, j0 = bpairs[0]
    i_base = i0
    flag = False
    if i0 >= 0:
        ynew = "." + ynew + "." # extend the outmost pair with two unpaired bases
        i_base -= 1
        greenplot = f"{i0+1-i_base} {j0+1-i_base} GREEN Fomark "
        flag = True
    else:
        assert len(m['y_star']) == len(m['y_sub'])
        greenplot = f"{1} {len(ynew)} GREEN Fomark "
        i_base = 0
    whiteplot = ""
    for i, j in bpairs[1:]:
        whiteplot += f"{i+1-i_base} {j+1-i_base} WHITE Fomark "
    pairplot = ""
    if bpairs[0][0] >= 0:
        pairplot += f"{i0+1-i_base} {j0+1-i_base} 0.667 0.5 colorpair "
    for i, j in bpairs[1:]:
        pairplot += f"{i+1-i_base} {j+1-i_base} 0.667 0.5 colorpair "

    ymask = list(ynew)
    for i, j in bpairs[1:]:
        for k in range(i, j+1):
            ymask[k-i_base] = "*"

    ymask = "".join(ymask)

    pairs_rival = extract_pairs_list(ymask)

    ipairs_rival = []

    for i, j in pairs_rival:
        if (i, j) != (i0-i_base, j0-i_base):
            ipairs_rival.append((i, j))

    # print(ymask)
    # print(ipairs_rival)
    # print(i0, j0)
    # print(i_base)
    # print()

    for i, j in ipairs_rival:
        pairplot += f"{i+1} {j+1} 0.1667 1.0 colorpair "
    poststr = ""
    for i, j in bpairs[1:]:
        poststr += f"{i+1-i_base+1} {j+1-i_base-1} 10 WHITE omark "
    if flag:
        poststr += f"1 1 10 WHITE omark "
        poststr += f"{len(ynew)} {len(ynew)} 10 WHITE omark "
    plotstr = ",".join([m['rival_id'], ynew, greenplot + whiteplot + pairplot, poststr])
    return plotstr


def main_plot_rival(motifs, mode='r'):
    counter = Counter()
    plot_lines = []
    for im, motif in enumerate(motifs):
        counter[motif['id_uniq']] += 1
        motif['motif_id'] = '_'.join([str(motif['id_uniq']), "motif"+str(counter[motif['id_uniq']])])
        for ir, rival in enumerate(motif['y_rivals']):
            motif['rival_id'] = '_'.join([motif['motif_id'], "rival"+str(ir+1)])
            plot_lines.append('"'+get_rival_motif_plotstr(motif, rival)+'"'+"\n")
    path_output = args.path + f'.{mode}plotstr'
    with open(path_output, 'w') as f:
        f.writelines(plot_lines)


def main_plot(motifs, mode='m'):
    counter = Counter()
    plot_lines = []
    motif2plotstr = get_mplotstr if mode == 'm' else get_yplotstr
    # print('motif2plotstr:', motif2plotstr)
    for i, motif in enumerate(motifs):
        counter[motif['id']] += 1
        if mode == 'm':
            motif['motif_id'] = '_'.join([str(motif['id']), "motif"+str(counter[motif['id']])])
            if 'id_uniq' in motif:
                motif['motif_id'] = str(motif['id_uniq']) # '_'.join(('uniq', str(motif['id_uniq'])))
        elif mode == 'y':
            motif['ymotif_id'] = '_'.join([str(motif['id']), "ymotif"+str(counter[motif['id']])])
            if 'id_in_structure' in motif:
                motif['ymotif_id'] = motif['id_in_structure']
        plot_lines.append('"'+motif2plotstr(motif)+'"'+"\n")
    path_output = args.path + f'.{mode}plotstr'
    with open(path_output, 'w') as f:
        f.writelines(plot_lines)
    print(path_output)


def main_time(motifs):
    time_list = []
    for motif in motifs:
        time_list.append(motif['time'])
    print('time_list:', time_list)
    print('max time:', max(time_list))
    print('min time:', min(time_list))
    print('average time:', sum(time_list)/len(time_list))


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("--path", '-p', type=str, default='9families_free/telomerase_ref.txt.pn.log.20240923174211.txt')
    parser.add_argument("--mode", '-m', type=str, default='m') # m: motif, y: original structure, r: rival, t: running time

    args = parser.parse_args()
    # print('args:')
    # print(args)

    json_motifs = get_motifs(args.path)
    if args.mode == 'r':
        main_plot_rival(json_motifs)
    elif args.mode == 't':
        main_time(json_motifs)
    else:
        main_plot(json_motifs, args.mode)
