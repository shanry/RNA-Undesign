import os
import argparse 

import json
import glob


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
    plotstr = ",".join([m['id']+"_m", ynew, greenplot + whiteplot + pairplot, poststr])
    return plotstr


def get_motifs(file):
    motifs_list = []
    for line in open(file):
        motifs_list.append(json.loads(line))
    return motifs_list


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("--path", '-p', type=str, default='9families_free/telomerase_ref.txt.pn.log.20240923174211.txt')

    args = parser.parse_args()
    print('args:')
    print(args)

    json_motifs = get_motifs(args.path)
    plot_lines = []
    for im, motif in enumerate(json_motifs):
        # print(im)
        plot_lines.append('"'+get_mplotstr2(motif)+'"'+"\n")
    path_output = os.path.basename(args.path).replace('log', 'plotstr')
    with open(path_output, 'w') as f:
        f.writelines(plot_lines)