import sys
import copy
import time
import subprocess

from utils.structure import check_compatible, pairs_match

def rna_eval(seq, ref):
    cmds = f"./bin/eval eval"
    rt = subprocess.run(cmds.split(), stdout=subprocess.PIPE, input="\n".join([seq, ref]).encode())
    lines = rt.stdout.decode('utf-8').strip().split('\n')
    return lines

def get_total_energy(seq, ref):
    outputs = rna_eval(seq, ref)
    assert "total energy" in outputs[-1]
    return float(outputs[-1].split(":")[-1])

def enum(diff_list, seq_list):
    if not diff_list:
        yield seq_list
    elif type(diff_list[-1]) is int:
        for seq_list_last in enum(diff_list[:-1], seq_list):
            for nuc in "ACGU":
                seq_list_new = copy.deepcopy(seq_list_last)
                seq_list_new[diff_list[-1]] = nuc
                yield seq_list_new
    else:
        assert type(diff_list[-1]) is tuple
        i, j = diff_list[-1]
        for seq_list_last in enum(diff_list[:-1], seq_list):
            for nuc_pair in ["GC", "CG", "AU", "UA", "GU", "UG"]:
                seq_list_new = copy.deepcopy(seq_list_last)
                seq_list_new[i] = nuc_pair[0]
                seq_list_new[j] = nuc_pair[1]
                yield seq_list_new
                 
def alg_1(target, struct, diffs, seqstr):
    start = time.time()
    X = []
    label = [" "]*len(target)
    for diff in diffs:
        if type(diff) is int:
            label[diff] = '^'
        else:
            label[diff[0]] = '^'
            label[diff[1]] = '^'
    for i, seq_list in enumerate(enum(diffs, list(seqstr))):
        seq = "".join(seq_list)
        if (i+1)%1000 == 0:
            print(f"{i+1:6d}", seq, f"{time.time()-start:.1f} seconds")
            print(f"      ", "".join(label))
        e_target = get_total_energy(seq, target)
        if check_compatible(seq, struct):
            e_struct = get_total_energy(seq, struct)
        else:
            e_struct = float("inf")
        if e_struct >= e_target:
            # print(seq)
            # print(f"{e_struct:.2f} >= {e_target:.2f}")
            X.append(seq)
    if X:
        print(f"len(X): {len(X)}")
        for i in range(min(len(X), 10)):
            print(X[i])
    else:
        print(f'the puzzle {target} is unsolvable.')
    return X

def diff_positions(ref1, ref2):
    cmds = f"./bin/eval critical"
    rt = subprocess.run(cmds.split(), stdout=subprocess.PIPE, input="\n".join([ref1, ref2]).encode())
    lines = rt.stdout.decode('utf-8').strip().split('\n')
    assert "critical positions" in lines[-1]
    return eval(lines[-1].split(":")[-1])

def get_diffs(ref1, ref2):
    cr_positions = diff_positions(ref1, ref2) #[(6, 38), (16, 28), 5, 7, 37, 39, 15, 29] # ( 7, 39) GC; ( 17, 29)
    pairs_dict = pairs_match(ref1)
    diffs = []
    for p in cr_positions:
        if p in pairs_dict:
            i, j = min(p, pairs_dict[p]), max(p, pairs_dict[p])
            if (i, j) not in diffs:
                diffs.append((i, j))
        else:
            diffs.append(p)
    return diffs

def count_enum(diffs):
    num_enum = 1
    for item in diffs:
        if type(item) is tuple:
            num_enum *= 6
        else:
            num_enum *= 4
    return num_enum

def test(seq, ref1, ref2, alg=None):
    # seq = "AAAAAAGGAAAAAAAAGCCCGAAAAGGGUGAAAAAAGACAGAAAAAAAAAAAAAAAAAAAA"
    # ref1 = "......(.........((((.....)))).........)......................"
    # ref2 = "................((((.....))))................................"
    
    # seqstr = "AAAAGCCGGGGGACGCAGGCACACACCGGCAAAA"
    # target = "....((((((((.(....)).).).)))))...."
    # struct = "....(((((((..(....)..).).)))))...."
    # diffs = [(13, 18), (11, 19), (10, 21), 12, 20]
    
    # seqstr = "AACGCUCUACGGCAACGGGAGGAAAAAAAACAAAAAAAAACACGCAAAGGGAAACCACCGGGAAACCUAAAACACGCGAAAAACCCCGUAGGGAAAA"
    # target = "....((((((((....((..(...................)..((...((....))...((....))........))......))))))))))...."
    # struct = "....((((((((....((..(...................).(((...((....))...((....))........))).....))))))))))...."
    # diffs = [(42, 77), (43, 76), 41, 78]
    
    diffs = get_diffs(ref1, ref2)
    print('diffs:', diffs)
    num_enum = count_enum(diffs)
    print(f'total number of enumerations: {num_enum}')
    if alg=="1" and num_enum < 10**8:
        alg_1(ref1, ref2, diffs, seq)

if __name__ == "__main__":

    alg = None if len(sys.argv) < 2 else sys.argv[1]
    print(f"alg: {alg}")
    while True:
        seq = sys.stdin.readline().strip()
        ref1 = sys.stdin.readline().strip()
        ref2 = sys.stdin.readline().strip()
        if not seq or not ref1 or not ref2:
            break
        test(seq, ref1, ref2, alg)