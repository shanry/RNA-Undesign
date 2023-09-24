
nuc_pair_all = ['AU', 'UA', 'CG', 'GC', 'UG', 'GU']

def extract_pairs(ss):
    pairs = list(range(len(ss)))
    stack = []
    for i, c in enumerate(ss):
        if c=='.':
            pass
        elif c=="(":
            stack.append(i)
        elif c==")":
            j = stack.pop()
            pairs[j] = i
            pairs[i] = j
        else:
            raise ValueError(f"wrong structure at position {i}: {c}")
    return pairs

def extract_pairs_list(ss):
    pairs = []
    stack = []
    for i, c in enumerate(ss):
        if c=='.':
            pass
        elif c=="(":
            stack.append(i)
        elif c==")":
            j = stack.pop()
            pairs.append((j, i))
        else:
            raise ValueError(f"wrong structure at position {i}: {c}")
    return pairs


def pairs_match(ss): # find the pairs in a secondary structure, return a dictionary
    assert len(ss) > 5
    pairs = dict()
    stack = []
    for i, s in enumerate(ss):
        if s==".":
            pass
        elif s=="(":
            stack.append(i)
        elif s==")":
            j = stack.pop()
            assert j < i
            pairs[j] = i
            pairs[i] = j
        else:
            raise ValueError(f'the value of structure at position: {i} is not right: {s}!')
    return pairs


def struct_dist(s1, s2):
    assert len(s1) == len(s2), f"len(s1)={len(s1)}, len(s2)={len(s2)}"
    pairs_1 = pairs_match(s1)
    pairs_2 = pairs_match(s2)
    union = len(pairs_1.keys()|pairs_2.keys())
    overlap = len(s1) - union
    for k in pairs_1:
        if k in pairs_2 and pairs_1[k]==pairs_2[k]:
            overlap += 1
    return len(s1) - overlap


def check_compatible(seq, ss, report=False):
    pairs_list = extract_pairs_list(ss)
    for (i, j) in pairs_list:
        if seq[i]+seq[j] not in nuc_pair_all:
            if report:
                print(i, j, ":", seq[i]+seq[j])
            return False
    return True


if __name__ == '__main__':
    ss = "..........((((....))))((((....))))((((...))))"
    pairs = extract_pairs(ss)
    print('structure:', ss)
    print('pairs:', pairs)
    pairs_list = extract_pairs_list(ss)
    print("pairs_list: ", pairs_list)
    
    seq = "AAAAGCCGGGGGACGCAGGGACACACCGGCAAAA"
    ss  = "....((((((((.(....)).).).)))))...."
    print(f"check compatible: {check_compatible(seq, ss)}")