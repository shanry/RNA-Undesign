# RNA-Undesign
Undesignability of RNA secondary structures

## Introduction

This project explores the unndesignability of RNA secondary structures, under the context of RNA design, building upon the research presented in the paper by T. Zhou et al. [1][2].

## References

[1] Zhou, T., Tang, W.Y., Mathews, D.H. and Huang, L. 
"Undesignable RNA Structure Identification via Rival Structure Generation and Structure Decomposition." 
arXiv preprint arXiv:2311.08339 (2023). To appear in Proceedings of RECOMB 2024. \
[2] Zhou, T., Tang, W.Y., Mathews, D.H. and Huang, L. 
"Scalable Identification of Minimum Undesignable RNA Motifs on Loop-Pair Graphs"

## Dependencies
- GCC 4.8.5 or above

## Build
### Linux
```
# Test with CentOS
$./make main   
$./make main_hp3  # only turn on special hairpin of triloop
```
### Mac
```
# Test with Apple Silicon
$./make main_mac   
$./make main_hp3_mac  # only turn on special hairpin of triloop
```
## Scalable Identification of Minimum Undesignable Motifs

### Loop + 1 Neighbor
```
$export OMP_NUM_THREADS=8 # parallel computing eabled by OpenMP
$./bin/main --alg n1
..((((((((.......(.((((((....)))))).(((((((....))))))).).......))))))))..
```
### Loop + 2 Neighbors
```
$export OMP_NUM_THREADS=8 # parallel computing eabled by OpenMP
$./bin/main --alg n2
((.(..(.(....).(....).)..).(....).))
```
### Loop with 3+ Neighbors
```
$export OMP_NUM_THREADS=8 # parallel computing eabled by OpenMP
$./bin/main --alg n3
((.(((..(((.........)))..)))(...)(...((......)).).))
```
## RIGEND

### Undesignability Alg1
```
$export OMP_NUM_THREADS=8 # parallel computing eabled by OpenMP
$./bin/main --alg 1
AUAAGCGGUAAAAAAAGUGCGAAAAGCAUGAAAAAAAACAGAAAAAAAAAAAAAAAAAAAA
......(.........((((.....)))).........)......................
................((((.....))))................................
```

### Undesignability Alg2
```
$export OMP_NUM_THREADS=8 # parallel computing eabled by OpenMP
$./bin/main --alg 2
AAAAUGAGCCCCACGAAAGGAGAGUGCUCACAAA
....((((((((.(....)).).).)))))....
....(((((((..(....)..).).)))))....
```

### Undesignability Alg2 (constrained)
```
$export OMP_NUM_THREADS=8 # parallel computing eabled by OpenMP, user-defined
$./bin/main --alg 2c
UUAAGGGAAAAUCUUAGCCGAGAAAUCGGAUCCAAAGCGGCAUAAAAAAGAAAGCGCCGAAAUUCGCAGAAAUGCGAGAAAGGCAAGCAAAGAAUUCGGCAGAAAAAAUGCCGACCGGGCAAUGAAAAUUCGCCCGUGGAGCCAAGCGGG
((((((.....)))))(((((....)))).)((...(((((............(((((....((((((....))))))...)))..)).......((((((.......))))))(((((((((....))).))))).)..)))..)))))
((((((.....))))).((((....))))..((...(((((............(((((....((((((....))))))...)))..)).......((((((.......))))))(((((((((....))).))))).)..)))..)))))
```
### Undesignability Alg3
```
$export OMP_NUM_THREADS=8 # parallel computing eabled by OpenMP, user-defined
$./bin/main --alg 3
ACUAAAUGGUGAGCAGACCCAGUGGAAACACACGCAGCCGAAAGGUACCCAUCCGAGAGGAAGUCAGGCGAAAGCUAACGGAAAGAACGUAGACAGGGAGCGAGGGACAAAGACUGCAAGGGAAAGUACACAAGACAAAGUAAAAAAAGGUGAGGCAGGGGAAACCCCGGGAAACCGGUCGAAAGACGCCAGCAAACCGCAGAAACAGCCACCCAGCGAGACAGACAAAAGCGGAUACGUAGUCGACGGAAACGUAGUCAGGGGAAACCCACGCAAUCGAAAGAUAGGGAGUCGGUGAAAACCAGAGAAAUCUACUCAAAAGAGGACAGGCAGCGGAACCCCUACACCGAAAAAA
.......((((.((((.(((.(((....))).(((.(((....))).(((.(((....))).(((.(((....))).(((.......))).))).))).))).))).(...).)))).((((...((.(((...((...))........))).(((.(((....)))(((....)))(((....)))))).))...((((.(...).(((.(((.(((.(((.(((....(((....))).))).(((....))).))).(((....))).))).(((....))).))).((((((....)))(((....))).(((....)))))).))).))))...)))).)))).......
```
## Utilities
### Energy Evaluation
``` 
$./bin/main --alg eval
CACACGCACUACAAAAUGUCCAAAGGAAAAGGCACCACCAGCAAAGCACCAAAGGUAAGGGGAAAAG
.....((.((.((...)).((...))...)))).((.((.((...)).((...))...)))).....
(output)total energy: -3.00
```
### Differential positions
```
$./bin/main --alg dp
..((((((((.......(.((((((....)))))).(((((((....))))))).).......))))))))..
..((((((((.........((((((....)))))).(((((((....))))))).........))))))))..
```
### Energy difference
```
$./bin/main --alg ed
GGGAGACCCAAAAAAAAGGGCAACUGCAAAAAGGAGACAGCACCCCGAAAAAAGACUGGAAAAAGGGCGAAAAGCUCGAAAAACACGACCAACGGAAAACAGGACGAAAGAGAACAAGCAAGCCAAAGGGAAACAGACUAAAAACGCGAAAGCGACUGCAAAGGGGGAGAAAAAGCGACCCUGAACGAAAAAGGGGCGAAAAAUUGGAACAAAAAAAGGAGGGGGGAAAGGAAAGUCAAAGACACUCGAAACGAGUGAGCGGGCAAAAAAAAAAACGGGGGAUGAAUAACGGACGGAAACGCGGCGGAAAGCGAAAAAAAGAAAAACGUCGUACGGACUACUGGGGUGCAAAAAAAAGGAGGGGCGCAAAAAGGAAAAAACAGGGUCCACUA
((..(((((........(.((..((.(.....(....).((((((((.........((......((((.....))))......))...((..((.......(..(......)..)..((..((.....(....).((((.....(((.....(..((.(...).))..).....))).((((...(.......(..(.((...)).)..).......)...))))........)))).....((((((...)))))).))..))...........))..)).........(..((....)((((((.....(........).....))))))..)..)...))))))))........).))..)).).....(.......).)))))..)).
((..(((((........(.((..((.(.....(....).((((((((.........((......((((.....))))......))...((..((.......(..(......)..)..((..((............((((.....(((.....(..((.(...).))..).....))).((((...(.......(..(.((...)).)..).......)...))))........)))).....((((((...)))))).))..))...........))..)).........(..((.....((((((.....(........).....)))))).))..)...))))))))........).))..)).).....(.......).)))))..)).
ref1 energy: -68.20
ref2 energy: -72.70
GGGAGACCCAAAAAAAAGGGCAACUGCAAAAAGGAGACAGCACCCCGAAAAAAGACUGGAAAAAGGGCGAAAAGCUCGAAAAACACGACCAACGGAAAACAGGACGAAAGAGAACAAGCAAGCCAAAGGGAAACAGACUAAAAACGCGAAAGCGACUGCAAAGGGGGAGAAAAAGCGACCCUGAACGAAAAAGGGGCGAAAAAUUGGAACAAAAAAAGGAGGGGGGAAAGGAAAGUCAAAGACACUCGAAACGAGUGAGCGGGCAAAAAAAAAAACGGGGGAUGAAUAACGGACGGAAACGCGGCGGAAAGCGAAAAAAAGAAAAACGUCGUACGGACUACUGGGGUGCAAAAAAAAGGAGGGGCGCAAAAAGGAAAAAACAGGGUCCACUA
ref1 energy = 7.20, ref2 energy = 2.70
pass test: true
e1 - e2: 4.50
delta  : 4.50
```
or
```
$./bin/main test_diff  < data/seq_refs.txt # batched input
```