# RNA-Undesign
Algorithms for Identifying the Undesignability of RNA secondary structures

## Introduction

This project explores the undesignability of RNA secondary structures, under the context of RNA design, building upon the research presented in the paper by T. Zhou et al. [1][2].

## References

[1] Zhou, T., Tang, W.Y., Mathews, D.H. and Huang, L. 
"Undesignable RNA Structure Identification via Rival Structure Generation and Structure Decomposition." 
(RECOMB 2024, [arXiv preprint](https://arxiv.org/abs/2311.08339)) \
[2] Zhou, T., Tang, W.Y., Mathews, D.H. and Huang, L. 
"Scalable Identification of Minimum Undesignable RNA Motifs on Loop-Pair Graphs" ([arXiv preprint](https://arxiv.org/abs/2402.17206))

## Dependencies
- GCC 4.8.5 or above

## Build
### Linux
```
# Test with CentOS
$ make main   
$ make main_nosh  # turn off special hairpins
```
### Mac
```
# Test with Apple Silicon
$ make main_mac   
$ make main_nosh_mac  # turn off special hairpins
```
## Rotation Invariance of Motifs
<p align="center">
<img src="figs/rotation_invariance.png" width="700">
</p>

----

<p align="center">
<img src="figs/rotation_invariance_86.png" width="640">
</p>

## Environment Variables
```
$export OMP_NUM_THREADS=8 # parallel computing eabled by OpenMP
$export PATH_UNDESIGNABLE_LIB=path/to/lib_undesignable.txt
$export PATH_DESIGNABLE_LIB/path/to/motifs/libs/lib_designable.txt
```
## FastMotif 
```
$echo ".((......((......))......((......((......))......((......))......))......))....." | ./bin/main --alg fastmotif
```
## RIGEND

### Undesignability Alg1
```
$export OMP_NUM_THREADS=8 # parallel computing eabled by OpenMP
$./bin/main --alg 1
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
