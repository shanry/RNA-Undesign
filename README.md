# RNA-Undesign
Undesignability of RNA secondary structures

## Build
```./make```

## Energy Evaluation
``` 
$./bin/eval eval
CACACGCACUACAAAAUGUCCAAAGGAAAAGGCACCACCAGCAAAGCACCAAAGGUAAGGGGAAAAG
.....((.((.((...)).((...))...)))).((.((.((...)).((...))...)))).....
(output)total energy: -3.00
```
## Critical positions
```
$ ./bin/eval critical
..((((((((.......(.((((((....)))))).(((((((....))))))).).......))))))))..
..((((((((.........((((((....)))))).(((((((....))))))).........))))))))..
multi_outside (9, 63) : 9 10 62 63 
multi_outside (17, 55) : 17 18 54 55 
internal (9, 63), (17, 55) : 9 10 16 17 55 56 62 63 
critical positions: 9, 10, 16, 17, 18, 54, 55, 56, 62, 63, 
```
## Energy difference
```
$ ./bin/eval test
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
## Undesignability Alg1
```
python main.py 1
CUAAGGACCACCCGGGAAACCAUAAGGGGCGAGAAAUCGAGGAUCAACAGCGCAGGGAAAACGAACCAUCCGAAAGGAAGCAAGCAAAAAAAGAAAAAAAAAAAA
.....((((.(((((....))....))).(((....))).)).))....((((.(((....)...)).(((....))).))..))....................
.....((((.(((((....))....))).(((....))).)).))....((((.((.........)).(((....))).))..))....................
diffs: [(55, 65), (56, 61), 57, 60, 64]
total number of enumerations: 2304
  1000 CUAAGGACCACCCGGGAAACCAUAAGGGGCGAGAAAUCGAGGAUCAACAGCGCAGAUGAACAGAUUCAUCCGAAAGGAAGCAAGCAAAAAAAGAAAAAAAAAAAA 4.4 seconds
                                                              ^^^  ^^  ^^                                       
  2000 CUAAGGACCACCCGGGAAACCAUAAGGGGCGAGAAAUCGAGGAUCAACAGCGCAGUCAAAUGGAUGCAUCCGAAAGGAAGCAAGCAAAAAAAGAAAAAAAAAAAA 8.8 seconds
                                                              ^^^  ^^  ^^                                       
the puzzle .....((((.(((((....))....))).(((....))).)).))....((((.(((....)...)).(((....))).))..)).................... is unsolvable.
```