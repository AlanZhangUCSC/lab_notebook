# Evaluate Alignments

[Difference between panMAN alignment MAFFT alignment](#Difference-between-panMAN-alignment-MAFFT-alignment)

[Quantify the prevalence of problematic block states](#quantify-the-prevalence-of-problematic-block-states)

## Difference between panMAN alignment MAFFT alignment
Inside evaluate_alignments, `parse_newick_to_pairs.py` is a script that would parse the newick tree and output the parent-child pairs in depth first order. The sequences of each parent-child pair will be compared.

Following commands will generate the pairs, run the alignment comparison, and plot the alignment differences.

```
python parse_newick_to_pairs.py data/hiv/hiv_panman.newick > data/hiv/hiv20000_pairs.tsv

sbatch evaluate_alignments/evaluate_alignments_hiv.sh

python evaluate_alignments/plot_alignment_diff.py evaluate_alignments/out/hiv20000_alignment_differences.tsv evaluate_alignments/data/hiv/hiv20000_pairs.tsv hiv evaluate_alignments/out/hiv20000_alignment_differences.png
```

RSV
![RSV](out/rsv4000_alignment_differences.png)

\
New RSV tree
![RSV](out/rsv4000_alignment_differences_new.png)


\
SARS
![SARS](out/sars20000_alignment_differences.png)

\
New SARS tree
![SARS](out/sars20000_alignment_differences_new.png)

\
HIV
![HIV](out/hiv20000_alignment_differences.png)

## Quantify the prevalence of problematic block states

[Examples](#example)

[Scripts](#scripts)

### Example
Use command below to identify examples of alignments where there could be potential problematic block states
```
awk 'BEGIN {OFS="\t"} NR==1 {print $0, "Difference"} NR>1 {print $0, $3-$4}' out/rsv4000_alignment_differences_new.tsv  > out/rsv4000_alignment_differences_new_diff.tsv

tail -n+2 out/rsv4000_alignment_differences_new_diff.tsv | sort -k 5,5 -gr | less
```

Take this entry for example
```
Parent  Child   PanMAN_Diff     MAFFT_Diff      Difference
node_3806       OM857267.1      1542    51      1491
```

PanMAN alignment difference between the two nodes is 1542 while the MAFFT alignment difference is only 51. This drastic difference indicates unopmtimized block states.

Now let's check where misalignments happen betwen the two sequences on a panMAN
```
mkdir example
cp data/rsv/rsv4000_fasta_aligned/{node_3806,OM857267_1}.fasta example/
python3 compare_aligned_sequences.py example/OM857267_1.fasta example/node_3806.fasta --print_ranges
```

```
length  range
2       3643-3644
1       3646-3646
1       3650-3650
7       3652-3658
3       3663-3665
4       3670-3673
3       3675-3677
2       3679-3680
2       3682-3683
39      3685-3723
146     249667-249812
120     249841-249960
146     251858-252003
120     252032-252151
1       284262-284262
1       337600-337600
20      349094-349113
33      349116-349148
160     349205-349364
20      349654-349673
33      349676-349708
160     349765-349924
87      363086-363172
104     363259-363362
15      363364-363378
48      363380-363427
87      364210-364296
104     364383-364486
15      364488-364502
48      364504-364551
10      367685-367694
```

These ranges seems like blocks that could be misaligned. Both 146-size to 120-size blocks have distance 29.
```
146     249667-249812
120     249841-249960
146     251858-252003
120     252032-252151
```

Let's see if their sequences look identical
```
(cat <(tail -n+2 example/OM857267_1.fasta | tr -d '\n') <(echo) <(tail -n+2 example/node_3806.fasta | tr -d '\n')) > example/merged_sequence.txt

cut -c 249668-249813 example/merged_sequence.txt 
cut -c 249842-249961 example/merged_sequence.txt 
cut -c 251859-252004 example/merged_sequence.txt 
cut -c 252033-252152 example/merged_sequence.txt 
```

Comparing the two ranges of size 146:
```
--------------------------------------------------------------------------------------------------------------------------------------------------
AAACACTATACTTGATGACTTCAAAGTGAGTCTAGAATCTATAGGTAGTTTGACACAAGAATTAGAATATAGAGGTGAAAGTCTATTATGCAGTTTAATATTTAGAAATGTATGGTTATATAATCAAATTGCATTACAACTTAAAA

AAACACTATACTTGATGACTTCAAAGTGAGTCTAGAATCTATAGGTAGTTTGACACAAGAATTAGAATATAGAGGTGAAAGTCTATTATGCAGTTTAATATTTAGAAATGTATGGTTATATAATCAAATTGCATTACAACTTAAAA
--------------------------------------------------------------------------------------------------------------------------------------------------
```

Coparing the two ranges of size 120:
```
------------------------------------------------------------------------------------------------------------------------
ATCATGCATTATGTAACAACAAATTATATTTGGATATATTAAAAGTTCTAAAACACTTAAAAACCTTTTTTAATCTTGATAACATTGATACAGCATTAACATTGTATATGAATTTGCCTA
                      :
ATCATGCATTATGTAACAACAAGTTATATTTGGATATATTAAAAGTTCTAAAACACTTAAAAACCTTTTTTAATCTTGATAACATTGATACAGCATTAACATTGTATATGAATTTGCCTA
------------------------------------------------------------------------------------------------------------------------
```

So they are indeed highly similar. The 146-long blocks are identical while the 120-long blocks are only 1 snp apart.

Now check the sequences between the ranges
```
cut -c 249815-249842 example/merged_sequence.txt
cut -c 252006-252033 example/merged_sequence.txt 
cut -c 249963-251859 example/merged_sequence.txt 
```

```
----------------------------
---------------------------A

---------------------------A
----------------------------

------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------A
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
```

The <code style="color : red">test</code> sequences between each 120-block and 146-block are identical. The sequence between two 140-120-block pairs is all gaps for one of the nodes. This allows use to swap the 120-146 block pairs for one of the nodes to fix the alignments.

### Scripts