# PanMANIA-notebook

This notebook will track the progress of the panMANIA project.

## 1/30/2025

To build the panMAN using Dockerfile, see `build_panman` directory.

PanMAN built using Dockerfile is broken. Currently using Summit's docker image.

## 1/31/2025 - 2/1/2025

Today I want to start by examining the alignments of sequences on a panMAN. I will write scripts that would compare the global alignments of sequences from the panMAN to the alignments of sequences using mafft. I expect to see that the mafft alignments should be better than from the panMAN because panMAN alignments are from MSA of all the sequences while mafft only compares two sequences at a time. Big differences, however, indicate particular problematic alignments, likely on the block level.

Inside evaluate_alignments, `parse_newick_to_pairs.py` is a script that would parse the newick tree and output the parent-child pairs in depth first order. The sequences of each parent-child pair will be compared.

Following commands will generate the pairs, run the alignment comparison, and plot the alignment differences.

```
python parse_newick_to_pairs.py data/hiv/hiv_panman.newick > data/hiv/hiv20000_pairs.tsv

sbatch evaluate_alignments/evaluate_alignments_hiv.sh

python evaluate_alignments/plot_alignment_diff.py evaluate_alignments/out/hiv20000_alignment_differences.tsv evaluate_alignments/data/hiv/hiv20000_pairs.tsv hiv evaluate_alignments/out/hiv20000_alignment_differences.png
```

For results, see `evaluate_alignments/README.md`.

