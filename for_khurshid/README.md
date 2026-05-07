# Tutorial on how to make sequence logo from taxonium file.

## Data

Download the taxonium jsonl file.

```bash
# get jsonl file
wget https://raw.githubusercontent.com/AngieHinrichs/viral_usher_trees/refs/heads/main/trees/Human_respiratory_syncytial_virus_A/tree.jsonl.gz
# get gbff file
wget https://raw.githubusercontent.com/AngieHinrichs/viral_usher_trees/refs/heads/main/trees/Human_respiratory_syncytial_virus_A/treetime_rerooted_NC_001803.1.gbff
```

The `tree.jsonl.gz` file has the mutation information we need for reconstructing the protein sequences for all the samples.

The `treetime_rerooted_NC_001803.1.gbff` has information about the protein sequences. For the case you need, I just found
the sequence for G protein and copied and pasted it into a fasta file `G_protein.fa`.

*Note: I decided not to use the usher tree file, `optimized.pb.gz`, because I need to some additional steps to get the*
*amino acid translations, which is annoying... I decided to just use the jsonl file instead.*

## Parse the jsonl file

I asked claude to write me a script (with close supervision) to reconstruct the sequences from the jsonl file. Here are the
main steps:

1. Read the jsonl file and use the node_id and parent_id information to reconstruct the tree structure, and store the mutation
information for each node.

2. For each leaf node (sample), walk upward from the node to the root. Reverse the path. Starting at the reference sequence,
from the root to the leaf node, apply the relevant mutations in the region of interest. At the end, I have a list of all
the amino sequences of the region of interest for all the samples.

3. Then it outputs the proportions of amino acids at each position, a fasta file containing the reconstructed sequences,
and a sequence logo plot in png and pdf format.

I used a random `logomaker` package to make the sequence logo plot but you can definitely use the other outputs to make your own
in R.

You can run the script like this:

```bash
python aa_logo_from_jsonl.py -i tree.jsonl.gz -G G -s 157 -e 198 --ref-fasta G_protein.fa  -o ccd
```

The arguments are:

```
-i: Input jsonl file
-G: Gene name as stored in the JSONL (e.g. G, F, N)
-s: Start position of the region of interest (1-based, inclusive)
-e: End position of the region of interest (1-based, inclusive)
--ref-fasta: Fasta file containing the reference sequence for the gene of interest
-o: Output prefix
```

Regarding the concern about indels causing a frameshift, I don’t think it’s a major issue. The CCD is highly conserved, so any
indels occurring upstream are likely to be in-frame. Even if an out-of-frame indel does occur, it would likely be non-viable
and quickly purged by natural selection.

However, if you want to be extra diligent, I would recommed getting all the RSV samples, translate the G proteins, and align them with
a traditional MSA tool then build a sequence logo out of that.




