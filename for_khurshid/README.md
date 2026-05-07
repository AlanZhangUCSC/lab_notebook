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

You might need to install some packages with conda. Chatgpt or claude are excellent tools to ask for help. You can run the
same script for rsvb once you have the jsonl file.

Actually, you might be able to get the rsv b from the UCSC genome browser: 

RSV_A: https://hgdownload.gi.ucsc.edu/hubs/GCF/002/815/475/GCF_002815475.1/UShER_RSV-A/

RSV_B: https://hgdownload.gi.ucsc.edu/hubs/GCF/000/855/545/GCF_000855545.1/UShER_RSV-B/

These are trees hosted by UCSC, which more or less have the same set of samples and may have a different reference from the 
tree that Russ sent you. It shouldn't affect the sequence logo.

You can prolly get the reference information here, which you will need to get the reference G protein sequence: 

https://genome.ucsc.edu/cgi-bin/hgPhyloPlace?hgsid=3967298784_oqFvwAaz7cYVaMWHKwaLHljjYq7c&hgpp_org=rsv


Regarding the concern about indels causing a frameshift, I don’t think it’s a major issue. The CCD is highly conserved, so any
indels occurring upstream are likely to be in-frame. Even if an out-of-frame indel does occur, it would likely be non-viable
and quickly purged by natural selection.

However, if you want to be extra diligent, I would recommed getting all the RSV samples, translate the G proteins, and align them with
a traditional MSA tool then build a sequence logo out of that.




