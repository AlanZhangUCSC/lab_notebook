# Lab notebook

This notebook will track the progress of my work in the lab.

## 2025

[Jan](#1302025) &nbsp; [Feb](#212025) &nbsp; [Mar](#392025) &nbsp; [Sep](#9192025) &nbsp; [Oct](#1012025) &nbsp;

[Nov](#1132025) &nbsp; [Dec](#1212025) &nbsp; 

## 2026

[Mar](#3302026) &nbsp; [Apr](#4172026) &nbsp; [May](#512026) &nbsp;


## Projects

[Population Populus PanMAT](#population-populus-panmat)

[Land plant-wide chloroplast PanMAN](#land-plant-wide-chloroplast-panman)

## 1/30/2025

### panMANIA

~~To build the panMAN using Dockerfile, see [panmania/build_panman/](panmania/build_panman/) directory.~~
To build the panMAN see notes on [2/2/2025](#222025).

PanMAN built using Dockerfile is broken. Currently using Summit's docker image.

## 1/31/2025

### panMANIA

Today I want to start by examining the alignments of sequences on a panMAN. I will write scripts that would compare the
global alignments of sequences from the panMAN to the alignments of sequences using mafft. I expect to see that the
mafft alignments should be better than from the panMAN because panMAN alignments are from MSA of all the sequences while
mafft only compares two sequences at a time. Big differences, however, indicate particular problematic alignments,
likely on the block level.

Inside evaluate_alignments, `parse_newick_to_pairs.py` is a script that would parse the newick tree and output the
parent-child pairs in depth first order. The sequences of each parent-child pair will be compared.

Following commands will generate the pairs, run the alignment comparison, and plot the alignment differences.

```
python parse_newick_to_pairs.py data/hiv/hiv_panman.newick > data/hiv/hiv20000_pairs.tsv

sbatch evaluate_alignments/evaluate_alignments_hiv.sh

python evaluate_alignments/plot_alignment_diff.py \
  evaluate_alignments/out/hiv20000_alignment_differences.tsv \
  evaluate_alignments/data/hiv/hiv20000_pairs.tsv \
  hiv \
  evaluate_alignments/out/hiv20000_alignment_differences.png
```

For results, see [panmania/evaluate_alignments](panmania/evaluate_alignments).

## 2/1/2025

### panMANIA

I realized that there are newer versions of panmans availble. I'm going to build the latest version of panman and look
at the alignments comparison again for newer panmans.

## 2/2/2025

### panMANIA

Build using Sumit's Docker image.

Set up
```
cd /private/groups/corbettlab/alan/panmania/
git clone https://github.com/TurakhiaLab/panman.git

cd panman
```

Pull the PanMAN docker image from DockerHub
```
docker pull swalia14/panman:latest
```

Run the docker container and mount the panman directory
```
docker run -it -v .:/local_panman swalia14/panman:latest
```

Edit the CMakeLists.txt by adding the following lines between `TARGET_LINK_LIBRARIES(...)` and
`target_include_directories(...)`

```
add_custom_command(TARGET panmanUtils POST_BUILD
    COMMAND mkdir -p ${CMAKE_BINARY_DIR}/lib
    COMMAND cp /usr/local/lib/libcapnp*.so* ${CMAKE_BINARY_DIR}/lib/
    COMMAND cp /usr/local/lib/libkj*.so* ${CMAKE_BINARY_DIR}/lib/
    COMMAND cp /usr/lib/x86_64-linux-gnu/libboost*.so* ${CMAKE_BINARY_DIR}/lib/
    COMMAND cp /usr/lib/x86_64-linux-gnu/libprotobuf*.so* ${CMAKE_BINARY_DIR}/lib/
    COMMAND cp ${CMAKE_BINARY_DIR}/tbb_cmake_build/tbb_cmake_build_subdir_release/libtbb*.so* ${CMAKE_BINARY_DIR}/lib/
)
```

Run the install script
```
./install/installationUbuntu.sh
```

Add this line to .bashrc and now panmanUtils can be ran outside the docker container
```
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/private/groups/corbettlab/alan/panmania/panman/build/lib
```

Will develop panMANIA using this container for now.

## 2/3/2025

### To-do

- [ ] Make a figure that shows the correlation between the pseudo-chaining score and the alignment score.
- [ ] Fix `position_info_by_haplotype` bug in `heuristic_v6.py`.
- [ ] Complete the alignment quality assessment for new panMANs.

### Logistics and general work

I have started writing the specific aim section of my thesis proposal. For aim 1, other than figures for accuracy
assessment, I think I can also have figure that shows the correlation between the pseudo-chaining score and the
alignment score.

### panMAMA

When I'm in the progress of formalizing the method for calling consensus, I noticed a pretty significant bug in
`heuristic_v6.py` that gives inconsistent `position_info_by_haplotype` for the same haplotype.
`position_info_by_haplotype` is a dictionary that stores other haplotypes' position info at a specific aligned position.

Will need to fix this before I can continue with the formalization.

### panMANIA

The panMAN gives `Exceeded message traversal limit` error for SARS and HIV trees. The reason is that the default message
traversal limit is too small.

To fix this, replace `uint64_t traversalLimitInWords = 8 * 1024 * 1024 * 128` in
`/home/capnproto-c++-1.0.2/src/capnp/message.h` with `uint64_t traversalLimitInWords = 8 * 1024 * 1024 * 256`.

Then run command below to rebuild capnproto.
```
make clean && ./configure && make -j && make install
```

\
As a continuation of `2/1/2025`, I need to assess the alignment quality of the new panMANs. Writing and split the fasta
files for nodes on the panMANs are currently taking too long for `SARS_20000` and `HIV20000`.

## 2/4/2025

### To-do

- [ ] Complete the alignment quality assessment for new panMANs.
- [ ] Make a figure that illustrates the effectiveness of pseudo-chaining scores at discerning between sequences.
- [x] ~~Make a figure that shows the correlation between the pseudo-chaining score and the alignment score.~~ *need to* 
*find a better way to illustrate this*
- [x] ~~Fix `position_info_by_haplotype` bug in `heuristic_v6.py`.~~ *Was redefining `pos` in debug print statements.*

### panMAMA

To make the pseudo-chaining score vs sequence similarity figure (refer to `2/3/2025`), I'm going to write a script that
would give all 150-bp kmers from a random node (preferably containing no ambiguous bases, `OM857280.1`) and generate 50
mutated kmers from each kmer, each containg 1 - 50 mutations. Then I will calculate align them using minimap2 and
pseudo-chaining and compare their scores. Scripts and data are placed in
`/private/groups/corbettlab/alan/lab_notebook/panmama/pseudo_chaining-vs-sequence_similarity`.

The final figure is saved in `panmama/pseudo_chaining-vs-sequence_similarity/pseudo_chaining_vs_seq_similarity.png`. It
doesn't look very good... A better to illustrate the effectiveness of pseudochaining is to show how it's able to discern
between haplotypes, since we only care about if it's sensitive enough to tell different sequences apart.

## 2/5/2025

### To-do

- [ ] Finish writing the first draft of the specific aim section of my thesis proposal.
- [ ] Really do finish formalizing the method for calling consensus.

### Logistics and general work

#### Meeting with Russ
I talked a little about panMANIA then the consensus calling step with Russ. Russ gave a really good sugesstion on
assuming there no back mutations. See notes below in panMAMA section.

### panMANIA

I just plotted the alignment differences for the new SARS and RSV panMANs (HIV is skipped because it's taking too long.
Seems like the aligned sequences in the new HIV panMAN are very long). The new SARS panMAN looks pretty identical to the
old SARs panMAN but **the new RSV panMAN looks **worse** than the old one**. See
[panmania/evaluate_alignments](panmania/evaluate_alignments) for more details.

Maybe it is worth to also fix the alignments in panMAN for Aim 2.

Notes below need to be beautified and executed. This is a place holder so I don't forget

*I have an idea to identify problematic block states. Will try to implement and see how it goes.*
- *Record positions of mismatches in the panMAN alignment. Merge positions into ranges. Identify nearby ranges with* 
  *identical lengths of mismatches.*
- *Nearby ranges with identical lengths of mismatches are likely to be problematic block states.*
- *Try to align those two ranges together and see how well they align.*


### panMAMA

Notes below need to be beautified and executed. This is a place holder so I don't forget

*Russ has a good idea to call consensus*
- *Assume there's not back mutation*
- *If there is a cryptic mutation, assign that cryptic mutation heuristically*
  - *If there are haplotypes where that cryptic mutation is assigned to that haplotype only, assign it.*
  - *If not use maximum likelihood.*
- *If no cryptic mutation, assign reference*

## 2/6/2025

### To-do

- [ ] Really do finish formalizing the method for calling consensus.
- [x] ~~Finish writing the first draft of the specific aim section of my thesis proposal.~~
      *See [thesis_proposal/](thesis_proposal/)directory for the draft*

### Logistics and general work

I finished the first draft of my specific aims section for my thesis proposal. Aim 3 still needs a lot of work. See
[thesis_proposal/](thesis_proposal/) section for the first draft.

### panMANIA

With the limited time I have today, I will write up the script for identifying problematic block states that cause large
chunks of misalignment. See [panmania/evaluate_alignments/](panmania/evaluate_alignments/) for details

## 2/7/2025

### To-do

- [ ] Really do finish formalizing the method for calling consensus.
- [ ] Quantify the prevalence of problematic block states

### panMANIA

Working on quantifying the prevalence of problematic block states

## 2/9/2025

### To-do

- [x] ~~Really do finish formalizing the method for calling consensus.~~
- [x] ~~Quantify the prevalence of problematic block states~~

### panMAMA

Today I will try to finish formalizing consensus calling...

Finally finished formalizing consensus calling. See [panmama/consensus_calling/](panmama/consensus_calling/) for sample
scripts.

Pseudocode:
```
for line in vcf_all:
  possible_alleles = all alleles with depth > abundance - tolerance

  split reads into groups by what other haplotypes they are assigned to

  if an allele-determining group (reads assigned to only the current haplotype) exists:
    assign the allele aligned in the allele-determining group
  else
    iterate through the groups of reads and calculate the likelihoods of the current haplotype having each possible
    allele
    
    split the groups of reads futher into groups by the number of haplotypes they are assigned to.
    
    within each group split by the number of haplotypes they are assigned to, select the subgroup with the highest
    likelihood ratio. This is the group's representative likelihood ratio

    accumulate the likelihood ratio for each possible allele across the groups

    if highest accumulated likelihood ratio - second highest accumluated likelihood ratio > threshold:
      assign allele with the highest accumulate likelihood ratio
    else:
      assign reference allele
```

### panMANIA


I fininally finished quantifying the prevalence of problematic block states (see
[panmania/evaluate_alignments/](panmania/evaluate_alignments/#quantify-the-prevalence-of-problematic-block-states)).
It looks like it's definitely worth fixing the misaligned blocks or representing them differently (a different layer of
block coordinate?). Will talk to Russ about it.

Fixing this should also substantially improve the runtime of `panmap` and `panMAMA`.

## 2/11/2025

### To-do

- [ ] Characterize errors in consensus calling in panMAMA
- [x] ~~Develop a formal pipeline for evaluating panMAN internal alignments~~

### panMAMA

Today I will manually look through the output of consensus calling and characterize the errors. See
[panmama/consensus_calling/](panmama/consensus_calling/)
for detail.

Idea: perhaps I can also calculate the average sequence similarity of the reference and reads supporting each allele. 

I think I'll also assess and plot the threshold for assigning alleles.

### panMANIA

Steps for developing a formal pipeline for testing panMAN internal alignments.
  
1. Add a function in panMAN to output the ranges of all blocks

2. Add a function in panMAN to randomly selectly N pairs of nodes and output their aligned sequences and unaligned
   sequences

3. Run my other scripts detailed in [panmania/evaluate_alignments/](panmania/evaluate_alignments/)

## 2/12/2025

### To-do

- [x] ~~Characterize errors in consensus calling in panMAMA~~

### panMAMA
See [panmama/consensus_calling/](panmama/consensus_calling/) for details.

## 2/13/2025

### To-do

- [ ] Improve panMAMA consensus calling

### panMAMA

I fixed several issues in the consensus calling step.

1. Implemented more clever way to handle when there is only one read group or all read groups support the same allele.

2. Partially incorporated read to reference sequence similairty into consensus calling. Currently support raw counts but
can take position into account as well (counting reads mapped to the same position only once).

3. Still need to implement the 2nd round of ambiguity resolution.

4. Need to implement a way to handle indels.

## 2/19/2025

### To-do

- [ ] Finalize panMAMA consensus calling

### panMAMA

For the past few days, I've been improving panMAMA consensus calling. I think I've improved it enough to start
finalizing the method.

Some of the trivial things to finish are:

- implement indel handling

- Add a likelihood ratio test for cases where more than one haplotype could possibly have the cryptic mutation

- Perhaps output a file that lists all possible variants and their likelihoods in different haplotypes

  Columns:
  
  1. #Type (cryptic, non-cryptic)
  2. REF
  3. ALT
  4. INFO (in the format Chrom:Pos:LikelihoodRatio)


~~I think I will also plot a PR curve for the consensus calling.~~ nvm I tried and it looked ugly.

## 2/26/2025

### To-do

- [ ] Finalize panMAMA consensus calling (see notes on [2/19/2025](#2192025))

- [ ] Optimize panMAMA runtime

- [ ] Evaluate panMAMA using real data

- [ ] Implement homopolymer compressed syncmers

- [ ] Implement sorted and merged syncmer/kminmer index if possible

- [x] ~~**Make a figure of read assignments for Russ**~~ *See
      [panmama/illustrating_figures_andor_diagrams/](panmama/illustrating_figures_andor_diagrams/)*

### panMAMA

I'm going to do some final optimization of panMAMA runtime. This will be done on silverbullet for efficiency.

Keeping the command and options here so I don't have to re-write them.

```
panmap example.panman example_R1.fastq example_R2.fastq --place-per-read --redo-read-threshold 0 --em-filter-round 2 \
  --remove-threshold 0.01 --rounds-remove 5 --preem-filter-method mbc --save-kminmer-binary-coverage
```

After meeting with Richard Durbin, we got some good ideas to further improve accuracy and speed. We can compress
homopolymers when sketching syncmers, and this could reduce syncmer errors due to homopolymer related sequencing errors.
Additionally, we can implement a sorted and merged syncmer/kminmer index as direct hash lookup, while has a O(1) lookup
time, is much slower than linear lookup of adjacent hash values.

## 2/27/2025

### To-do

- [ ] Optimize panMAMA runtime (high priority)

- [ ] Modify panmap build and src files to read new panMANs (high priority)

- [ ] Evaluate panMAMA using real data (high priority)

- [ ] Finalize panMAMA consensus calling (see notes on [2/19/2025](#2192025)) (high priority)

- [ ] Implement homopolymer compressed syncmers (low priority)

- [ ] Implement sorted and merged syncmer/kminmer index if possible (low priority)

### panMAMA

#### Optimize panMAMA runtime

I'm going to graph the distribution of the types of kminmer changes. Specifically, the number of kminmer changes for a
read and the type of kminmer changes (Whether the change had added a kminmer to the chain, removed a kminmer from the
chain, or didn't change the chain). See [panmama/run_time_optimization/](panmama/run_time_optimization/) for details.

## 2/28/2025

### To-do

- [ ] Optimize panMAMA runtime and memory usage (high priority)

- [ ] Modify panmap build and src files to read new panMANs (high priority)

- [ ] Evaluate panMAMA using real data (high priority)

- [ ] Finalize panMAMA consensus calling (see notes on [2/19/2025](#2192025)) (high priority)

- [ ] Implement homopolymer compressed syncmers (low priority)

- [ ] Implement sorted and merged syncmer/kminmer index if possible (low priority)

### panMAMA

#### Optimize panMAMA runtime and memory usage

I will install and test out a profiler (linaro map) today to see where panMAMA is spending most of its time. See
[panmama/run_time_optimization/](panmama/run_time_optimization/) for more details.


## 3/9/2025

### To-do

- [x] Optimize panMAMA runtime and memory usage (high priority)

- [ ] Modify panmap build and src files to read new panMANs (high priority)

- [ ] Evaluate panMAMA using real data (high priority)

- [ ] Finalize panMAMA consensus calling (see notes on [2/19/2025](#2192025)) (high priority)

- [ ] Implement homopolymer compressed syncmers (low priority)

- [ ] Implement sorted and merged syncmer/kminmer index if possible (low priority)

### panMAMA

I think I've finished optimizing panMAMA runtime and memory usage for the time being. Now, panMAMA scores 1 million
reads against 20,000-sample SARS tree in 5 minutes using 9GB of memory. EM step takes ~14 minutes.

Now to evaluate panMAMA using real data. See [panmama/real_data/](panmama/real_data/) for details.

## 3/18/2025

### To-do

- [ ] Evaluate panMAMA using real data (high priority)

- [ ] Implement coverage filter for panMAMA (high priority)

- [ ] Modify panmap build and src files to read new panMANs (high priority)

### panMAMA

I have been trying different ways to analyze the RSV data that Marc sent me. They are looking a little better now that I
ran reads generated from RSVA priemrs and RSVB primers separately. Assuming Marc's assignment of mixture/isolated is
correct, the false positives from panMAMA are because of low depth over all or high depth of likely contaminants at one
or two loci. This is actually very useful information, as it's probably a good idea to consider coverage in the tool
during demixing (set a filter to ignore haplotypes with coverage less than some threshold, [at least for amplicon
samples]).

See my [notes](https://docs.google.com/spreadsheets/d/14gnGfWoBKk4Jzjnv00v8XbslElWe5N_IhrAk8pv9s_k/edit?usp=sharing) for
details.

There are some cases I'm not sure when the samples would have 40-50% coverage with low depth in both RSV A and RSV B, I
have these highlighted in bright yellow in the notes column. I also found two instances of "false negatives" where I had
assigned them to be unmixed while you assigned them to be mixed:

```
RSV00174
RSV00220
```

After manual inspection of these two samples on the IGV, I actually found them to be unmixed (their alignments on IGV
look like the false positives where they have high depth at only one or two loci). They are also highlighted in bright
yellow.

It seems like a good idea for me to also incorporate kminmer overlap into prefiltering the nodes before EM (at least for
amplicon samples). Most of the false positives are also cases where contaminant read stacks map to one or two loci on
one of the RSV types.

Will stare at the google sheet for a while to see if I can find any other patterns.

## 3/23/2025

### To-do

- [ ] Evaluate panMAMA using real data (high priority)

- [ ] Implement coverage filter for panMAMA (high priority)

- [ ] Look over Faith's code and add more plots (high priority)

- [ ] Modify panmap build and src files to read new panMANs (high priority)

### panMAMA

I discussed with Marc our reasonings for determining whether the RSV samples are mixed or unmixed. My logic is that if
there are truly both RSVA and RSVB present in sample, then the reads generated from RSVA primers should map well to RSVA
ref and cover most of the genome, and the reads generated from RSVB primers should map well to RSVB ref and cover most
of the genome. On the other hand, samples where reads only map well to a few regions on the genome likely contain
contaminants that get amplified in the PCR process. Marc, on the other hand, offers an alternative interpretation that
we should believe what the data is telling us:

>The target genomes in any of the samples (not just the mixed/co-infected ones) may not always be full-length (due to
>sample degradation over time(?)), and may not always be present in equimolar amounts (along the genome for one subtype).
>For me [Marc], I find it much harder to understand how RSVA or RSVB contaminants could be introduced in an essentially
>stochastic way among 192 samples all collected on different dates in a variety of locations and processed on different
>dates.  You might want to discuss these alternatives (and others) with Russ.

I will discuss this with Russ.

### panMANIA

Thanks to Faith Okamoto's amazing rotation project, I have some really good preliminary data and code to work with. I
copied over Faith's google docs to my google drive in case they some how disappear from the internet.

I will generate some additional figures as preliminary data for my advancement.

See [panmania/imputation/](panmania/imputation/) for details.


## 3/23/2025

### To-do

- [x] ~~Look over Faith's code and add more plots (high priority)~~ *See [panmania/imputation/](panmania/imputation/)*
      *for details*

- [ ] Evaluate panMAMA using real data (high priority)

- [ ] Implement coverage filter for panMAMA (high priority)

- [ ] Modify panmap build and src files to read new panMANs (high priority)

### panMANIA

#### Imputation

I've generated some additional plots for my advancement. See [panmania/imputation/](panmania/imputation/) for details.
In a nutshell, I made 2 figures, 1 to show the prevalence of missing data in different panMANs in average number of Ns
per genome, and 1 to show the distribution of the sizemissing data neighborhoods. Missing data neighorhoods are
coordinates with missing data that occur in contiguous nodes.

#### Annotation

I will also begin exploring existing annotation tools to get a sense of how to start with panMAN annotation.

See [panmania/annotation/](panmania/annotation/) for details.

## 9/19/2025

It's been a while since I've added to this notebook. I think I should probably
start again.

Today, I am going to continue to explore ways to address the memory problem when
running panMAMA with the 8M tree.

I tried out sorting the reads using a debruijn graph yesterday. The idea is that reads whose kminmers are close to each
other on the graph will tend to be updated together and have similar updates. This way, it will both improve the cache
locality during place and allows clever tricks to reduce the memory of score-annotated index.

Memory can be reduced by having a 32-bit startIndex, 16-bit scoreDelta, and 64-bit trailingScoresDelta that stores the
differences in scoreDelta with respect to the 16-bit scoreDelta in the trailing 16 sorted reads, using 4 bits per reads.
This offers ~6.58 bits per read (in best case scenario) as opposed to 64-bit per read I had previously, with 32-bit
readIndex and 32-bit scoreDelta.

Here is a table of estimated memory reduction under different threads using 1 million SARS reads and 20K-SARS tree. The
old method, which uses 64 bits per read, will use ~1,459 Mb. I've also attached the time for **placement** for
additional information.

*Note: 20K-SARS tree doesn't have any memory issues. I'm using it here to measure relatively how much memory will be*
*reduce. The problem is the 8M-SARS tree. 8M-SARS tree currently is F'ed-up and being fixed by Sumit.*

| Threads | Time (s) | Memory (Mb) | Fraction of original | Median Thread Time | Min Thread Time | Max Thread Time | (Max - Min) Thread Time | 
| :------ | :------- | :---------- | :------------------- | :----------------- | :-------------- | :-------------- | :---------------------- |
| 1 | 183 | 241 | 0.164 | 1 | 1 | 1 | 0 |
| 2 | 155 | 241 | 0.165 | 155 | 132 | 155 | 22 | 
| 4 | 77 | 241 | 0.164 | 46 | 29 | 77 | 47 |
| 8 | 89 | 242 | 0.165 | 33 | 21 | 89 | 56 | 
| 16 | 52 | 242 | 0.165 | 16 | 13 | 52 | 39 | 
| 32 | 65 | 242 | 0.165 | 10 | 8 | 65 | 57 |
| 64 | 39 | 241 | 0.164 | 6 | 5 | 40 | 35 |

It seems like we will get ~6x smaller memory from the read scores, which is not bad. However, dividing up the reads
evenly after sorting has caused some unbalanced workload between the threads. Table below provides runtime information
with shuffled reads for reference (the original version before sorting).

| Threads | Time (s) | Median Thread Time | Min Thread Time | Max Thread Time | (Max - Min) Thread Time | 
| :------ | :------- | :----------------- | :-------------- | :-------------- | :---------------------- |
| 1 | 265 | 265 | 265 | 265 | 0 |
| 2 | 110 | 110 | 104 | 110 | 5 |
| 4 | 58 | 56 | 54 | 57 | 3 |
| 8 | 33 | 30 | 29 | 33 | 4
| 16 | 19 | 18 | 17 | 19 | 3 |
| 32 | 14 | 11 | 11 | 14 | 3 |

The unbalanced workflow between threads increases the runtime too much for multiple threads. 7X reduced memory is good
but might not be worth it in smaller trees, such as the SARS-20K, RSV-4K, and HIV-20K. After talking with Russ, we
decided to implement a low-memory mode, which will be recommended to turn on for the 8M-SARS tree but not for the
smaller trees.

Implementation in progress: `5368e898b7da68a3b7f1c358df2b484c1141fc8f`

## 9/22/2025

It's probably a good idea to also store the number of tailing deltas stored to avoid redundant operations when getting
read scores for EM. I can also add a `uint16_t` in the struct without adding memory used because the 2 bytes of memory
was part of the padding in the struct that will now be used for the new `uint64_t`.

I didn't take into account the cases when there are reads without any updates in the trailing delta. Now trailing 
scoreDeltas are within -7 to +7 to the first scoreDelta (was -8 to +7 in the implementation from [yesterday](#9192025)).
That saves me one integer, -8, to represent skipped scoreDelta in a group. This means the memory improve will be
slightly less than the estimation from [9/19/2025](#9192025).

```c++
struct readScoreDeltaLowMemory {
  uint64_t trailingDelta = 0;
  uint32_t readIndex;
  uint16_t numTrailing = 0;
  int16_t  scoreDelta;

  void encodeTrailingDelta(int16_t scoreDeltaDiff, uint32_t indexToEncode);
  int16_t decodeTrailingDelta(uint32_t offset);
};
```

Now, to score 1 million SARS reads against the SARS-20K tree panMAMA uses ~~~241~~ ~250 Mb of memory, which is ~5.9
times smaller than the original approach to score read scoreDeltas individually.

This is implemented in commit: `928a7d15a53afcfa8296b38a554d3c457c0c4401`

## 9/23/2025

Seems like Sumit has fixed the SARS-8M panMAN. I will take a look at it today.

Hmmm... I got this error when trying to build the MGSR index for the new SARS-8M panMAN

```
stepping right over more than one block from 0 to -1
```

After many debug lines and long time of waiting, I found that I had a hidden bug that seemed tohave never cause a
problem in RSV-4K, SARS-20K, HIV-20K, and the old SARS-8M.

After fixing the bug, I rebuilt the RSV-4K and SARS-20K index and compared them to the old index before I fixed the bug.
Indeed they are identical. Just to make sure, I will compare them to brute force again tomorrow.

MGSR index finally built for the new SARS-8M tree. I will take a look at it tomorrow.

## 9/24/2025

Last night, I scored 100K reads against the SARS-8M tree using both the low-memory mode and the normal mode.

Using 8 threads:

Normal mode took 7352 secs (122.5 mins), using max res 220 Gbs

Low-mem mode took 7115 secs (118.5 mins), using max res 81 Gbs

There are 82,640 nodes with significant kminmer overlap coefficient (using the default settin to include the top 1000
rank).

Inside silverbullet:/scratch1/alan/goodmap/panmap/build, the commands used are:

```
./bin/panmap ../panmans/sars_8M.panman \
  ../test_data/sars/rep1/sars20000_5hap-a_100000_rep1_R1.fastq \
  ../test_data/sars/rep1 sars20000_5hap-a_100000_rep1_R2.fastq \
  -m sars_8M.new.pmai --cpus 8
```

```
./bin/panmap ../panmans/sars_8M.panman \
  ../test_data/sars/rep1/sars20000_5hap-a_100000_rep1_R1.fastq \
  ../test_data/sars/rep1 sars20000_5hap-a_100000_rep1_R2.fastq \
  -m sars_8M.new.pmai --cpus 8 --low-memory
```

This doesn't look very good tbh... I scheduled a meeting on 9/29/2025 with Yatish, Russ, and Alex to go over potential
strategies. I will make a slide deck to go over the panMAMA program step-by-step.


### Meanwhile, I think I will explore possible ways to better subset probable nodes.

Pranav's scoring function looks quite interesting. I already store the maximum parsimony score of each read, and it's
trivial to also store the EPP (number of genomes that have the maximum parsimony with the read).

For the sake of testing, I only scored the nodes that are already selected using the kminmer overlap coefficient method.
Since the overlap coefficient method is quite sensitive, it will pick up a lot of nodes, and I want to see how well the
node scores can filter out some of the noises. Similar to Pranav's scoring function, each read's score depends on its 
maximum parsimony across the entire tree and how many nodes share the paximum parsimony score.

Essentially the score for a specific read i is

`Si = (Ri.kminmers.size() - Ri.maxScore + 1) * pow(Ri.epp, 2)`

And the score for a specific node j is the sum of read scores of reads that maximally map to the node, or in pseudocode:

```
curNodeScores
nodeScore = 0
for (i = 0 ... len(reads)):
  if (curNodeScores[i] == reads[i].maxScore && reads[i].maxScore > 0):
      nodeScore += (reads[i].kminmers.size() - reads[i].maxScore + 1) * pow(reads[i].epp, 2)
```

I haven't implemented the regularization by fraction of sites covered but it already looks pretty good as it is.

I made some minor changes in panmap (see commit ~~`2ae3b7c85528b05b51b80955437e823218b50c28`~~
~~`97aea060e789dafde7cdad4e855d0ea7e6f3c10a`~~ `fff14e03f12b609f1877d845299862cb9b71b12d`) to output some stats on the
node scores, and wrote a bash script
[panmama/node_scores/run_panmap_for_test_node_scores.sh](panmama/node_scores/run_panmap_for_test_node_scores.sh)
to run it on some sample data

```bash
# In termial
reps=(rep1 rep2 rep3 rep4 rep5 rep6 rep7 rep8 rep9 rep10)
for rep in "${reps[@]}"; do
  bash run_panmap_for_test_node_scores.sh /scratch1/alan/goodmap/panmap/build/bin/panmap /scratch1/alan/goodmap/panmap/panmans/sars_20000_optimized.panman /scratch1/alan/goodmap/panmap/build/sars_20000.pmai /scratch1/alan/goodmap/panmap/test_data/sars/$rep/ &
done
```

I will take a look at the results tomorrow.

## 9/25/2025

The bash script seems to have run correctly, and the output files are all there.

I will write a script to ouptut and plot whether true haplotypes are selected by kminmer overlap coefficients and if
they are, how their node scores rank within the selected nodes.

Actually, I think I will rebuilt and rerun the panmap to output the node scores with more precision. Using commit
~~`97aea060e789dafde7cdad4e855d0ea7e6f3c10a`~~ `fff14e03f12b609f1877d845299862cb9b71b12d`.

Here is the script to output stats on true haplotypes kminmer overlap coefficients and node scores
[panmama/node_scores/getRank.py](panmama/node_scores/getRank.py).

To run it, for example

```bash
python3 getRank.py \
  /scratch1/alan/goodmap/panmap/test_data/sars/panman_outputs/rep1/sars20000_5hap-a_abundance_rep1.tsv \
  /scratch1/alan/goodmap/panmap/test_data/sars/rep1/sars20000_5hap-a_200000_rep1.testScores.txt \
  <(sort -k3,3 -gr /scratch1/alan/goodmap/panmap/test_data/sars/rep1/sars20000_5hap-a_200000_rep1.testScores.txt)
```

The output (after piping it to `column -t`):

```
Denmark/DCGC-636002/2022|OY803981.1|2022-12-05                                      0.5   1.0     0  393  342.41096  0    1
Germany/Molecular_surveillance_of_SARS-CoV-2_in_Germany/2021|OV351922.1|2021-09-17  0.2   1.0     0  393  13.38391   3    4
Scotland/QEUH-3C94804/2022|OW532842.1|2022-03-27                                    0.15  1.0     0  393  47.15048   2    3
USA/MA-CDCBI-CRSP_AP7MZ6THJQWE6Q7J/2023|OQ727617.1|2023-03-18                       0.1   1.0     0  393  121.36997  1    2
England/OXON-AD71F/2020|OX589494.1|2020-04-04                                       0.05  0.9998  1  396  0.94889    359  578
```

I then wrote a script [panmama/node_scores/run_getRank.sh](panmama/node_scores/run_getRank.sh) to run `getRank.py` for my samples.

```bash
reps=(rep1 rep2 rep3 rep4 rep5 rep6 rep7 rep8 rep9 rep10)
for rep in "${reps[@]}"; do
  bash run_getRank.sh getRank.py /scratch1/alan/goodmap/panmap/test_data/sars/$rep/ /scratch1/alan/goodmap/panmap/test_data/sars/panman_outputs/$rep &
done
```

I found that all of my samples with 80K reads actually contain different haplotypes from the true haplotypes? I probably
forgot to delete them before re-simulation from a long time ago. Anyway, I just removed them.

I think I will actually output kminmer overlap coefficient and node scores of all the nodes on the SARS-20K tree. Using
commit `fff14e03f12b609f1877d845299862cb9b71b12d` instead.

After manual inspection, it seems like Prana's node scoring scheme is pretty promising if I were to pre-select based on
kminmer overlap coefficients. I do think it's worth it regularize the node scores using kminmer coverage, which I think
will improve both sensitivity and specificity and might free us from pre-selection by kminmer overlap coefficients
entirely

Although, I found something weird... I was expecting that, at least for high coverage, in samples with only a single haplotype, 
the correct haplotype should have the highest node score. However, for `rep10` of my samples, the correct haplotype does
not have the highest node score.

```
==> sars20000_1hap-a_100000_rep10_rank_stats.tsv <==
England/MILK-9A8502/2020|OA972423.1|2020-09-02  1.0  1.0  0  7  34.7493051276  2  3  2  3

==> sars20000_1hap-a_10000_rep10_rank_stats.tsv <==
England/MILK-9A8502/2020|OA972423.1|2020-09-02  1.0  1.0  0  1  3.1636666008  1  2  1  2

==> sars20000_1hap-a_200000_rep10_rank_stats.tsv <==
England/MILK-9A8502/2020|OA972423.1|2020-09-02  1.0  1.0  0  9  65.2952724646  3  4  3  4

==> sars20000_1hap-a_20000_rep10_rank_stats.tsv <==
England/MILK-9A8502/2020|OA972423.1|2020-09-02  1.0  1.0  0  1  5.9588101007  1  2  1  2

==> sars20000_1hap-a_2000_rep10_rank_stats.tsv <==
England/MILK-9A8502/2020|OA972423.1|2020-09-02  1.0  0.998377611  2  3  0.6550883658  38  39  1  2

==> sars20000_1hap-a_40000_rep10_rank_stats.tsv <==
England/MILK-9A8502/2020|OA972423.1|2020-09-02  1.0  1.0  0  1  12.4419547827  1  2  1  2
```

I might have made some sense if only one sample doesn't have the correct haplotype ranked highest, due to sheer luck 
that a sequencing error caused a read or two to map better to the wrong node, but all the samples are like this.

Running `for file in sars20000_1hap-a_*_rep10.testScores.txt; do  echo ">${file}"; sort -k 3,3 -gr $file | head -n 5; done`

```
>sars20000_1hap-a_100000_rep10.testScores.txt
  England/CAMC-1172A52/2021|OD972409.1|2021-01-25 0.9943170286 35.2495987447
  England/QEUH-9CFC6A/2020|OA994276.1|2020-09-15 0.9967539055 35.0967812923
  England/MILK-9A8502/2020|OA972423.1|2020-09-02 1 34.7493051276
  England/PHEC-14D2C4/2020|OX658315.1|2020-11-12 0.9939135727 34.7215691619
  England/QEUH-B419FE/2020|OC996183.1|2020-11-08 0.9957412290 34.3945398909
>sars20000_1hap-a_10000_rep10.testScores.txt
  England/QEUH-9FFDAA/2020|OA998403.1|2020-09-30 0.9979716024 4.1104912017
  England/MILK-9A8502/2020|OA972423.1|2020-09-02 1 3.1636666008
  England/CAMC-1172A52/2021|OD972409.1|2021-01-25 0.9912725797 3.1612628447
  England/PHEC-14D2C4/2020|OX658315.1|2020-11-12 0.9939135727 3.1610533713
  England/QEUH-96E052/2020|OA967653.1|2020-08-18 0.9989858012 3.1107656447q
>sars20000_1hap-a_200000_rep10.testScores.txt
  England/CAMC-1172A52/2021|OD972409.1|2021-01-25 0.9943170286 66.2484996600
  England/QEUH-9CFC6A/2020|OA994276.1|2020-09-15 0.9997971191 66.0141781118
  England/PHEC-14D2C4/2020|OX658315.1|2020-11-12 0.9951308582 65.7404822238
  England/MILK-9A8502/2020|OA972423.1|2020-09-02 1.0000000000 65.2952724646
  England/QEUH-B419FE/2020|OC996183.1|2020-11-08 0.9959440276 65.1060260050
>sars20000_1hap-a_20000_rep10.testScores.txt
  England/QEUH-9CFC6A/2020|OA994276.1|2020-09-15 0.9975654291 6.8180396934
  England/MILK-9A8502/2020|OA972423.1|2020-09-02 1 5.9588101007
  England/CAMC-1172A52/2021|OD972409.1|2021-01-25 0.9912725797 5.9534744836
  England/PHEC-14D2C4/2020|OX658315.1|2020-11-12 0.9939135727 5.9525848546
  England/ALDP-B24FC2/2020|OB994572.1|2020-11-04 0.9983779400 5.9253104840
>sars20000_1hap-a_2000_rep10.testScores.txt
  USA/NV-CDC-QDX26257699/2021|OK250130.1|2021-06-27 0.9629254457 1.1328065150
  England/MILK-B393BF/2020|OC996628.1|2020-11-01 0.9697400487 1.1200297175
  SouthAfrica/NHLS-UCT-GS-D051/2021|OM765744.1|2021-07-08 0.9484094617 1.0110808329
  Denmark/DCGC-641260/2023|OY797037.1|2023-02-01 0.9163256956 1.0102208789
  Germany/IMS-10116-CVDP-DFD838F5-4D85-4B5B-AE67-0724E1E8F1FA/2021|OU078334.1|2021-02-03 0.9691495839 1.0089380969
>sars20000_1hap-a_40000_rep10.testScores.txt
  England/CAMC-1172A52/2021|OD972409.1|2021-01-25 0.9924903592 12.7652863894
  England/MILK-9A8502/2020|OA972423.1|2020-09-02 1 12.4419547827
  England/PHEC-14D2C4/2020|OX658315.1|2020-11-12 0.9939135727 12.4300685190
  node_6130 0.9989860069 12.2138878733
  England/OXON-F888DC/2020|OY954521.1|2020-10-06 0.9939418417 12.2138875445
```

I will use the 10K read sample to investigate the issue. `England/QEUH-9FFDAA/2020|OA998403.1|2020-09-30` has a much
higher node score than `England/MILK-9A8502/2020|OA972423.1|2020-09-02`.

I modified `panmap` to output some info on the read scores that contribute to each of the two haplotypes.

I found that the read below, contributes 1.0 to `England/QEUH-9FFDAA/2020|OA998403.1|2020-09-30` but does not contribute
to `England/MILK-9A8502/2020|OA972423.1|2020-09-02` at all. For simplicity, I will refer to 
`England/QEUH-9FFDAA/2020|OA998403.1|2020-09-30` as `QEUH`, and `England/MILK-9A8502/2020|OA972423.1|2020-09-02` as
`MILK`.

```
@England_MILK_9A8502_2020_OA972423_1_2020_09_02_181_7/2
AGGGTTATGATTTTGGAAGCGCTCTGAAAAACAGCAAGAAGTGCAACGCCAACAATAAGCCATCCGAAAGGGAGTGAGGCTTGTATCGGTATCGTTGCAGTAGCGCGAACAAAATCTGAAGGAGTAGCATCGTTGATTTCACCTTGCTTCA
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF-FFFFFF-FFFFFFF8FFF-FF8FFFFFFFFFFFFFFFFF8FFFFFFFF-FFFFFFFFFF---FFFFFFFFFFFFFF-F-
```

I then used `minimap2` to align this read to both haplotypes.

`minimap2 -a --MD sars_original/England_MILK_9A8502_2020_OA972423_1_2020_09_02.unaligned.fasta ../test_data/sars/exp.fastq`

```
England_MILK_9A8502_2020_OA972423_1_2020_09_02_181_7/2	16	England/MILK-9A8502/2020|OA972423.1|2020-09-02	25436	60	151M	*	0	0	TGAAGCAAGGTGAAATCAACGATGCTACTCCTTCAGATTTTGTTCGCGCTACTGCAACGATACCGATACAAGCCTCACTCCCTTTCGGATGGCTTATTGTTGGCGTTGCACTTCTTGCTGTTTTTCAGAGCGCTTCCAAAATCATAACCCT	-F-FFFFFFFFFFFFFF---FFFFFFFFFF-FFFFFFFF8FFFFFFFFFFFFFFFFF8FF-FFF8FFFFFFF-FFFFFF-FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF	NM:i:1	ms:i:296	AS:i:296	nn:i:0	tp:A:P	cm:i:24	s1:i:129	s2:i:0	de:f:0.0066	MD:Z:19G131	rl:i:0
```

`minimap2 -a --MD sars_original/England_QEUH_9FFDAA_2020_OA998403_1_2020_09_30.unaligned.fasta ../test_data/sars/exp.fastq`

```
England_MILK_9A8502_2020_OA972423_1_2020_09_02_181_7/2	16	England/QEUH-9FFDAA/2020|OA998403.1|2020-09-30	25436	60	151M	*	0	0	TGAAGCAAGGTGAAATCAACGATGCTACTCCTTCAGATTTTGTTCGCGCTACTGCAACGATACCGATACAAGCCTCACTCCCTTTCGGATGGCTTATTGTTGGCGTTGCACTTCTTGCTGTTTTTCAGAGCGCTTCCAAAATCATAACCCT	-F-FFFFFFFFFFFFFF---FFFFFFFFFF-FFFFFFFF8FFFFFFFFFFFFFFFFF8FF-FFF8FFFFFFF-FFFFFF-FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF	NM:i:0	ms:i:302	AS:i:302	nn:i:0	tp:A:P	cm:i:28	s1:i:150	s2:i:0	de:f:0	MD:Z:151	rl:i:0
```

It seems the sequencing error did cause `QEUH` to map better than `MILK`...

Could it just be that the `MILK` genome is so close to its relatives that, more than usually, sequencing errors can
cause reads to map better to other nodes. I will run `panmap` with `no-single` mode to remove reads with obvious
sequencing errors. 

Indeed, `MILK` is ranked first in node scores when `panmap` was run with `no-single` turned on. 

```
$ sort -k3,3 -gr panmap.testScores.txt | column -t
  England/MILK-9A8502/2020|OA972423.1|2020-09-02   1             2.8943746447
  England/CAMC-1172A52/2021|OD972409.1|2021-01-25  0.9912725797  2.8922007995
  England/PHEC-14D2C4/2020|OX658315.1|2020-11-12   0.9939135727  2.8919996823
  England/QEUH-96E052/2020|OA967653.1|2020-08-18   0.9989858012  2.8472196151
  node_6123                                        0.9989768774  2.8472193844
```

I will also try with the 100K and 200K read samples.

100K sample:

```
$ sort -k3,3 -gr panmap.testScores.txt | column -t
  England/MILK-9A8502/2020|OA972423.1|2020-09-02   1             32.8271963293
  England/CAMC-1172A52/2021|OD972409.1|2021-01-25  0.9912725797  32.8044433507
  England/PHEC-14D2C4/2020|OX658315.1|2020-11-12   0.9939135727  32.8012638459
  England/ALDP-A2E3AF/2020|OA991257.1|2020-10-05   0.9955348082  32.2821446972
  England/QEUH-96E052/2020|OA967653.1|2020-08-18   0.9989858012  32.2653056183
```

200K sample:
```
$ sort -k3,3 -gr panmap.testScores.txt | column -t
  England/CAMC-1172A52/2021|OD972409.1|2021-01-25  0.9924903592  62.8683882201
  England/PHEC-14D2C4/2020|OX658315.1|2020-11-12   0.9943193346  62.6107755427
  England/MILK-9A8502/2020|OA972423.1|2020-09-02   1             62.4128674700
  England/QEUH-9CFC6A/2020|OA994276.1|2020-09-15   0.9963481436  62.1861284750
  England/ALDP-A6E6AA/2020|OB982816.1|2020-10-15   0.9975659229  61.7590675371
```

It seems like 200K sample has a little bit more sequencing error because of its large size. Nevertheless, it also seems 
that all the top ranked node scores are `MILK`'s close relatives, as they are all sequenced from England with very close 
dates to each other. I think I will actually also regularize the scores using coverage. I cannot use the existing 
kminmer overlap coefficient because a kminmer is considered covered as long it's in the read sample. For the purose of 
regularization, a kminmer should only be considered covered if it's present in a read that has maximum parsimony at that 
node. Hope that makes sense.....

## 9/26/2025

Today I will come up with and develop a way to calculate the kminmer coverage. I will do a second light-weight tree
traversal to dynamically calculate it.

I will put the kminmer coverage function in `ThreadsManager` for now.

## 9/29/2025

I added a `computeKminmerCoverage()` function in `ThreadsManger` class. This function will only count a k-min-mer being
covered at a node when it has at least one read that has a maximum score with the node. It's quite difficult to engineer
a way to efficiently do this. Nevertheless, it seems like this coverage actually doesn't differ from the overlap
coefficient too much. I think I will hold off this approach for now.

For book-keeping purpose, `computeKminmerCoverage()` is implemented in commit:
`2e18748a0f23ee7ec69b2630985c5a9458e78837`

## 9/30/2025

### Meeting with Yatish, Russ, and Alex was fruitful.

Here are some keypoints from it:
* PanMAMA's runtime is actually not bad and quite similar to WEPP
* Yatish ~~suggested~~ mentioned the option of alternative scoring and deconvolution method.
  * They tried both EM and Freyja and found that Freyja consistently performed better than EM.
  * Alternative deconlution method might allow us to avoid chaining.
* Can pre-filter the probable haplotypes (like the current overlap coeffcient method), then **condense the tree for only
  the selected nodes**
  * Still perform score update but on a much smaller scale.
* Can group reads by their initial mapping position on the root (or even pseudo-root) node, then "build" smaller trees,
each representing a contiguous small section of the genome on which the read group maps to.

Russ and Yatish think that I'd be solid if I can improve my performance by 5-10X.

I will first try to make a subtree that contains only nodes that have mutations relevent to my read sample. I created
a new branch `mgsr-optimization` and will work from there.

I need an efficient way to collapse the k-min-mer mutations if I want to collapse the nodes. To do
this, I think I need the indexed seed insertions, deletions, and substitutions to be stored together instead of

separately. I would also need to sort them by their positions. I will start with that.

## 10/1/2025

### Modifying structure of seed updates in panMAMA index

I finally finished this update. Now I'm storing a single `seedDeltas<SeedDelta>` list per node in the index. Each
`SeedDelta` contains the `seedIndex` and `isDeleted`. This list is sorted by the `startPos` of `seedInfos[seedIndex]`.
This way, I know that if two consecutive `SeedDeltas` have the same `startPos`, they are substitutions. As before, this
allows me to backtrack the reference kminmers using the same index...

### I just realized that I also need an efficient way to collapse the gapMap changes...

This is a bit trickier because my current implementation requires that my `coordDeltas` at node to be processed 
sequentially as indexed. This means that I cannot just sort them to make the collapse easy. I might need to change how 
gaps are represented or updated... 

During building, a `gapRunUpdates` vector was processed to produce the `coordDeltas` vector. I didn't use
`gapRunUpdates` directly for the index because I want to leave the messy computing in indexing so that I can just apply
the `gapMap` changes according to `coordDeltas`. Unlike `coordDeltas`, `gapRunUpdates` doesn't need to be processed in a
specific order. To be able to collapse the nodes, I might have to store `gapRunUpdates` in the index instead. I don't
think the performance will be impacted more than trivially.

## 10/2/2025

### Switching from `coordDeltas` to `gapRunDeltas`

I just implemented the new gapMap update, see [10/1/2025](#1012025). Everything seems to be working correctly.

### Finally I can start working on collapsing the nodes.

I implemented a `MgsrLiteTree` class and a `MgsrLiteNode` class for better organization of the panMAMA index data.

Now `MgsrLiteTree` owns the `seedInfos` instead of `ThreadsManger` or `mgsrPlacer`. And `MgsrLiteNodes` own
their respective indexed seed updates, including `seedDeltas`, `gapRunDeltas`, and `invertedBlocks`, in MgsrLiteTree.

This will set up a better code structure for collaping the tree for optimization.

## 10/3/2025

### Also moved the ownership `scoreDeltas` from `threadsManager` and `mgsrPlace` to their respective nodes.

Same as the changes made on [10/2/2025](#1022025), This will set up a better code structure for collaping the tree for
optimization.

### Investing how we can split the tree

I tried splitting the tree by segment size of `[1000, 5000, 10000, 50000, 100000, 500000]`, and each tree segment is 
collapsed to only include nodes that have at least one k-min-mer mutation in its range. I also give it either a `500` or 
`0.1 * segmentSize` overlap offset so that we can include reads at the boundries of the segment when splitting the reads 
to each segment. 

#### SARS 4K

```
binSize  overlap  numTrees  totalNodes  avgTreeSize  medianTree  minTree  minTreeRange     maxTree  maxTreeRange
1000     500      3446      658944      191          10          0        21500,22500      75       1157000,1158000
1000     100      1915      364805      190          10          0        21600,22600      74       1152900,1153900
5000     500      383       239527      625          123         134      886500,891500    401      1156500,1161500
5000     500      383       239527      625          123         134      886500,891500    401      1156500,1161500
10000    500      182       195632      1074         417         612      1064000,1074000  681      1149500,1159500
10000    1000     192       204953      1067         408         636      1062000,1072000  698      1152000,1162000
50000    500      35        121376      3467         1639        1761     990000,1040000   1944     1138500,1188500
50000    5000     39        130549      3347         1543        21723    1710000,1760000  1930     1125000,1175000
100000   500      18        96420       5356         2878        25280    1691500,1791500  2987     1094500,1194500
100000   10000    20        105363      5268         3246        25542    1710000,1810000  3314     1080000,1180000
500000   500      4         65668       16417        24597       27053    1498500,1998500  24597    999000,1499000
500000   50000    4         65182       16295        21293       29192    1350000,1850000  21293    900000,1400000
```

#### SARS 8M

```
binSize  overlap  numTrees  totalNodes  avgTreeSize  medianTree  minTree  minTreeRange   maxTree  maxTreeRange
1000     500      285       200645202   704018       647251      662308   75500,76500    118214   1500,2500
1000     100      159       110509365   695027       644707      661926   75600,76600    110471   1800,2800
5000     500      32        51089087    1596533      1541060     1355099  45000,50000    417769   0,5000
5000     500      32        51089087    1596533      1541060     1355099  45000,50000    417769   0,5000
10000    500      15        37001282    2466752      2279767     1522025  19000,29000    1099764  0,10000
10000    1000     16        39027983    2439248      2417613     1353098  18000,28000    1090724  0,10000
50000    500      3         17701937    5900645      6053592     6053592  49500,99500    7013652  99000,149000
50000    5000     4         21570598    5392649      6053592     7219395  135000,185000  6053592  90000,140000
100000   500      2         14357690    7178845      7354519     7354519  99500,199500   7003171  0,100000
100000   10000    2         14718565    7359282      7364046     7354519  0,100000       7364046  90000,190000
500000   500      1         8560783     8560783      8560783     8560783  0,500000       8560783  0,500000
500000   50000    1         8560783     8560783      8560783     8560783  0,500000       8560783  0,500000
```

### Investigating where read kminmers can match to on the tree

To roughly estimate this, I get all the global reference positions that each read's k-min-mers can map to, then 
calcualte the read mapping range between the min and max its reference match positions. This measures how far apart on
the reference a read can potentially map to. 

From working directory, `/scratch1/alan/goodmap/panmap/build`, I used 
`../test_data/sars/rep1/sars20000_5hap-a_100000_rep1_R*.fastq` as a test sample.

#### SARS 20K

![sars_20k_read_match_ranges](figures/read_match_positions_plot.png)

#### SARS 8M

![sars_8M_read_match_ranges](figures/8M_read_match_positions_plot.png)

#### RSV 4K

![rsv_4k_read_match_ranges](figures/rsv_read_match_positions_plot_plot.png)

The results show that most of the reads have very close mapping range (which is what we want) on the SARS 20K tree.
However, that doesn't seem like the case for the **SARS 8M** tree and RSV 4K tree.



## 10/6/2025

It seems like a non-trivial amount of reads would not map to a small segument of the SARS 8M tree, so grouping the reads
by their mapping positions then splitting the tree would not work very well.

Therefore, we might have to resolve to reducing the search space early. Today, I will write some scripts to run tests on
how well k-min-mer overlap coefficient is able to pre-filter read samples of various compositions:

<ol>
  <li>Original and mutated haplotype
    <ol>
      <li>All original haplotypes</li>
      <li>All mutated haplotypes</li>
      <li>Mix of original and mutated haplotypes</li>
    </ol>
  </li>
  <li>Sequencing type
    <ol>
      <li>shot-gun</li>
      <li>tailed amplicon</li>
    </ol>
  </li>
  <li>Sequencing depth</li>
</ol>

I will make all combinations of (`#SNPs`, `#Haplotypes`, `%Mutated`, `SeqType`, `Depth`)

`#SNP: 0 1 5 10 20`

`#haplotypes: 1 5 10 20 50 100`

`%Mutated: 0% 20% 50% 70% 100%`

`SeqType: shot-gun tiled-amplicon`

`Depth: 10X 50X 100X 500X 1000X`

### Output haplotypes and simulate SNPS

Added `--dump-sequences` and `--simulate-snps` options to `panmap`. Just for the purpose of testing if correct nodes
will be selected, SNPs are uniformally simulated, e.g. base A can mutate into C, G, and T with equal probabilities. 

Added `--dump-random-nodeIDs` option to output a specified number of nodes. This is for selecting random haplotypes to simulate
the mixed samples

Added `--random-seed` option to specify the seed for the RNG. This is for replicating the results for future reference.

This is implemented in commit `b96f016a79afca83c8e8b45b86f339cd367062e3`

I just finished writing a script, 
[simulate_and_test_overlap_coefficient.sh](panmama/overlap_coefficients/simulate_and_test_overlap_coefficient.sh), to
simulate reads with parameters described above. Currently, for `%Mutated`, I am only simulate either 0% or 100%. I will
deal with mixed mutations tomorrow. See below for how to run the script... I pulled `#SNP` parameter out for simpler
parallelization on `silverbullet`.

```
for num_snp in 0 1 5 10 20; do
  bash simulate_and_test_overlap_coefficient.sh \
    /scratch1/alan/goodmap/panmap/build/bin/panmap \
    /scratch1/alan/goodmap/panmap/panmans/sars_20000_optimized.panman \
    /scratch1/alan/goodmap/panmap/build/sars_20000.pmai \
    /home/alan/tools/SWAMPy/src/simulate_metagenome.py \
    /home/alan/tools/SWAMPy/primer_sets/nimagenV2.bed \
    /home/alan/tools/SWAMPy/ref/MN908947.3.fasta \
    /home/alan/tools/jvarkit/dist/jvarkit.jar \
    /scratch1/alan/lab_notebook/panmama/overlap_coefficients/out \
    delete_intermediate \
    $num_snp > /dev/null 2>&1 &
done
```

This is currently running on `silverbullet`. I will take a look at the results tomorrow.

## 10/7/2025

Here are some stats of the results (for full results, see `panmama/overlap_coefficients/shotgun.results` and
`panmama/overlap_coefficients/amplicon.results`):

### Shotgun sequencing

```
num_hap  num_snp  seq_type  depth  min_rank  max_rank  median_rank
1        0        0         500    0         0         0
1        1        0         500    0         0         0
1        5        0         500    0         0         0
1        10       0         500    0         0         0
1        20       0         500    0         0         0
5        0        0         500    0         89        0
5        1        0         500    0         0         0
5        5        0         500    0         133       12
5        10       0         500    0         21        0
5        20       0         500    0         27        0
10       0        0         500    0         119       0.5
10       1        0         500    0         25        0
10       5        0         500    0         119       0
10       10       0         500    0         27        0
10       20       0         500    0         1199      0
20       0        0         500    0         491       0
20       1        0         500    0         190       0
20       5        0         500    0         321       0
20       10       0         500    0         509       0
20       20       0         500    0         417       51.5
50       0        0         500    0         728       0
50       1        0         500    0         252       0
50       5        0         500    0         1014      0
50       10       0         500    0         508       0
50       20       0         500    0         905       0
100      0        0         500    0         1340      0
100      1        0         500    0         1771      0
100      5        0         500    0         1490      18
100      10       0         500    0         1241      82.5
100      20       0         500    0         1617      7
1        0        0         1000   0         0         0
1        1        0         1000   0         0         0
1        5        0         1000   0         0         0
1        10       0         1000   0         0         0
1        20       0         1000   0         0         0
5        0        0         1000   0         23        0
5        1        0         1000   0         4         0
5        5        0         1000   0         0         0
5        10       0         1000   0         26        0
5        20       0         1000   0         249       0
10       0        0         1000   0         189       0
10       1        0         1000   0         164       0
10       5        0         1000   0         49        0
10       10       0         1000   0         39        0
10       20       0         1000   0         303       37.5
20       0        0         1000   0         75        0
20       1        0         1000   0         62        0
20       5        0         1000   0         297       0
20       10       0         1000   0         65        0
20       20       0         1000   0         481       0
50       0        0         1000   0         172       0
50       1        0         1000   0         570       0
50       5        0         1000   0         292       0
50       10       0         1000   0         203       0
50       20       0         1000   0         519       0
100      0        0         1000   0         482       0
100      1        0         1000   0         695       0
100      5        0         1000   0         593       0
100      10       0         1000   0         937       0
100      20       0         1000   0         1460      0
```

### Tiled amplicon sequencing

```
num_hap  num_snp  seq_type  depth  min_rank  max_rank  median_rank
1        0        1         500    23        23        23
1        1        1         500    0         0         0
1        5        1         500    17        17        17
1        10       1         500    0         0         0
1        20       1         500    47        47        47
5        0        1         500    232       1233      906
5        1        1         500    255       2729      1713
5        5        1         500    56        643       439
5        10       1         500    253       1877      1262
5        20       1         500    58        2315      586
10       0        1         500    62        5714      692
10       1        1         500    149       2062      564.5
10       5        1         500    129       2039      968
10       10       1         500    227       2828      1245.5
10       20       1         500    37        3582      1883.5
20       0        1         500    134       5274      998.5
20       1        1         500    96        2915      1407.5
20       5        1         500    105       5027      1116
20       10       1         500    140       4749      2079
20       20       1         500    43        3918      1274.5
50       0        1         500    66        6687      1686
50       1        1         500    52        6232      2124
50       5        1         500    14        9003      2020.5
50       10       1         500    132       7184      1785
50       20       1         500    101       6261      1905.5
100      0        1         500    111       8461      1982.5
100      1        1         500    165       8915      1910
100      5        1         500    29        5538      1839
100      10       1         500    70        8112      1938
100      20       1         500    128       8900      1708.5
1        0        1         1000   7         7         7
1        1        1         1000   0         0         0
1        5        1         1000   0         0         0
1        10       1         1000   2         2         2
1        20       1         1000   212       212       212
5        0        1         1000   455       6224      1780
5        1        1         1000   139       2042      1268
5        5        1         1000   452       1295      724
5        10       1         1000   260       1910      475
5        20       1         1000   118       877       367
10       0        1         1000   27        1847      1217
10       1        1         1000   70        2403      256.5
10       5        1         1000   72        3154      459.5
10       10       1         1000   267       2365      890
10       20       1         1000   93        2856      1199
20       0        1         1000   117       6584      1576.5
20       1        1         1000   200       2605      1064
20       5        1         1000   339       6378      1577
20       10       1         1000   41        4245      1566.5
20       20       1         1000   37        3512      1590.5
50       0        1         1000   31        4911      1356
50       1        1         1000   310       5102      1636
50       5        1         1000   70        7058      1792
50       10       1         1000   26        7466      1518
50       20       1         1000   4         6954      1507
100      0        1         1000   32        8745      1696
100      1        1         1000   18        6383      1608
100      5        1         1000   284       5937      3426.5
100      10       1         1000   82        4691      2008
100      20       1         1000   92        7673      2611
```

I also just realized that I completely forgot to run it with replicates, so I will make sure to include that in my 
subsequent runs.

I'm **not** very confident about the results, especially for the tiled amplicon sequencing. The median true haplotype 
overlap coefficient ranks are pretty low, indicating that true haplotypes are not being selected very well. While it
performs a little bit better for the shotgun sequencing samples, a LOT of other nodes would also be selected. I think
I need to come up with a better way to reduce the search space.

I started trying WEPP's node scoring approach again with the tiled amplicon data. I tried it manually (running commands
on terminal) and found that it's very good with some haplotypes for also very bad for others. After reading the WEPP
paper on biorxiv, I found that there's one step that I missed for scoring the nodes: from highest scoring node, after
removing the top node, I need to update other node scores by subtracting from them read scores that are also
parsimonious for the node just removed. I will work on this tomorrow.

I suddenly realized that Pranav was simulating and working with tiled amplicon data, meaning that much more reads can be
deduplicated. I then simulated 100K amplicon data and ran my current `panMAMA` implementation on it. As expected,  
runtime is much shorter for the amplicon samples than shotgun sequencing samples. Using 8 threads, normal mode took 
*24.6*  minutes, while low-memory mode took *29.5* minutes. In comparison, 100K shotgun samples took *120* minutes on 8 
threads. Maybe, what we currently have is pretty good afterall? I will discuss this with Russ tomorrow.

Nevertheless, I still need to implement a better way to reduce the search space. One WIP is the WEPP scoring approach.

## 10/8/2025

After discussing with Russ, he also thinks that my current implementation might be pretty good already. I will simulate
1.5 million amplicon and shotgun sequencing reads and run my current `panMAMA` implementation on them. I will run them 
on the SARS 8M tree on 64 threads.

```
Mode        SeqType    Runtime
--------------------------------------------
Normal      Amplicon   2342.3s  (39.04 mins)
LowMem      Amplicon   2801.53s (46.68 mins)
LowMem      Shotgun    16246.8s (4.512 hrs)
```

This is actually not bad compared to WEPP. This should also be further improved after reducing the search space before
the scoring step. Shotgun sequencing samples took much longer than amplicon samples because shotgun reads are much less
collaspible than amplicon reads. 1,500,000 shotgun reads were collapsed to 374,400 unique k-min-mer sets (~25% of
original reads), while 1,401,619 amplicon reads were collapsed to 57,922 unique k-min-mer sets (~4.13% of original
reads).

## 10/9/2025

### Found some major bugs in calculating WEPP node scores...

1. `epp` was not updated correctly. I was only incrementing it when a read score was updated. This caused `epp` to be 
   undercounted.
2. I didn't account for when nodes are identical or indistinguishable by their k-min-mer sets. I was treating them as
   separate nodes when identical nodes should be treated as the same node.

### Some nodes were just not selected by WEPP scores for some reason...

I arbitrarily used `test_data/sars/rep3/sars20000_10hap-a_100000_rep3_R*.fastq` as a test sample. And true haplotypes'
rank for their node scores were:

```
Rank  NodeID                                                                              Score
-------------------------------------------------------------------------------------------------------
3675  USA/CA-CHLA-PLM46323971/2020|MZ722413.1|2020-12-09                                  0.1845768297
5244  USA/IL-CDC-QDX40817725/2022|OP411780.1|2022-08-28                                   0.0764590824
2     England/MILK-344FEB3/2022|OV817379.1|2022-01-26                                     57.8659772952
1     Denmark/DCGC-599763/2022|OY836217.1|2022-10-17                                      85.0571468377
3     USA/CO-CDPHE-41411769/2023|PP031687.1|2023-09-25                                    54.3541009750
4     USA/CO-CDPHE-2007061387/2020|OK659067.1|2020-07-06                                  25.4734981362
7     USA/WA-S20280/2022|ON660803.1|2022-05-20                                            8.3954038314
8     USA/FL-CDC-LC0637436/2022|ON608372.1|2022-05-12                                     8.1637073062
10    Germany/Molecular_surveillance_of_SARS-CoV-2_in_Germany/2021|OV342561.1|2021-09-10  6.2916502660
5     USA/CA-CDC-STM-XAH3WETJM/2022|OP911649.1|2022-11-14                                 19.0574859274
```

It's very weird that the the first 2 haplotypes are ranked very low. I then filtered out the reads for only
`MZ722413.1`. Surprisingly, even with reads that are only from `MZ722413.1`, it still ranked **636** for its node score.
And it's not even selected if I were to progressively exhaust the reads during the selection process. The same happened
for `OP411780.1`, which ranked **1016** for its node score and was also not selected if I were to progressively exhaust
the reads during the selection process. As a positive control, I filtered out the reads for only `OV342561.1`. As
expected, it was ranked **1** for its node score.

I then investigated how the reads map to `MZ722413.1`. Of all reads that were simulated from `MZ722413.1`, after
collapsing them, 10,528 unique reads map parismoniously to `MZ722413.1`, while 1,277 map parsimoniously elsewhere,
meaning that they were sequencing errors. Additionally, among the 1,277 reads that map parsimoniously elsewhere, the
vast majority of them have very low epp value, meaning that they are parsimonious for very few nodes, which boost up
their score weights for other nodes. This was also observed for `OP411780.1`. `OV817379.1`, on the other hand, while
also have reads that map parsimoniously elsewhere and have very low epp values, some reads that do map parsimoniously to
`OV817379.1` also have very low epp values, meaning that they are parsimonious for many nodes, which downgrades their
score weights for other nodes.

I then tried pruning the reads that have obvious errors using the `--skip-singleton` option. This boosted up the rank of
`MZ722413.1` to **39** for its node score. The same thing happened for `OP411780.1`, whose node score rank was improved
to **39**.

What if I used perfect reads with no errors? I simulated 100,000 perfect reads for the same set of 10 haplotypes, and 
here are their ranks for their node scores:

```
Rank  NodeID                                                                              Score
-------------------------------------------------------------------------------------------------------
118  USA/CA-CHLA-PLM46323971/2020|MZ722413.1|2020-12-09                                  0.2122666691
337  USA/IL-CDC-QDX40817725/2022|OP411780.1|2022-08-28                                   0.0570905193
3    England/MILK-344FEB3/2022|OV817379.1|2022-01-26                                     56.9054583197
1    Denmark/DCGC-599763/2022|OY836217.1|2022-10-17                                      94.1630006842
2    USA/CO-CDPHE-41411769/2023|PP031687.1|2023-09-25                                    64.9913636546
4    USA/CO-CDPHE-2007061387/2020|OK659067.1|2020-07-06                                  34.9137689944
8    USA/WA-S20280/2022|ON660803.1|2022-05-20                                            6.9469028815
13   USA/FL-CDC-LC0637436/2022|ON608372.1|2022-05-12                                     3.6754517996
7    Germany/Molecular_surveillance_of_SARS-CoV-2_in_Germany/2021|OV342561.1|2021-09-10  8.0952134279
5    USA/CA-CDC-STM-XAH3WETJM/2022|OP911649.1|2022-11-14                                 17.6248774466
```

After filtering the reads separately for `MZ722413.1` and `OP411780.1`, each filtered sample has the true haplotype
ranked **1** for its node score. I also noticed that these two haplotypes do not have parsimonious reads that have low
epp values (lowest epp for `MZ722413.1` is 25, for `OP411780.1` is 33, `OV817379.1` as epp value of 1 for 24 collapsed
reads for comparison), meaning that all the reads for these two haplotypes are parsimonious for a non-trivial amount of
nodes, meaning their node scores can get "drowned out" by other true nodes that have more "unique" reads on the tree. 

Note that this is all done on the SARS 20K tree, which is much less balanced than the SARS 8M tree. I will run more
tests on the SARS 8M tree. I will also need to generate some consistent datasets on phoenix.

### Make a docker for panmap/panMAMA to build on phoenix

[Dockerfile](panmama/10_9_2025/Dockerfile) and [CMakeLists.txt](panmama/10_9_2025/CMakeLists.txt) can be found in 
[panmama/10_9_2025](panmama/10_9_2025).

#### Start on Docker build

To build the docker image
```
docker build -t panmap-dev .
```

To enter the interactive shell
```
docker run --rm -it \
  -v $(pwd):/panmap \
  -w /panmap \
  panmap-dev
```

To build `panmap` fresh from inside the container
```
mkdir -p /panmap/build && cd /panmap/build
cmake ..
make -j 32
```

#### Running things inside the container from the host

A simple example:
```
docker run --rm \
  -v $(pwd):/panmap \
  -w /panmap \
  --user $(id -u):$(id -g) \
  panmap-dev \
  bash -c "echo rawrrr > tiger"
```
*Using --user tag for `tiger` to be owned by user instead of root*

<br/>

To run `panmap` from the host
```
docker run --rm \
  -v $(pwd):/panmap \
  -v /scratch1/alan/goodmap/panmap/panmans:/panmap/panmans \
  -v /scratch1/alan/goodmap/panmap/test_data:/panmap/test_data \
  -w /panmap \
  --user $(id -u):$(id -g) \
  panmap-dev \
  bash -c "/panmap/build/bin/panmap /panmap/panmans/sars_20000_optimized.panman \
          /panmap/test_data/sars/rep1/sars20000_5hap-a_100000_rep1_R*.fastq \
          -m /panmap/panmans/sars_20000.pmai \
          --cpus 4"
```

To rebuild `panmap` from the host
```
docker run --rm \
  -v $(pwd):/panmap \
  -w /panmap \
  panmap-dev \
  bash -c "cd /panmap/build && make -j 32"
```
*Not using --user tag to avoid permission issues*

Nice. It also built on Phoenix. Now I can run more tests on Phoenix.

## 10/10/2025

:exclamation: Reminder: Evaluation scripts and figures for genotyping with mutation spectrum is in `panmap`'s `main`
branch (commit `f239d06`) in `dev/denotype_eval`[^1].

[^1]: This is a reminder created before I deleted `bzhan146@emerald:/private/groups/corbettlab/alan/panmap` on phoenix.

Today I will write a script on phoenix to generate benchmark data. This script should be generalizable to any tree with
customizable parameters.

### Borrowed from [10/6/2025](#1062025)

<ol>
  <li>Original and mutated haplotype
    <ol>
      <li>All original haplotypes</li>
      <li>All mutated haplotypes</li>
      <li>Mix of original and mutated haplotypes</li>
    </ol>
  </li>
  <li>Sequencing type
    <ol>
      <li>shot-gun</li>
      <li>tailed amplicon</li>
    </ol>
  </li>
  <li>Sequencing depth</li>
</ol>

I will make all combinations of (`#SNPs`, `#Haplotypes`, `%Mutated`, `SeqType`, `Depth`)

`#SNP: 0 5 10 20`

`#haplotypes: 1 5 10 50 100`

`%Mutated: 20% 50% 70% 100%` *If #SNP > 0*

`SeqType: shot-gun(0) tiled-amplicon(1)`

`#Reads: 100000 500000 1000000 2000000`

*Each combination has 5 replicates*

Use `panmama/benchmark/gencomb.py` to generate combinations of parameters.

```
python3 panmama/benchmark/gencomb.py \
        --snps 0 5 10 20 \
        --haplotypes 1 5 10 50 100 \
        --percent-mutated 0.2 0.4 0.8 1.0 \
        --seq-types shotgun amplicon \
        --num-reads 100000 500000 1000000 2000000 \
        --num-rep 5 | head | column -t
```

```
seq_type  haplotype  snp  percent_mutated  num_reads  rep
shotgun   1          0    0                100000     0
shotgun   1          0    0                100000     1
shotgun   1          0    0                100000     2
shotgun   1          0    0                100000     3
shotgun   1          0    0                100000     4
shotgun   1          0    0                500000     0
shotgun   1          0    0                500000     1
shotgun   1          0    0                500000     2
shotgun   1          0    0                500000     3
```


Use `panmama/benchmark/genreads.sh` to generate reads.

```
bash panmama/benchmark/genreads.sh \
  --seqtype shotgun \
  --numhap 5 \
  --numsnp 10 \
  --permut 0.4 \
  --numreads 1000 \
  --rep 0 \
  --cpus 4 \
  --panmap /private/groups/corbettlab/alan/panmap/ \
  --panman /private/groups/corbettlab/alan/panmap/panmans/sars_optimized.panman \
  --pmi /private/groups/corbettlab/alan/panmap/panmans/sars_optimized.panman.pmi \
  --random-seed random \
  --out-prefix panmama/benchmark/test
```

or to simulate `amplicon` reads

```
bash panmama/benchmark/genreads.sh \
  --seqtype amplicon \
  --numhap 5 \
  --numsnp 10 \
  --permut 0.4 \
  --numreads 1000 \
  --rep 0 \
  --cpus 4 \
  --panmap /private/groups/corbettlab/alan/panmap/ \
  --panman /private/groups/corbettlab/alan/panmap/panmans/sars_optimized.panman \
  --pmi /private/groups/corbettlab/alan/panmap/panmans/sars_optimized.panman.pmi \
  --random-seed random \
  --swampy /private/home/bzhan146/tools/SWAMPy/src/simulate_metagenome.py \
  --reference-primer-bed-file /private/home/bzhan146/tools/SWAMPy/primer_sets/nimagenV2.bed \
  --reference-fasta-file /private/home/bzhan146/tools/SWAMPy/ref/MN908947.3.fasta \
  --jvarkit /private/home/bzhan146/tools/jvarkit/dist/jvarkit.jar \
  --out-prefix panmama/benchmark/test
```

## 10/13/2025

I just implemented a fast node scoring function `scoreReads()` and `scoreReadsHelper()`. This will traverse the tree and
only count the number of matching k-min-mers in each read instead of pseudo-chaining them. This should make the read and
node scoring step much faster as a pre-filter step. The plan is to still do full chaining after filtering out the
probable nodes and collapsing the tree.

### Quick follow up on [10/9/2025](#1092025)

Getting sidetracked a bit but I wonder if the two true nodes, `MZ722413.1` and `OP411780.1`, have low scores  because 
reads that parismoniously  map to them also parsimoniously map to many inferred internal nodes, which drives up the epp 
and reduces the read weights that map to the true nodes.  

Not counting internal nodes when counting epp for reads, the node score ranks for the true haplotypes from sample
`test_data/sars/rep3/sars20000_10hap-a_100000_rep3_R*.fastq`:

```
2292  USA/CA-CHLA-PLM46323971/2020|MZ722413.1|2020-12-09                                  0.3971806677
3833  USA/IL-CDC-QDX40817725/2022|OP411780.1|2022-08-28                                   0.1865463517
2     England/MILK-344FEB3/2022|OV817379.1|2022-01-26                                     58.1048931881
1     Denmark/DCGC-599763/2022|OY836217.1|2022-10-17                                      85.2544216228
3     USA/CO-CDPHE-41411769/2023|PP031687.1|2023-09-25                                    54.7220291402
4     USA/CO-CDPHE-2007061387/2020|OK659067.1|2020-07-06                                  27.1159401016
8     USA/WA-S20280/2022|ON660803.1|2022-05-20                                            8.3984476356
6     USA/FL-CDC-LC0637436/2022|ON608372.1|2022-05-12                                     10.3158523194
9     Germany/Molecular_surveillance_of_SARS-CoV-2_in_Germany/2021|OV342561.1|2021-09-10  6.8616225312
5     USA/CA-CDC-STM-XAH3WETJM/2022|OP911649.1|2022-11-14                                 19.1232183505
```

It looks a little bit better but still not very good. Never mind then...

### Use fast-mode as default

Since the fast-mode, which doesn't pseudo-chain the k-min-mer matches, outputs very similar results as the normal mode,
Russ and I decided to have the fast mode as the default option for faster speed. I will reserve a whole node on phoenix
to do some benchmark on the runtime.

```
sbatch /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/fast-mode.sh --array=0-3
```

## 10/14/2025

I did a little bit of optimization and now scoring 1.5 million amplicon short reads on the SARS 8M took ~8 minutes,
using 32 threads.

### Search space reduction methods options

As discussed on [10/9/2025](#some-nodes-were-just-not-selected-by-wepp-scores-for-some-reason), some nodes are not
selected using WEPP scoring methods because they are very similar to other nodes on the tree, meaning that all of the 
reads that map pasimoniously to them have very high EPP values, thus making them have less scores and more susceptible
to sequencing errors.

Russ and I discussed some other options:

1. Find local optimas of sum read scores
2. Find local optimas of WEPP node scores
3. Regularize EPP by giving it a "span value", which measures the genetic distance between the parsimonious nodes...
e.g. if all the parsimonious nodes of a read are neighbors, give the read weight a positive correction to increase the 
node scores of those nodes.
4. Combine option 2 and 3

I will be trying out these options in the next couple of dates

### Efficient way to compute EPPs

I can actually compute the EPP during scoring with minimum overhead. During scoring, if a read has a new max score, I 
can set a `read.lastParsimoniousNode` to the current node. Downstream, when the read's score is changed from max to
a low score, I can calculate the number of parsimoniously mapped nodes by
`epp += curNode.dfsIndex - read.lastParsimoniousNode.dfsIndex`. If a new max score is reached, I can simply reset the
`epp`.

As implemented:
```c++
auto& lastParsimoniousNode = curRead.lastParsimoniousNode;
if (score > maxScore) {
  // reset epp and assign current node as lastParsimoniousNode
  curRead.epp = 0;
  lastParsimoniousNode = node;
} else if (score < maxScore) {
  // collect epp and reset lastParsimoniousNode
  if (lastParsimoniousNode != nullptr) {
    curRead.epp += node->dfsIndex - lastParsimoniousNode->dfsIndex;
    lastParsimoniousNode = nullptr;
  }
} else {
  if (lastParsimoniousNode == nullptr) {
    lastParsimoniousNode = node;
  }
}
```

This method would require me to handle identical nodes on the tree more carefully, which is a WIP right now.

## 10/15/2025

Today I finished collapsing the tree.

WIP: Efficient counting of read EPP during scoring. It's actually a bit trickier than described on
[10/14/2025](#efficient-way-to-compute-epps).

## 10/16/2025

After roughly a day, I finally correctly implememnted an efficient counting of read EPPs during scoring. Now the program
is able to score the SARS 8M tree and simultaneously count the read EPPs in ~6.5 minutes using 32 threads.

For bookkeeping, this is implemented in commit `ccc3bcc73b31ba686b4a43d8fedef5180d1f1582`.

Now I will experiment with various methods to reduce search space, as described on
[10/14/2025](#search-space-reduction-methods-options).

Just implemented calculation of the sum of raw read matches for each node on the tree (commit
`eff3a812d9d76f69f0f04e78fb98221b4998d6c9`), which took ~180s for 1.5 million amplicon short reads on the SARS 8M tree.
Will work on the WEPP score tomorrow.

## 10/17/2025

Just implemented the calculate of WEPP node scores in commit `2a327e40c1a6787b56971a148a989e7ac00a9d16`. Dynamic
calculating of sum raw read matches and WEPP scores are confirmed to be correct after checking with brute foce on the
SARS 20K tree.

## 10/18/1025


### Generate test data on phoenix

For a given combination of num_snps, num_haps, percent_mutated, and replicate_number, generate the same set
of haplotypes and abundances for amplicon and shotgun reads with various depth.

```
sbatch /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/gendata.sh \
  /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/data_sars_20K \
  /private/groups/corbettlab/alan/panmap/panmans/sars_optimized.panman \
  /private/groups/corbettlab/alan/panmap/panmans/sars_20000.pmai
```

and 

```
sbatch /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/gendata.sh \
  /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/data_sars_8M \
  /private/groups/corbettlab/alan/panmap/panmans/sars_8M.panman \
  /private/groups/corbettlab/alan/panmap/panmans/sars_8M.pmai
```

I compressed the outputs into `data_sars_20K.tar.gz` and  `data_sars_8M.tar.gz` and removed the original output dirs for
now. I also scp'ed the `*tar.gz` to `alan@silverbullet:/scratch1/alan/lab_notebook/panmama/benchmark`.

## 10/19/2025

### Oops, I generated the wrong data...

I forgot to change the panman path when I copied over the bash command to run `gendata.sh` and ended up generating 
samples of SARS 20K haplotypes in the `data_sars_8M/` dir. Anyway, regenereting data rn.

### Trying different scoring schemes

I wrote several node scoring schemes for selecting probable nodes.

```
sbatch /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/run_panmap_score_nodes_sars_20K.sh \
  /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/node_scores_out \
  /private/groups/corbettlab/alan/panmap/ \
  /private/groups/corbettlab/alan/panmap/panmans/sars_optimized.panman \
  /private/groups/corbettlab/alan/panmap/panmans/sars_20000.pmai \
  /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/data_sars_20K
```

```
sbatch /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/run_panmap_score_nodes_large_trees.sh \
  /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/node_scores_out \
  /private/groups/corbettlab/alan/panmap/ \
  /private/groups/corbettlab/alan/panmap/panmans/sars_8M.panman \
  /private/groups/corbettlab/alan/panmap/panmans/sars_8M.pmai \
  /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/data_sars_8M
```

## 10/20/2025 & 10/21/2025

I just implemented the full WEPP scoring scheme, which includes updating the rest of the node scores every time a top
scoring haplotype was selected. This was not implemented in the previous commit described on [10/17/2025](#10172025).

The nodes scoring step and especially the node score updating step were not as fast as I'd like. So I spent the last two
days making them fast enough to my satisfaction.

Hmmm.. It seems like the penalty of `1/epp^2`, which is used in WEPP, might be too aggressive. I tried `1/epp` and it
actually doesn't look too bad. I will send this to phoenix and do a couple rounds of testing.

## 10/22/2025

### processing nodeScores output

Use `grep_true.sh` to view the node score ranking of an output `nodeScores.tsv` file.

```
cd /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark && \
bash grep_true.sh node_scores_out/sars_optimized_amplicon_10_0_0_1500000_1.nodeScores.tsv data_sars_20K/
```

Then run it on all the outputs

```
cd /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark && \
bash grep_true_all.sh node_scores_out/ > node_scores_ranks.txt
```

### Genereating data for mixed HIV and RSV samples

Only generating shotgun sequencing for now.

```
sbatch /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/gendata.sh \
  /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/data_hiv \
  /private/groups/corbettlab/alan/panmap/panmans/hiv_optimized.panman \
  /private/groups/corbettlab/alan/panmap/panmans/hiv_optimized.pmai
```

```
sbatch /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/gendata.sh \
  /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/data_rsv \
  /private/groups/corbettlab/alan/panmap/panmans/rsv_optimized.panman \
  /private/groups/corbettlab/alan/panmap/panmans/rsv_optimized.pmai
```


and testing node scoring schemes on pheonix

```
sbatch /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/run_panmap_score_nodes_hiv_20K.sh \
  /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/node_scores_out_hiv \
  /private/groups/corbettlab/alan/panmap/ \
  /private/groups/corbettlab/alan/panmap/panmans/hiv_optimized.panman \
  /private/groups/corbettlab/alan/panmap/panmans/hiv_optimized.pmai \
  /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/data_hiv
```

### Identical nodes found that far apart on the tree while their own neighbors are very different

On the SARS 20K tree I found some nodes that have identical sequences but are non-neighbors, and their actual respective
neighbors are quite different from themselves. These are not collapsed during the initial tree collapse because they are
not neighbors.

One (hopefully easy and not computationally intensive) fix would be to first sort the nodes by their raw read seed
matches and raw maximum-placement read seed matches. Groups of nodes with identical scores are potential idenical nodes.
I can then find their LCA and apply read score deltas to each of them to see if they have identical nodes. I hope this
wouldn't be too computationally intensive because I expect identical score node groups would be small.

I will do pairwise comparison of all collapsed leaf and internal nodes
(including leaf node vs leaf node, leaf node vs internal node, BUT NOT internal node vs internal node)... Using
collapsed leaf nodes to avoid comparing neighboring nodes, which are not unexpected if they were identical.

Actually, there would be too many comparisons... Comparing all leaf node pairs would have `(N * (N - 1)) / 2`
comparisons. That's 199,990,000 comparisons!

I have a way to reduce the search space I will use a random node score files and first find potential identical
groups that have identical raw read seed matches and identical raw maximum-placement read seed matches.

```
cd /private/groups/corbettlab/alan/lab_notebook/panmama/pairwise_comparison &&
python3 search_potential_identical.py panmap.nodeScores.tsv  > potential_identical_pairs_sars20k.tsv
```

```
sbatch pairwise_comparison.sh \
  potential_identical_pairs_sars20k.tsv \
  /private/groups/corbettlab/alan/lab_notebook/panmama/pairwise_comparison/sars_fasta \
  pairwise_comparison_out_sars20k
```

Pairs of potential identical nodes are evenly distributed for each task array, and distance stats are printed to each
task array's individual output. After the job is complete, I will `cat` all the outputs together.

## 10/23/2025

All the pairwise comparisons are complete and I concatenated each task array's output to a final all_comparisons.tsv

```
cd /private/groups/corbettlab/alan/lab_notebook/panmama/pairwise_comparison &&
python3 get_identical_groups.py pairwise_comparison_out_sars20k/all_comparisons.tsv > sars_20k_identical_group_indices.tsv
```

It outputs a tsv file where first column is the sequence ids and the second column is the identical group ids. For now,
I skipped over groups containing only internal nodes. Nodes are considered identical if their snps, snps_ambiguous, gaps,
and gaps_edge_corrected are all 0 (this can be modified to be more lenient, such as allowing snps_ambiguous or gaps as
long as gaps_edge_corrected is 0).

I did find many groups of identical leaf nodes on the SARS 20K tree. I also did the same for RSV 4K and HIV 20K. RSV 4K
has much less while HIV 20K actually has a lot more.

## 10/28/2025

I will simulate SARS data that are more realistic. Before I had randomly selected 100 nodes from the tree.
Realistically, wastewater data should contain many similar haplotypes. So I'm going to selected random clusters of nodes
to simulate reads from.

Using `panmap` commmit `c0aeb3a631aa0ab9a1dfb2d149d6d6362a316c58` to selected the random clusters of nodes.

For each cluster, I randomly select a node then use Dijkstra's algorithm to selected its closet N neighbors. Not sure
if I should use weighted branch (and how to weights the branch length) or not but using branch length weighted by the
number of k-min-mer updates for simlicity for now.

```
sbatch /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/gendata_clustered_nodes.sh \
  /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/data_sars_20K_clustered \
  /private/groups/corbettlab/alan/panmap/panmans/sars_20000_twilight_dipper.panman \
  /private/groups/corbettlab/alan/panmap/panmans/sars_20000_twilight_dipper.pmai
```

```
sbatch /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/gendata_clustered_nodes.sh \
  /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/data_sars_8M_clustered \
  /private/groups/corbettlab/alan/panmap/panmans/sars_8M.panman \
  /private/groups/corbettlab/alan/panmap/panmans/sars_8M.pmai
```

What to do tomorrow:
  - Check on read simulation. Resubmit jobs if needed
  - Try out how the new read samples looks for node selection
  - Keep thinking of new ways to select probable nodes


## 10/30/2025

Results look much better using reads simulated from clusters of nodes. On the sars_20K tree, I was able to select 95%
to 100% of the true nodes from which the reads are simulated from, for both shotgun and amplicon reads. However, it
doesn't look very good on the 8M tree, because for a random sample I looked at, only 9 out of 50 of the true nodes were
selected. Nevertheless, I will do a formal evaluation with replicates on phoenix.

I also deviced another scoring method to use seed weights instead of read weights, where a seed weight is proportional 
to its rarity on the the tree, similar to the rarity of a read's parsimony placement. By the first look of it, it's not
bad but also not as good as the read weighted scoring method on the 20K tree. I will also include this as part of the 
evaluation described above.

```
sbatch /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/run_panmap_score_nodes_sars_20K.sh \
  /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/node_scores_out_sars20K \
  /private/groups/corbettlab/alan/panmap/ \
  /private/groups/corbettlab/alan/panmap/panmans/sars_20000_twilight_dipper.panman \
  /private/groups/corbettlab/alan/panmap/panmans/sars_20000_twilight_dipper.pmai \
  /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/data_sars_20K_clustered
```

```
sbatch /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/run_panmap_score_nodes_sars_8M.sh \
  /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/node_scores_out_sars8M \
  /private/groups/corbettlab/alan/panmap/ \
  /private/groups/corbettlab/alan/panmap/panmans/sars_8M.panman \
  /private/groups/corbettlab/alan/panmap/panmans/sars_8M.pmai \
  /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/data_sars_8M_clustered
```

## 10/31/2025

### Seed-weighted node scores update is too slow on the 8M tree

It seems like the node seed scores update step was taking too long on the 8M tree. All of the jobs for the seed weighted
node scoring scheme ran out of time. I will make it faster. 


### Seed-weighted node scores performed better than I thought

I was wrong about seed-weighted node scoring on [10/30/2025](#10302025). After doing a formal evaluation, it actually
performed better than read-weighted node scoring in more cases than not.

```
python3 eval_node_scores.py node_scores_out_sars20K/ > node_scores_out_sars20K_all.tsv
(head -n1 node_scores_out_sars20K_all.tsv && tail -n+2 node_scores_out_sars20K_all.tsv | sort -t '_'  -k2,2n -k6,6n -k5,5n -k 1,1) > tmp
mv tmp node_scores_out_sars20K_all.tsv
```

Out of 60 simulated cases:
- 15 cases of 100,000 shotgun reads
  - equal: 7
  - seed-weight: 3
  - read-weight: 5
- 15 cases of 1,500,000 shotgun reads
  - equal: 6
  - seed-weight: 7
  - read-weight: 2
- 15 cases of 100,000 amplicon reads
  - equal: 2
  - seed-weight: 8
  - read-weight: 5
- 15 cases of 1,500,000 amplicon reads
  - equal: 8
  - seed-weight: 5
  - read-weight: 2



### Simulating perfect reads

I think I will also simulate some **PERFECT** reads to see how PCR and sequencing errors affect the scoring and node
selection. If I get much better node selection on error-free reads, I will spend a little more time on identifying and
removing likely errors.

```
sbatch /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/gendata_clustered_nodes_perfect.sh \
  /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/data_sars_20K_clustered \
  /private/groups/corbettlab/alan/panmap/panmans/sars_20000_twilight_dipper.panman \
  /private/groups/corbettlab/alan/panmap/panmans/sars_20000_twilight_dipper.pmai
```

```
sbatch /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/gendata_clustered_nodes_perfect.sh \
  /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/data_sars_8M_clustered \
  /private/groups/corbettlab/alan/panmap/panmans/sars_8M.panman \
  /private/groups/corbettlab/alan/panmap/panmans/sars_8M.pmai
```

### Regarding my trouble with the SARS 8M

It seems like the SARS 8M tree has some fundamental limitaitons in terms of its phylogenetic correctness. This arises
primarily due to edge sequences. Two very similar genomes (a few SNPs apart) can be positioned very far apart on the 
tree when one genome includes edges in its assembly while the other does not. Edge sequences can alse create more seed
hits for genomes that include them. Relying on seed hits, per our method, can also have reduced specificity.

In terms of what to do next. I think I will do some **final testing** on the 8M tree to see if perfect reads without any
sequencing or PCR errors can improve the node selection performance. If it does, then I know it's primarily due to 
sequencing or PCR errors that create false parsimonious matches, which jack up the node scores of other false nodes to 
be selected over the true nodes.

### What to focus on instead

After talking with Russ. I will instead start focusing on more divergent genomes, such as HIV and bacterial genomes
(TB? etc.). We may also need to build our own panMANs of other divergent pathogen species.

To select out relevent reads from shotgun sequencing, use Kraken or discard reads without significant hits during read
scoring in panMAMA.



## 11/3/2025

### To-do

**Last bit of attemp/investigation on 8M tree**
- [x] Finish the new faster seed-weighted node score update.
- [x] Run node selection on error-free reads.
- [ ] Evaluate the results

**More divergent tree**
- [ ] Find SARS/HIV/etc. isolate samples and create artificial mixtures to test for demixing and genome deconvolution
- [ ] Find real mixed samples for divergent genomes, such as HIV and TB
- [ ] Build our own panMANs if necessary

### Re-running node scores with fast, debugged seed-weighted node scres

I just implemented a much faster way to score seed-weighted scores. I also found a bug in the original code so I will
rerun the node scores for all the simualted samples.

```
sbatch /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/run_panmap_score_nodes_sars_20K.sh \
  /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/node_scores_out_sars20K \
  /private/groups/corbettlab/alan/panmap/ \
  /private/groups/corbettlab/alan/panmap/panmans/sars_20000_twilight_dipper.panman \
  /private/groups/corbettlab/alan/panmap/panmans/sars_20000_twilight_dipper.pmai \
  /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/data_sars_20K_clustered
```
job id: `21214699` `21217044`

```
sbatch /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/run_panmap_score_nodes_sars_8M.sh \
  /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/node_scores_out_sars8M \
  /private/groups/corbettlab/alan/panmap/ \
  /private/groups/corbettlab/alan/panmap/panmans/sars_8M.panman \
  /private/groups/corbettlab/alan/panmap/panmans/sars_8M.pmai \
  /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/data_sars_8M_clustered
```
job id: `21214720` `21217989`

Also running with perfect samples


```
sbatch /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/run_panmap_score_nodes_sars_20K_perfect.sh \
  /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/node_scores_out_sars20K \
  /private/groups/corbettlab/alan/panmap/ \
  /private/groups/corbettlab/alan/panmap/panmans/sars_20000_twilight_dipper.panman \
  /private/groups/corbettlab/alan/panmap/panmans/sars_20000_twilight_dipper.pmai \
  /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/data_sars_20K_clustered
```
job id: `21215120` `21217046`

```
sbatch /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/run_panmap_score_nodes_sars_8M_perfect.sh \
  /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/node_scores_out_sars8M \
  /private/groups/corbettlab/alan/panmap/ \
  /private/groups/corbettlab/alan/panmap/panmans/sars_8M.panman \
  /private/groups/corbettlab/alan/panmap/panmans/sars_8M.pmai \
  /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/data_sars_8M_clustered
```
job id: `21215142`

## 11/4/2025

### To-do

**Last bit of attemp/investigation on 8M tree**
- [x] Run node selection on error-free reads.
- [ ] Evaluate the results

**More divergent tree**
- [ ] Find SARS/HIV/etc. isolate samples and create artificial mixtures to test for demixing and genome deconvolution
- [ ] Find real mixed samples for divergent genomes, such as HIV and TB
- [ ] Build our own panMANs if necessary

### Rerunning error-free reads

It doesn't seem like the amplicon reads were error-free. Gonna take a look rn.

Yup, just adding `--errfree` in `art_runner.py` is not enough. I needed to convert the error-free SAM ouput from `art`
to fastqs.

Rerunning scripts to generate real error free amplicon reads.

```
sbatch \
  --mem=30000 \
  --cpus-per-task=8 \
  /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/gendata_clustered_nodes_perfect.sh \
  /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/data_sars_20K_clustered \
  /private/groups/corbettlab/alan/panmap/panmans/sars_20000_twilight_dipper.panman \
  /private/groups/corbettlab/alan/panmap/panmans/sars_20000_twilight_dipper.pmai
```
job id:  `21225800`

```
sbatch \
  --mem=60000 \
  --cpus-per-task=8 \
  /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/gendata_clustered_nodes_perfect.sh \
  /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/data_sars_8M_clustered \
  /private/groups/corbettlab/alan/panmap/panmans/sars_8M.panman \
  /private/groups/corbettlab/alan/panmap/panmans/sars_8M.pmai
```
job id: `21226337`


And rerunning panMAMA with the perfect reads

```
sbatch /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/run_panmap_score_nodes_sars_20K_perfect.sh \
  /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/node_scores_out_sars20K \
  /private/groups/corbettlab/alan/panmap/ \
  /private/groups/corbettlab/alan/panmap/panmans/sars_20000_twilight_dipper.panman \
  /private/groups/corbettlab/alan/panmap/panmans/sars_20000_twilight_dipper.pmai \
  /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/data_sars_20K_clustered
```
job id: `21226880`

```
sbatch /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/run_panmap_score_nodes_sars_8M_perfect.sh \
  /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/node_scores_out_sars8M \
  /private/groups/corbettlab/alan/panmap/ \
  /private/groups/corbettlab/alan/panmap/panmans/sars_8M.panman \
  /private/groups/corbettlab/alan/panmap/panmans/sars_8M.pmai \
  /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/data_sars_8M_clustered
```
job id: `21226921` `21228336`

### Searching real samples

#### SARS

I think I will use the 10 SARS samples Nico found.

#### HIV

Here are a list of potential samples (to be updated)

| Project/Sample ID | Strategy | Layout | #Samples | Data/Paper link |
| ----------------- | -------- | ------ | -------- | ---------- |
| PRJNA1118440 | HIV-1 Probe-capture WGS | paired | 82 | [Coldbeck-Shackley et al.](https://pubmed.ncbi.nlm.nih.gov/38924832/) |
| PRJNA644953 | HIV-1 partial pol gene (Protease, Reverse Transcriptase, and Integrase) RT-PCR | single | 99 | [data](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA644953) | 

There's not many data available for HIV-2

### Downloading fasta sequences from ncbi

Use commands below to retrieve fasta sequences from ncbi

**Download a single sequence**

`efetch -db nuccore -id ON284258.1 -format fasta > ON284258.1.fasta`

**Download multiple sequences from a file (one accession per line)**

`cat accessions.txt | efetch -db nuccore -format fasta > sequences.fasta`

**Or download multiple specific accessions**

`efetch -db nuccore -id ON284258.1,OL672836.1,MW580573.1 -format fasta > sequences.fasta`

## 11/6/2025 - 11/7/2025

### I spent the majority of my time yesterday (11/5/2025) preparing for a BME seminar talk.

### Amplicon sample error detection using primer-stack size

Perfect amplicon reads (without sequencing/PCR errors) actually did quite well during node selection step. So I do think
I will spend a bit more time handling amplicon errors.

I used Claude to write a watered-down python version if `ivar trim` to output a tsv file, in which the first column is
the read name and second column is the name of the primer that it's assigned to.

During query initialization step, I will first process reads by groups of their primer template and mask probable errors
using the depth of the primer stack as reference.

```
sbatch \
  --mem=30000 \
  --cpus-per-task=8 \
  /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/gendata_clustered_nodes_amplicon_stack.sh \
  /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/data_sars_20K_clustered \
  /private/groups/corbettlab/alan/panmap/panmans/sars_20000_twilight_dipper.panman \
  /private/groups/corbettlab/alan/panmap/panmans/sars_20000_twilight_dipper.pmai
```
job id:  `21271306`

```
sbatch \
  --mem=60000 \
  --cpus-per-task=8 \
  /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/gendata_clustered_nodes_amplicon_stack.sh \
  /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/data_sars_8M_clustered \
  /private/groups/corbettlab/alan/panmap/panmans/sars_8M.panman \
  /private/groups/corbettlab/alan/panmap/panmans/sars_8M.pmai
```
job id: `21326790`

### Running panMAMA with amplicon error removal

For fair comparison, I will also rerun panmap on the same newly simulated datasets with and without error removal:

```
sbatch /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/run_panmap_score_nodes_sars_20K.sh \
  /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/node_scores_out_sars20K \
  /private/groups/corbettlab/alan/panmap/ \
  /private/groups/corbettlab/alan/panmap/panmans/sars_20000_twilight_dipper.panman \
  /private/groups/corbettlab/alan/panmap/panmans/sars_20000_twilight_dipper.pmai \
  /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/data_sars_20K_clustered
```
job id: `21327409`

```
sbatch /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/run_panmap_score_nodes_sars_8M.sh \
  /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/node_scores_out_sars8M \
  /private/groups/corbettlab/alan/panmap/ \
  /private/groups/corbettlab/alan/panmap/panmans/sars_8M.panman \
  /private/groups/corbettlab/alan/panmap/panmans/sars_8M.pmai \
  /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/data_sars_8M_clustered
```
job id: `21327420`

```
sbatch /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/run_panmap_score_nodes_sars_20K_error_detect.sh \
  /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/node_scores_out_sars20K \
  /private/groups/corbettlab/alan/panmap/ \
  /private/groups/corbettlab/alan/panmap/panmans/sars_20000_twilight_dipper.panman \
  /private/groups/corbettlab/alan/panmap/panmans/sars_20000_twilight_dipper.pmai \
  /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/data_sars_20K_clustered
```
job id: `21368167`

```
sbatch /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/run_panmap_score_nodes_sars_8M_error_detect.sh \
  /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/node_scores_out_sars8M \
  /private/groups/corbettlab/alan/panmap/ \
  /private/groups/corbettlab/alan/panmap/panmans/sars_8M.panman \
  /private/groups/corbettlab/alan/panmap/panmans/sars_8M.pmai \
  /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/data_sars_8M_clustered
```
job id: `21368178` `21368977`

## 11/12/2025

### Read-seed-weighted node scoring function

```
sbatch /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/run_panmap_score_nodes_sars_20K.sh \
  /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/node_scores_out_sars20K \
  /private/groups/corbettlab/alan/panmap/ \
  /private/groups/corbettlab/alan/panmap/panmans/sars_20000_twilight_dipper.panman \
  /private/groups/corbettlab/alan/panmap/panmans/sars_20000_twilight_dipper.pmai \
  /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/data_sars_20K_clustered
```
job id: `21476247`

```
sbatch /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/run_panmap_score_nodes_sars_8M.sh \
  /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/node_scores_out_sars8M \
  /private/groups/corbettlab/alan/panmap/ \
  /private/groups/corbettlab/alan/panmap/panmans/sars_8M.panman \
  /private/groups/corbettlab/alan/panmap/panmans/sars_8M.pmai \
  /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/data_sars_8M_clustered
```
job id: `21476259`

```
sbatch /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/run_panmap_score_nodes_sars_20K_error_detect.sh \
  /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/node_scores_out_sars20K \
  /private/groups/corbettlab/alan/panmap/ \
  /private/groups/corbettlab/alan/panmap/panmans/sars_20000_twilight_dipper.panman \
  /private/groups/corbettlab/alan/panmap/panmans/sars_20000_twilight_dipper.pmai \
  /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/data_sars_20K_clustered
```
job id: `21476275`

```
sbatch /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/run_panmap_score_nodes_sars_8M_error_detect.sh \
  /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/node_scores_out_sars8M \
  /private/groups/corbettlab/alan/panmap/ \
  /private/groups/corbettlab/alan/panmap/panmans/sars_8M.panman \
  /private/groups/corbettlab/alan/panmap/panmans/sars_8M.pmai \
  /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/data_sars_8M_clustered
```
job id: `21476287` `21492432` `21497204`

```
sbatch /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/run_panmap_score_nodes_sars_20K_perfect.sh \
  /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/node_scores_out_sars20K \
  /private/groups/corbettlab/alan/panmap/ \
  /private/groups/corbettlab/alan/panmap/panmans/sars_20000_twilight_dipper.panman \
  /private/groups/corbettlab/alan/panmap/panmans/sars_20000_twilight_dipper.pmai \
  /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/data_sars_20K_clustered
```
job id: `21476302`

```
sbatch /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/run_panmap_score_nodes_sars_8M_perfect.sh \
  /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/node_scores_out_sars8M \
  /private/groups/corbettlab/alan/panmap/ \
  /private/groups/corbettlab/alan/panmap/panmans/sars_8M.panman \
  /private/groups/corbettlab/alan/panmap/panmans/sars_8M.pmai \
  /private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/data_sars_8M_clustered
```
job id: `21476314`

## 11/13/2025

So read-seed-weighted scoring method didn't work out well. While it performed better than read-weighted method on the
SARS 20K tree, it's comparably worse on the SARS 8M tree. During my discussion with Russ and Alex, weighting the reads
additionally with phylogenetic information (similar to autolin... information theory... how much phylogenetic
information does a read's parsimony tells me). I did think about something like this but decided not to pursue because
the trees are not quite good enough.

Wait... Now thinking back on it, since SARS 8M tree was built using UShER MAT, the tree topology might be good enough to
do this.

### Quick plan on what to do next

- [ ] Investigate the phylogeny of nodes that reads are parsimonious to on the 8M tree
  - [ ] If it looks good: pursue phylogenetic information approach mentioned above
- [ ] Restore EM function on panMAMA
- [ ] Test on more divergent genomes!
  - [ ] Simulate RSV amplicon reads and compare panMAMA to WEPP
    - [ ] Figure out from `art_illumina` page how to simulate amplicon reads
    - [ ] Write a script to simulate amplicon reads for RSV and other more divergent genomes
  - [x] Modify WEPP to also use shotgun reads
  - [x] Test if WEPP can run with HIV shotgun reads and compare to panMAMA

## 11/14/2025

### Modify WEPP to also use shotgun reads

After pulling the WEPP docker...

Enter docker container: 

`docker run -it pranavgangwar/wepp:latest` 

Use modified scripts for shotgun sequencing:

```
docker run -it \
  -v /home/alan/tools/WEPP/src/WEPP/qc_preprocess.py:/WEPP/src/WEPP/qc_preprocess.py \
  -v /home/alan/tools/WEPP/config/config.yaml:/WEPP/config/config.yaml \
  -v /home/alan/tools/WEPP/workflow/rules/qc.smk:/WEPP/workflow/rules/qc.smk \
  -v /home/alan/tools/WEPP/shotgun_reads:/data \
  -v /scratch1/alan/hiv_usher:/hiv_usher \
  pranavgangwar/wepp:latest
```

Run with `SHOTGUN=True` tag inside container:

```
snakemake --config DIR=shotgun_reads FILE_PREFIX=test_run PRIMER_BED=RSVA_all_primers_best_hits.bed TREE=hiv_K03455.pb REF=K03455.fasta  CLADE_IDX=-1 SHOTGUN=True --cores 32 --use-conda
```

### Evaluation WEPP and MAMA's performance on HIV shotgun reads

```
python3 calculate_distance_score.py \
  --tree-in test/hiv_K03455.pb \
  --true-in test/true_abundance.tsv \
  --estm-in test/estm_abundance_wepp.tsv \
  --tool WEPP \
  --fasta-dir /scratch1/alan/panmania/backup/panman/panmans/info/hiv_original/ \
  --prefix test
```

I'm actually quite unsure about how to compare panMAMA and WEPP... It's easy if an estimated haplotype is a sampled
haplotype because I can either do `mafft` or use the usher MAT to calculate the distance between the estimated haplotype
and its closest true haplotype. The problem is when an estimated haplotype is an internal node because panMAN and MAT
do not have the same internal nodes. panMAN also has the full genomic sequences of the internal nodes while MAT only has
the SNPs with respect to the reference. So I'm not sure how to deal with gaps....

I also noticed that `bte.get_haplotype()`, which "returns the complete set of mutations (haplotype) the indicated node
has with respect to the reference," also has mutations at positions where a node has a gap with respect to the
reference. I think that this might be due to how `bte.get_haplotype()` gets the mutations. It might accumulate
all the branch mutations from the root the target node, and some of the branch mutations on the path are infered
mutations on the gap coordinates of the target node. Additionally, when I place the same sample onto the tree, after 
changing the sample ID, the new sample's mutations from `bte.get_haplotype()` don't have any mutations at the gapped
positions.... Very weird


## 11/15/2025 - 11/16/2025

### Lab meeting presentation planner

1. Remind everybody the background and purpose of panMAMA
2. Same as the seminar talk
3. Transition to the problems I've solved and am still solving
    - New index building function: show a quick overview on method and a graph to show the time it takes to index each
      genome
    - Improve scoring time
      - de-duplicating reads
      - Multi-threading
      - collapsing nodes
      - removing seed mutations within each node
      - fast mode
      - remove most of the unnecessary block mutations
    - Memory problem
      - Even with score-annotated index, memory still seems to be a problem for sars 8M tree
      - Use de bruijn graph to group reads with similar seeds together and store neighboring reads deltas
    - Select probable nodes
      - Why overlap-coefficients no long work
      - Alternative methods
        - WEPP method
        - Can we make it better?
          - seed-weighted method -> nope
          - read-seed-weighted method -> better than seed-weighted.. but still nope
          - what to do next?
    - Various panMAN-related problems
      - panGraph sucks
        - explain what panGraph is and how panMAN uses panGraph
          - F'd up block mutations
            - How I some-what fixed them
          - Identical samples are not close to each other
            - Ways to fix them
          - Should-be-sister nodes are not close because their genome edges are also included in the construction

## 11/20/2025

Lab meeting yesterday was a success.

### Shift focus on HIV

I just wrote a script on silverbulelt that makes it easier to run the WEPP

I will run it on all the hiv samples I have on phoenix and compare the results between WEPP and panMAMA

## 11/21/2025

I ran all of the simualed HIV samples I have on WEPP and compared it to MAMA. Weighted Haplotype/Peak distances were
calculated using [calculate_distance_score.py](panmama/WEPP/scripts/calculate_distance_score.py).

Gathered all the measurements:

```
cd /private/groups/corbettlab/alan/lab_notebook/panmama/WEPP/WEPP_results/results
(echo -e "hap\trep\tmama_whd\tmama_wpd\twepp_whd\twepp_wpd" && for hap in 1 2 3 5 10; do for rep in 1 2 3 4 5 6 7 8 9 10; do wepp_distances=$(tail -n+2 "hiv20000_${hap}hap-a_60000_rep${rep}/hiv_${hap}hap_rep${rep}.wepp.distance_sum.tsv" | cut -f 2,3); mama_distances=$(tail -n+2 "hiv20000_${hap}hap-a_60000_rep${rep}/hiv_${hap}hap_rep${rep}.mama.distance_sum.tsv" | cut -f 2,3);  echo -e -n "${hap}\t${rep}\t${mama_distances}\t"; echo -e "$wepp_distances"; done; done;) > ../../distances.tsv 
```

|hap|rep|mama_whd|mama_wpd|wepp_whd|wepp_wpd|
|---|---|--------|--------|--------|--------|
|1|1|0.000|0.000|0.000|0.000|
|1|2|0.000|0.000|0.000|0.000|
|1|3|0.000|0.000|0.000|0.000|
|1|4|0.000|0.000|0.000|0.000|
|1|5|0.000|0.000|0.000|0.000|
|1|6|0.000|0.000|0.000|0.000|
|1|7|0.000|0.000|154.000|154.000|
|1|8|0.000|0.000|0.000|0.000|
|1|9|0.000|0.000|0.000|0.000|
|1|10|0.000|0.000|0.000|0.000|
|2|1|0.000|0.274|0.000|0.290|
|2|2|0.000|0.000|0.000|0.016|
|2|3|0.000|0.000|0.000|0.855|
|2|4|0.000|0.000|0.000|10.817|
|2|5|0.000|0.000|0.000|0.000|
|2|6|0.000|0.000|11.400|13.376|
|2|7|0.000|0.000|0.000|5.442|
|2|8|0.000|0.157|0.000|0.740|
|2|9|0.000|0.000|0.000|0.000|
|2|10|0.000|0.000|0.000|0.000|
|3|1|0.000|0.000|0.000|9.486|
|3|2|0.000|0.000|0.300|7.708|
|3|3|0.000|0.000|11.100|13.487|
|3|4|0.000|0.000|49.600|24.703|
|3|5|0.000|0.109|5.700|9.067|
|3|6|0.000|0.000|3.400|13.945|
|3|7|0.000|0.000|0.600|38.168|
|3|8|0.000|0.053|0.200|21.974|
|3|9|0.000|0.000|0.000|26.607|
|3|10|0.000|0.147|0.000|0.068|
|5|1|0.000|0.000|0.000|37.111|
|5|2|0.250|0.250|0.050|14.396|
|5|3|0.000|0.000|0.000|21.012|
|5|4|0.000|0.000|15.400|75.145|
|5|5|0.000|0.000|4.300|45.224|
|5|6|0.000|0.000|0.000|78.194|
|5|7|0.000|0.000|0.000|32.938|
|5|8|0.000|0.000|0.500|21.572|
|5|9|0.000|0.000|0.000|14.293|
|5|10|0.000|0.000|2.000|46.489|
|10|1|1.800|1.754|3.150|47.646|
|10|2|0.100|0.122|1.850|67.739|
|10|3|0.000|0.556|3.700|68.996|
|10|4|0.100|0.151|5.950|60.707|
|10|5|0.000|0.000|1.200|109.556|
|10|6|0.000|0.000|3.750|71.918|
|10|7|0.000|0.000|0.000|81.815|
|10|8|0.000|0.000|0.300|102.034|
|10|9|0.050|0.133|300.900|38.395|
|10|10|0.000|0.389|0.050|37.351|

Results show that panMAMA outperforms WEPP for HIV samples, likely because HIV genomes are too divergent to be demixed 
using SNP information alone. PanMAMA shows better performance in both the WHD and WPD metrics, which measure sensitivity 
and specificity, respectively. The improvement is particularly significant for specificity.

## 11/24/2025

Yatish and Sumit sent us the verbetrate mito panMAT. I generated a metadata file containing samples' taxonomic. Out of
the total 15,655 samples, there are 7,893 unique species, 854 families, and 157 orders.

I visually inspected the panMAT using Taxonium and found that, while most samples cluster as expected, there are a
non-trivial number of outliers, and several expected clades are fragmented into multiple separate clusters.

1. The highlighted nodes represent order Primates. There are two separate small primate clusters that are positioned 
distinctly from the major primate cluster. Additionally, non-human great apes are also quite far away from the human 
samples.

![primates](panmama/mito_panmat/image-2.png)

(For validation, I aligned a human mitochondrial sample (J01415.2) to a pileated gibbon sample (AB504749.1) from the 
nearest clade to the human sample and to a chimpanzee sample (D38113.1) using MAFFT. The genetic distance between human 
and chimpanzee mitochondrial sequences is about half that between human and gibbon sequences.)

2. The highlighted nodes belong to order Anura, which contains all frogs and toads. The Anura samples are distributed in 
several distinct clusters. There are also multiple non-Anura species, such as Caecilians, reptiles, rodents, 
interspersed within the Anura-major clades.

![anura](panmama/mito_panmat/image-3.png)

I then compared tree topology between the mito panMAT and the [TimeTree](https://timetree.org/). I was able to identify 
6,231 overlapping unique  species between the two trees. After removing duplicate samples and excluding non-overlapping 
species, I calculated a  normalized Robinson-Foulds distance of 0.302 (unweighted RF distance: 7,380; maximum RF 
distance: 24,456).

I wonder if these discrepancies are due to differences in the recorded starting positions of the circular mitochondrial 
genome. In many MAFFT alignments I have inspected, both samples have very similar sequence length, but one sample often 
has an extended gap sequence at the beginning of the aligned sequence while the other sample has an extended gap 
sequence at the end, which is what you'd expect if two circular genomes are linearized at different positions. 

I did a quick experiment on the human mitochondrial sample (J01415.2) and the chimpanzee sample (D38113.1) mentioned 
above. MAFFT alignment of the original NCBI fastas showed that the two fastas differ by 1380 SNPs and 1160 gaps. I then
manually rotated the human mitochondrial sequence by moving the region at the beginning where the chimpanzee sequence
contains a stretch of gaps to the end of the sequence. MAFFT alignment of the rotated human fasta showed 1,456 SNPs and
26 gaps, or ~8.6% difference, which is consistent with the [literature](https://pubmed.ncbi.nlm.nih.gov/8919866/). 

I also confirmed from the v_mtdna.aln.gz file that neither human mito sample nor the chimpanzee sample were rotated for
the construction of the panMAT.

## 11/25/2025

### To-do

- [ ] Investigate the phylogeny of nodes that reads are parsimonious to on the 8M tree
  - [ ] If it looks good: pursue phylogenetic information approach mentioned above
- [ ] implement syncmer only index
- [x] ~~Restore EM function on panMAMA~~
  - [x] ~~Modify WEPP to also use shotgun reads~~
  - [x] ~~Test if WEPP can run with HIV shotgun reads and compare to panMAMA~~

### Restore EM function on panMAMA

I restored the EM function on panMAMA. For now, I will the option to either use overlap coefficients or read-weighted
node scores to select probable nodes.

## 11/26/2025

### To-do

- [ ] Investigate the phylogeny of nodes that reads are parsimonious to on the 8M tree
  - [ ] If it looks good: pursue phylogenetic information approach mentioned above
- [x] ~~implement syncmer only index~~
- [x] ~~Restore EM function on panMAMA~~
  - [x] ~~Modify WEPP to also use shotgun reads~~
  - [x] ~~Test if WEPP can run with HIV shotgun reads and compare to panMAMA~~

### Fixing l parameter

Index building broke when I set l to 1 or 2... Due to some conceptual misunderstanding, I thought l would have to be at
least 3, so I coded the index building function to assume l >= 3...

Now index also builds correctly for l = 1 and l = 2.  l = 1 is the equivalent as implenting syncmer only index.

## 11/30/2025

I will make new to-do list for the next couple of weeks

- [ ] Formalize evaluation metrics
- [x] ~~Update seed sketching function to take in any seed parameters~~
- [ ] Compare accuracy of different combinations of l, k, s using the evaluation metric
- [ ] Implement t offsets for syncmers
- [ ] Work on SAPP application
- [ ] Look into ancient DNA read error correction

## 12/1/2025


### Build indices with different combinations of l, k, s parameters

I will write a skeleton slurm script that will generate combinations of indices for panMAMA.

### Implement t offsets for syncmers

I have also just modified the seed sketching function to sketch syncmers with all l, k, s, t, and open values.

## 12/2/2025

### Running panmap on simulated mito reads

Simulate reads

```
sbatch /private/groups/corbettlab/alan/lab_notebook/panmama/v_mtdna/gendata.sh \
  /data/tmp/ \
  /private/groups/corbettlab/alan/lab_notebook/panmama/v_mtdna/input_data/v_mtdna.panman \
  /private/groups/corbettlab/alan/lab_notebook/panmama/v_mtdna/input_data/v_mtdna.pmai
```

Run panmap

```
sbatch /private/groups/corbettlab/alan/lab_notebook/panmama/v_mtdna/run_panmap.sh \
  /private/groups/corbettlab/alan/lab_notebook/panmama/v_mtdna/out \
  /private/groups/corbettlab/alan/panmap/ \
  /private/groups/corbettlab/alan/lab_notebook/panmama/v_mtdna/input_data/v_mtdna.panman \
  /private/groups/corbettlab/alan/lab_notebook/panmama/v_mtdna/input_data/v_mtdna.pmai \
  /private/groups/corbettlab/alan/lab_notebook/panmama/v_mtdna/sim_reads
```

I've also just added a `simulate_and_run.py` (real path:
`/private/groups/corbettlab/alan/lab_notebook/panmama/benchmark/simulate_and_run.py`) that will simulate reads and run
panMAMA then delete all intermediate files. 

Run 

```
python3 ../../benchmark/simulate_and_run.py \
  ../input_data/v_mtdna.panman \
  ../input_data/v_mtdna.pmai \
  1 \
  --selection node_scores
```

### Investigating ancient DNA samples from Bianca

Out of the 652 samples that we got from Bianca, 520 samples do not align to the mammoth mito at all.

I grabbed a random sample `ar1_10.mammoth.fq` and tried it on panMAMA.

I first aligned it to a mammath sample `JF912199.1` on the tree and found that only two reads align to it but the 
alignment quality is very good. One read perfectly matches with the reference while the other one only has one SNP
difference.

However, the node scores and the abundance estimation doesn't look too good... Either the sample is very mixed with 
other haplotypes or the program is not performing correctly. I grabbed the most abundant haplotype from the output,
`MW232459.1`, and used `bwa` to align all the reads to it, and none of the reads mapped to it. I then aligned the reads
to other haplotypes and most reads that do align appear to have very low complexity, or highly repetitive, aka junk.

I think I will do some pre-processing to remove low-complexity reads. Additionally, I should also toss out reads whose
maximum parsimony is below a certain threshold.

## 12/3/2025

I will try `bbduk` to remove reads with low complexity

Using `bbduk` and setting entropy threshold to `0.7` seem good enough:

```
~/tools/BBTools/bbduk.sh in=ar1_10.mammoth.fq out=filtered.fq outm=low_complexity.fq entropy=0.7 
```

I will set a paremeter to throw away reads whose maximum parsimony score is less than a specific threshold (default to
0.5 of the total seeds).

Now with `ar1_10.mammoth.fq`, two reads remain after removing reads low maximum parismony score, and they are reads from
the mammoth genomes

I will try with another sample, `cr5_10.mammoth.fq`... Also worked. Very good. I will now go on phoenix and run the 
pipeline with all the reads from Bianca.

`panmama_ancient_mito.sh` maps all the read samples to the vertebrate mito tree. Outputs are in `ancient_mito_out/`. Job
array id: `22686096`.

`panmama_ancient_mammoth_mito.sh` maps all the read samples to the mammoth mito tree. Outputs are in
`ancient_mammoth_mito_out/`. Job array id: `22687109`.

While many samples were estimated to only contain mammoth mito, there still seem to be some contaminants not filtered 
out. I will also have panmama output the assigned reads to each haplotype estimated to be present, as well as most 
information on the parsimony scores of the reads remaining after the low-score removal step.

Rerun `panmama_ancient_mito.sh` (`22756229`) and `panmama_ancient_mammoth_mito.sh` (`22756235`)

## 12/4/2025

### Investigate ancient_mito_out/

I may have to manually look at each read that's assigned outside of the *Elephantidae* family. I will first filter out
the good outputs where all reads are placed to the *Elephantidae* family.

```
cd /private/groups/corbettlab/alan/lab_notebook/panmama/v_mtdna
python3 gen_stats_mito_out.py ancient_mito_out
```

#### Closer look at sample ar1_2

When I was manually investigating samples with reads placing parsimoniously outside of the mammoth "clade", I found this
sample `ar1_2` that has some reads that place to other non-mammoth Elepahntidae nodes, specifically the genus
*Loxodonta*. This may not may not indicate that this sample contains reads from a mixture of haplotype because the reads
might all map better to an ancestor of both *Loxodonta* and *Mammathus*. Will look into this further when we have a
better tree.

## 12/5/2025

### Investigate ancient_mammoth_mito_out/

I visually inspected the outputs of several samples, and panMAMA and pathphynder have either the same placement or
very closer placement. Since pathphynder labels the node differently from panMAN, I will need to first convert the node
names before I can write a script to compare the output of all the samples.

## 12/8/2025

Sumit just sent us the newly improved vertebrate mito panman. I also found and fixed a bug in the early exit of panMAMA
when no reads have significant placement score. Rerunning panMAMA on the mammoth tree (`22920510`) and the new verberate
tree (`22920457`).

### Mammoth tree results

Among cases where both PanMAMA and Pathphynder produced results, PanMAMA and Pathphynder identified identical nodes in
approximately 70% of samples, and nodes differing by a single branch in 14% of samples.

![compare_panmama_pathphynder_mammoth.png](panmama/v_mtdna/compare_panmama_pathphynder_mammoth.png)

I will also try to rn pathphynder myself to see how fast it is.

### Vertebrate tree results

Refering to `panmama/v_mtdna/mito_assignment_stats.tsv`, I still see what seem to be either mixed or unfiltered read
contaminations in the panMAMA results using the vertebrate mito tree. I think the latter is probably true. I will take a
deeper look at what the reads assigned outside of the mammoth "clade" look like.

### Real data to look at

Look for key words like: eDNA, metagenome, metagenomics, and commonly targeted mitochondrial regions like 16S, COI,
Cyt b and 12S

COI, 12S, Cyt b can better discriminate metazoan species than 16S as 16S can also target bacteria and swampy the data
with bacterial reads. On the other hands, we can try to use panMAMA to filter 16S reads for metazoan or organisms of
interest.

From Claude:

> Fish-focused surveys: MiFish primers (12S) are the clear winner
>
> General vertebrate diversity (fish, amphibians, mammals, birds): 12S with broader primers (e.g., Riaz primers, Vert-12S)
>
> Terrestrial vertebrates: 12S or 16S can work; some researchers prefer 16S for mammals
>
> High-resolution species identification with tissue samples: COI or Cyt b
>
> Multi-kingdom surveys: 18S (nuclear, eukaryotic universal) or combined 16S bacterial + 12S/18S eukaryotic approaches

##### SRX27728057: Metagenomic sequencing of Homo sapiens: blood microbiome in HIV patients (pre-ART)

Are we able to filter out the HIV reads from the rest of the reads we don't care able?

##### SRX23959018: metabarcoding of fisheries samples: eDNA in catch water to estimate species composition

Can look at the entire bioproject. Filter out junk and estimate what's in the fisheries?

## 12/9/2025

I tried different combinations of `--discard` parameter in panMAMA (discard reads with placement scores lower than a
threshold) and `entropyk` value in `bbduk` (kmer window size to measure entropy) and see which one produces better
read assignments. I've also added more information on the intersection of assigned reads between panMAMA and
Pathphynder.

I will come with a metric to measure the performance and try out different hyperparameters for maximum performance on 
both the vertebrate tree and the mammoth tree.

## 12/10/2025

I calculated the precision, recall, and F1 scores of the different `--discard` parameters in panMAMA and `entropyk`
values in `bbduk` 

I tried `entropyk` value at 4, 5, and 6 for `--dicard` threshold at 0.4, 0.5, 0.6, 0.7, 0.8.


![hyperparameter](panmama/v_mtdna/mito_assignment_stats_comparison.png)

It seems like an `entropyk=5` is the obvious winner. Between a discard threshold 0.6 and 0.8, there's a tradeoff
between precision and recall, which is not surprising. **F1 score shows that a discard threshold of 0.7 has a slight
edge.**

Precision-recall curve using `entropyk=5`

![prc](panmama/v_mtdna/mito_assignment_entropyk5_prc.png)

Below is the raw data

```
MPS  Entropyk  Precision            Recall               F1
0.4  4         0.11351030110935023  0.315702479338843    0.16698236922628587
0.4  5         0.29283069673510603  0.9613259668508287   0.44891640866873067
0.4  6         0.07950463225562489  0.4631057268722467   0.13571082781991287
0.5  4         0.260543580131209    0.3068432671081678   0.2818043588443994
0.5  5         0.6777950310559007   0.9635761589403974   0.7958067456700091
0.5  6         0.15223097112860892  0.4476295479603087   0.22719641857862338
0.6  4         0.3025302530253025   0.30353200883002207  0.30303030303030304
0.6  5         0.8676028084252758   0.9547461368653422   0.9090909090909091
0.6  6         0.24492931776275353  0.4403314917127072   0.3147709320695103
0.7  4         0.3093858632676709   0.2950276243093923   0.3020361990950226
0.7  5         0.9519230769230769   0.9298342541436464   0.9407490217998882
0.7  6         0.3572080291970803   0.432357813362783    0.3912065950537097
0.8  4         0.31840796019900497  0.2830292979546711   0.2996780801872988
0.8  5         0.9728395061728395   0.871199557766722    0.9192184310294547
0.8  6         0.4380069524913094   0.417910447761194    0.4277227722772277
```

```
MPS  Entropyk  Precision             Recall              F1
0.3  5         0.057291321903753224  0.9526692350027518  0.10808279479254472
0.4  5         0.29283069673510603   0.9613259668508287  0.44891640866873067
0.5  5         0.6777950310559007    0.9635761589403974  0.7958067456700091
0.6  5         0.8676028084252758    0.9547461368653422  0.9090909090909091
0.7  5         0.9519230769230769    0.9298342541436464  0.9407490217998882
0.8  5         0.9728395061728395    0.871199557766722   0.9192184310294547
0.9  5         0.9792592592592593    0.7340366463076069  0.8390986988257696
```


There might be some reads that could be in the Elephantidae family but are not mappable if a Mammuthus primigenius
reference is used. I will need to look at the sample individually. 

I will run the same analysis again when Sumit sent the a vertebrate tree with more primigenius reads merged in.

I'm starting to look at the sample individually to look for reads potentially from other Elephantidae family members.

| Interesting samples |
| ------------------- |
| ar1_2 |


Actually, another feature I can add for read filtering is throwing away reads mapping to spurious nodes on the tree with
node phylogenetic siginificance. Will try it out maybe tomorrow. And maybe skip over the EM step entirely for filtering
and assigning reads to divergent trees like this here?

## 12/11/2025

I think I will try to make a pretty figure today. Strategy:

1. Condense the tree by families
2. subsample the family tree
3. plot the family tree and highlight the Elephantidae family
4. plot the detailed Elephantidae family tree and annotate with reads assigned

## 12/12/2025

### Running and comparing WEPP to panMAMA

Command to enter WEPP docker

```
docker run -it   -v /private/groups/corbettlab/alan/lab_notebook/panmama/WEPP/WEPP/src/WEPP/qc_preprocess.py:/WEPP/src/WEPP/qc_preprocess.py   -v /private/groups/corbettlab/alan/lab_notebook/panmama/WEPP/WEPP/config/config.yaml:/WEPP/config/config.yaml   -v /private/groups/corbettlab/alan/lab_notebook/panmama/WEPP/WEPP/workflow/rules/qc.smk:/WEPP/workflow/rules/qc.smk   -v /private/groups/corbettlab/alan/lab_notebook/panmama/WEPP/WEPP/src/WEPP/config.hpp:/WEPP/src/WEPP/config.hpp   -v /private/groups/corbettlab/alan/lab_notebook/panmama/WEPP/for_pranav:/data   pranavgangwar/wepp:latest
```

Use this command to get the results comparing panMAMA's and WEPP's performance on mixed HIV samples with 50 haplotypes:

```
for rep in 0 1 2 3 4 5 6 7 8 9; do echo $rep; sed 's/WEPP/MAMA/g'  panmama/hiv20000_50hap-a_300000_rep${rep}.mama.distance_sum.tsv; log="wepp/results/hiv20000_50hap-a_300000_rep${rep}/hiv_50hap_rep${rep}_run.txt"; whd=$(grep 'true_node' $log | cut -f3 -d ' ' | awk '{sum += $1 * 0.02} END {print sum}'); wpd=$(grep 'PEAK' $log | cut -f 6 -d ' ' | awk '{sum += $1} END {print sum}'); echo -e "WEPP\t${whd}\t${wpd}"; echo ""; done
```

## 3/30/2026

We just submitted Panmap to Nature Genetics!

I haven't updated the lab notebook in a while. Now that I will start a new project,I will try to update it more
consistently from now on.

### Panmap

I updated `/scratch1/alan/panmap/build/_deps/panman-src/src/panman.cpp` to use the new PanMAN format.

#### Ancient DNA

I will run panmap on the ancient DNA for selected Salicaceae chloroplast genomes.

I just implemented a new function for 
fast assignment of reads to their LCA nodes on the tree. Now panmap also outputs the LCA read counts for each node.

The Salicaceae output is a bit confusing. Will need to investigate further.

For tomorrow: FIX BUG IN PLOTTING SCRIPT. NCXXXXX nodes are not being correctly labeled.

## 3/31/2026

### Panmap

I updated the panmap CMakeLists.txt to build on Phoenix.

A backup of the new CMakeLists.txt is in `/backup/3_31_2026/CMakeLists.txt`.

### Ancient DNA

I rerooted the salicaceae panman tree and ran Zihao's data again.


## 4/17/2026

### Chloroplast panMAN

cpstools IR had some trouble identifying the regions for a few of the genomes. For one, it flagged the assembly as
incorrect; cutting the first 50 lines of its fasta and pasting them to the end fixed it. For 18 others, there were
ambiguous bases that cpstools refused to handle, so I replaced each with the lowest lexicographical canonical base it
represents, solely for the purpose of region identification. I then ran cpstools Seq to orient the original fastas using
the resulting region coordinates. I've also attached information on the ambiguous sites so you can decide whether to
discard genomes with too much missing data or with ambiguous sites too close to the region boundaries.

I used `OP650215.1` as the seed to orient the regions in the same direction using minimap2 and seqkit.

Used mafft to align the regions. I'm using mafft for now because it's fast and easy to use. I decided to do both mafft
auto mode and ginsi mode (--globalpair) to see which one performs better. Ginsi mode will take much longer so I will go
ahead and build the tree using the auto mode results for now.

Trim the gappy boundries using `trimal` and also get the stats.

```bash
cd /scratch1/alan/lab_notebook/panmama/salicaceae/data/salicaceae/salicaceae_regions/split_regions/alignments
mkdir -p stats
for aln in *mafft-auto.aln; do
  prefix=$(basename $aln .aln)
    
  trimal -in $aln -out ${prefix}.trimmed.aln -automated1

  trimal -in $aln                  -sident > stats/${prefix}.sident.txt
  trimal -in ${prefix}.trimmed.aln -sident > stats/${prefix}.trimmed.sident.txt

  python3 /scratch1/alan/lab_notebook/panmama/salicaceae/data/salicaceae/salicaceae_regions/plot_identity_heatmap.py \
    stats/${prefix}.sident.txt \
    --metadata /scratch1/alan/lab_notebook/panmama/salicaceae/data/salicaceae/meta.tsv \
    --output stats/${prefix}.sident.heatmap.png

  python3 /scratch1/alan/lab_notebook/panmama/salicaceae/data/salicaceae/salicaceae_regions/plot_identity_heatmap.py \
    stats/${prefix}.trimmed.sident.txt \
    --metadata /scratch1/alan/lab_notebook/panmama/salicaceae/data/salicaceae/meta.tsv \
    --output stats/${prefix}.trimmed.sident.heatmap.png
done
```

In the `alignments` dir, I have the MSA aln files of ingroup, ingroup&outgroup (aligned together), and ingroup+outgroup 
(outgroup added to the ingroup alignment). I will use the ingroup&outgroup alignment to build the tree.

Oh.. forgot to add, the outgroup is a Viola biflora chloroplast genome (`OM177182.2`).

First, concatenate the alignments.

```bash
cd /scratch1/alan/lab_notebook/panmama/salicaceae/data/salicaceae/salicaceae_regions/concatenated/regions
python3 /scratch1/alan/lab_notebook/panmama/salicaceae/data/salicaceae/salicaceae_regions/concat_aln_regions.py \
  -i merged.LSC.oriented.with_outgroup.mafft-auto.trimmed.aln \
     merged.SSC.oriented.with_outgroup.mafft-auto.trimmed.aln \
     merged.IRa.oriented.with_outgroup.mafft-auto.trimmed.aln \
  -o ../concat.fasta \
  -p ../partition.txt
```

Then run `iqtree`
```
iqtree2 -s concat.fasta -p partition.txt -m MFP+MERGE -B 1000 -alrt 1000 -T AUTO &
```

Seems like 13 genomes failed the initial gap/ambiguity assessment... I will go ahead and build the tree with them 
included. I will try building the tree without them later.

```
grep failed partition.txt.log  | grep -v TOTAL
```

<details>
<summary>Failed genomes</summary>

```bash
   2  CM018592.1            14.70%    failed      0.00%
   8  MK341060.1            13.88%    failed      0.00%
   9  MK722343.1            16.35%    failed      0.00%
  10  MK748469.1            14.94%    failed      0.00%
  13  MN078141.1             3.59%    failed      1.40%
  21  MW147595.1             2.76%    failed      0.44%
  33  MW801246.1             3.57%    failed      1.43%
  38  NC_028350.1           14.83%    failed      0.00%
  68  NC_044462.1           11.28%    failed      0.00%
  71  NC_045919.1           16.37%    failed      1.70%
  72  NC_046687.1           11.26%    failed      0.00%
  87  OL622106.1            10.45%    failed      0.00%
 116  PQ400062.1             3.68%    failed      2.59%
```

</details>

The results are back

```bash
iqtree -rf partition.txt.treefile partition.txt.contree
column -t partition.txt.contree.rfdist
```

```
1      1
Tree0  0
```

Excellent both the ML tree (`.treefile`) and the bootstrap consensus tree (`.contree`) are identical.


Now reroot the tree to the outgroup then remove the outgroup from the tree.

```
nw_reroot iqtree_results/partition.txt.treefile OM177182.2.oriented > iqtree_results/partition.txt.rerooted.nwk
nw_prune iqtree_results/partition.txt.rerooted.nwk OM177182.2.oriented > iqtree_results/partition.txt.root_pruned.nwk
```

## 4/19/2026

Now I'm going to investigate the genomes that failed the initial gap/ambiguity assessment. Some of them have a very
small IR region and "abnormally" large SSC and LSC regions, while their total genome sizes are consistent with other
genomes... Did cpstools identify the regions incorrectly?

```bash
$ sort -k3,3 -g region_sizes_info.tsv | column -t | head
MK722343.1   Salix             7701   50551  89868   155821
NC_045919.1  Homalium          8135   21454  119392  157116
CM018592.1   Salix             9338   51477  85451   155604
NC_028350.1  Salix             9389   46603  91438   156819
MK748469.1   Salix             9880   50333  84933   155026
MK341060.1   Populus           10606  50453  84629   156294
NC_044462.1  Populus           13886  43905  84513   156190
OL622106.1   Salix             14618  41869  84432   155537
NC_046687.1  Flacourtia        14908  41078  85329   156223
NC_033876.1  Populus           16416  16806  105717  155355
```

Seems like `cpstools` did not identify the regions correctly for some of the genomes, at least for `MK722343.1` that I
checked manually using `mummer`. I manually rotated its fasta file for easier inspection. Here is its mummer coordinates
output:

```
[S1]    [E1]    [S2]    [E2]    [LEN 1] [LEN 2] [% IDY] [LEN R] [LEN Q] [COV R] [COV Q] [TAGS]
1       155821  1       155821  155821  155821  100.00  155821  155821  100.00  100.00  MK722343.1.fasta        MK722343.1.fasta
61780   89319   133095  105562  27540   27534   99.95   155821  155821  17.67   17.67   MK722343.1.fasta        MK722343.1.fasta
105562  133095  89319   61780   27534   27540   99.95   155821  155821  17.67   17.67   MK722343.1.fasta        MK722343.1.fasta
```

vs the cpstools output

```
MK722343.1.fasta        LSC:130419-155821,1-64465       IRb:64466-72166 SSC:72167-122717        IRa:122718-130418
```

While cpstools identified the IR regions to be 7701 bp, `mummer` identified the IR regions to be ~27540 bp, which is a 
much more expected size... Same is true for `NC_045919.1`.

mummer output:
```
[S1]    [E1]    [S2]    [E2]    [LEN 1] [LEN 2] [% IDY] [LEN R] [LEN Q] [COV R] [COV Q] [TAGS]
1       157116  1       157116  157116  157116  100.00  157116  157116  100.00  100.00  NC_045919.1.fasta       NC_045919.1.fasta
78750   106383  150671  123038  27634   27634   99.91   157116  157116  17.59   17.59   NC_045919.1.fasta       NC_045919.1.fasta
123038  150671  106383  78750   27634   27634   99.91   157116  157116  17.59   17.59   NC_045919.1.fasta       NC_045919.1.fasta
```

vs cpstools output

```
NC_045919.1.fasta       LSC:133573-157116,1-95848       IRb:95849-103983        SSC:103984-125437       IRa:125438-133572
```

I'm going to write a alignment (prolly nucmer) based script to identify the IR regions..

```
Warning: rotated genome still has IR regions near the edges for NC_044462.1_intermediates/rotated_NC_044462.1.fa
Warning: rotated genome still has IR regions near the edges for OK505606.1_intermediates/rotated_OK505606.1.fa
Warning: rotated genome still has IR regions near the edges for OL622106.1_intermediates/rotated_OL622106.1.fa
Warning: rotated genome still has IR regions near the edges for OQ791207.1_intermediates/rotated_OQ791207.1.fa
```

## 4/20/2026

I wrote a `iden_region_mummer.py` script that self-aligns a circular plastid genome FASTA with MUMmer's nucmer, rotates
it so the inverted repeats (IRs) aren't split across the boundary, and writes the four quadripartite regions (LSC, IRa,
SSC, IRb) as separate FASTAs. If the IRs are non-identical or multiple inversion pairs are detected, it logs a warning
for manual inspection.

Only two genomes require manual inspection: `NC_028350.1` and `NC_032368.1`.

```
$ grep 'Found more than a single inversion pair' *log
NC_028350.1.log:Found more than a single inversion pair in ../../salicaceae_fastas/original_fastas/NC_028350.1.fasta
NC_032368.1.log:Found more than a single inversion pair in ../../salicaceae_fastas/original_fastas/NC_032368.1.fasta
```

**NC_028350.1**

nucmer output:
```
[S1]    [E1]    [S2]    [E2]    [LEN 1] [LEN 2] [% IDY] [LEN R] [LEN Q] [COV R] [COV Q] [TAGS]
1       156819  1       156819  156819  156819  100.00  156819  156819  100.00  100.00  NC_028350.1     NC_028350.1
45258   59091   116817  102998  13834   13820   99.24   156819  156819  8.82    8.81    NC_028350.1     NC_028350.1
59196   72879   102885  89153   13684   13733   98.52   156819  156819  8.73    8.76    NC_028350.1     NC_028350.1
89153   102885  72879   59196   13733   13684   98.52   156819  156819  8.76    8.73    NC_028350.1     NC_028350.1
102998  116817  59091   45258   13820   13834   99.24   156819  156819  8.81    8.82    NC_028350.1     NC_028350.1
103745  116817  58344   45258   13073   13087   99.19   156819  156819  8.34    8.35    NC_028350.1     NC_028350.1
```

Seems like there's a small gap at 59092-59195 and 102886-102997, ~100 bps, while the rest of the alignment is very high
identity. The sizes of the inverted repeats are also consistent with the expected IR regions after merging... I will 
manually merge the IR regions to 45258-72879 and 89153-116817. 

```
NC_028350.1_final.fa  LSC:116818-156819,1-45257  IRb:45258-72879 SSC:72880-89152 IRa:89153-116817
```

**NC_032368.1**

nucmer output:
```
[S1]    [E1]    [S2]    [E2]    [LEN 1] [LEN 2] [% IDY] [LEN R] [LEN Q] [COV R] [COV Q] [TAGS]
1       158591  1       158591  158591  158591  100.00  158591  158591  100.00  100.00  NC_032368.1     NC_032368.1
44634   72301   118591  90924   27668   27668   100.00  158591  158591  17.45   17.45   NC_032368.1     NC_032368.1
70230   72301   72302   74373   2072    2072    100.00  158591  158591  1.31    1.31    NC_032368.1     NC_032368.1
72302   74373   70230   72301   2072    2072    100.00  158591  158591  1.31    1.31    NC_032368.1     NC_032368.1
72302   74373   92995   90924   2072    2072    100.00  158591  158591  1.31    1.31    NC_032368.1     NC_032368.1
90924   118591  72301   44634   27668   27668   100.00  158591  158591  17.45   17.45   NC_032368.1     NC_032368.1
90924   92995   74373   72302   2072    2072    100.00  158591  158591  1.31    1.31    NC_032368.1     NC_032368.1
```

Seems like it's got a moderate sized tandem repeat (~2kb) at the SSC and IRb boundry that's also inverted repeat... I
think I will set the tandem repeat to be the SSC region. This is also consistent with `cpstools` output:

```
NC_032368.1_final.fa    LSC:118592-158591,1-44633       IRb:44634-72301 SSC:72302-90923 IRa:90924-118591
```

### Rerunning alignment and iqtree using regions identified by `iden_region_mummer.py`

```
grep failed partition.txt.log  | grep -v TOTAL
```

<details>
<summary>Failed genomes</summary>

```bash
  13  MN078141.1     4.75%    failed      1.38%
  21  MW147595.1     4.04%    failed      0.23%
  33  MW801246.1     4.72%    failed      1.58%
 117  PQ400062.1     4.60%    failed      4.26%
```

</details>

## 4/21/2026

### Salicaceae tree

The alignment files and the final tree seem to be better with the regions identified by `iden_region_mummer.py`. The
3/4 genomes that failed the initial gap/ambiguity assessment are all from `Casearia` genus and the other one, MW147595,
is from `Dianyuea` genus and is the only one genome im the `Dianyuea`.

Seems possible that the `Casearia` genus might have some distinct IR boundry shifts. I will try to build a `Casearia` 
only tree with the outgtoup to see if the tree is consistent with the `Casearia` subtree in the Salicaceae tree.

I just confirmed that the `Casearia` only tree is consistent with the `Casearia` subtree in the Salicaceae tree. Yay.

### Primate Alu tree

First get all the Alu family consensus sequences from `famdb`

```
famdb.py -i famdb/ families -f fasta_name --include-class-in-name --class SINE/Alu -ad 9443 > primate_alu.ad.fa
```

I used different combinations of `-ad` and `--curated` to get

```
$ ls -1 *fa
primate_alu.ad.curated.fa
primate_alu.ad.fa
primate_alu.curated.fa
primate_alu.d.curated.fa
primate_alu.d.fa
primate_alu.fa
```

I will start with the `primate_alu.ad.curated.fa` file for now.

I think I will use two different methods to get the actual Alu sequences. First, through blastn against the core_nt database. Second, through the annotated Alu's from UCSC genome browser's rmsk table/track. 

#### BLASTn against core_nt database

Running `blastn` on phoenix.

```
echo "Running blastn..."
time /private/home/bzhan146/tools/blast/ncbi-blast-2.17.0+/bin/blastn \
  -task dc-megablast \
  -query "$fasta" \
  -db "$tmp_dir/ncbi_core_nt/core_nt" \
  -taxids 9443 \
  -evalue 1e-10 \
  -dust no \
  -perc_identity 80 \
  -max_target_seqs 10000 \
  -outfmt '6 qseqid sseqid staxid ssciname pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand qseq sseq' \
  -num_threads 16 \
  -out "$out_file"
```

See `bzhan146@emerald.prism:/private/home/bzhan146/scripts/tes/run_blast.sh` for the slurm script. Output is saved to 
`bzhan146@emerald.prism:/private/groups/corbettlab/alan/lab_notebook/tes/primate_alu.ad.curated.blast_hits.txt`

#### UCSC rmsk track

Get all the UCSC assemblies.

```
curl -s "https://api.genome.ucsc.edu/list/ucscGenomes" | \
  jq -r '.ucscGenomes | to_entries[] | 
         [.key, .value.organism, .value.scientificName, .value.taxId, .value.description] 
         | @tsv' \
  > all_ucsc_assemblies.tsv
```

Run script below to get all the **primate assemblies** from the UCSC assemblies and save to
`/scratch1/alan/lab_notebook/tes/ucsc_primate_assemblies.tsv`

```
python3 /scratch1/alan/lab_notebook/tes/scripts/get_primate_ids.py
```

Then I manually selected the latest assembly for eachp primate species and saved to
`/scratch1/alan/lab_notebook/tes/ucsc_primate_assemblies.selected.tsv`:

```
calJac4 Marmoset        Callithrix jacchus      9483    May 2020 (Callithrix_jacchus_cj1700_1.1/calJac4)
chlSab2 Green monkey    Chlorocebus sabaeus     60711   Mar. 2014 (Chlorocebus_sabeus 1.1/chlSab2)
gorGor6 Gorilla Gorilla gorilla gorilla 9595    Aug. 2019 (Kamilah_GGO_v0/gorGor6)
hg38    Human   Homo sapiens    9606    Dec. 2013 (GRCh38/hg38)
hs1     Human   Homo sapiens    9606    Jan. 2022 (T2T CHM13v2.0/hs1)
macFas5 Crab-eating macaque     Macaca fascicularis     9541    Jun. 2013 (Macaca_fascicularis_5.0/macFas5)
micMur2 Mouse lemur     Microcebus murinus      30608   May 2015 (Mouse lemur/micMur2)
nasLar1 Proboscis monkey        Nasalis larvatus        43780   Nov. 2014 (Charlie1.0/nasLar1)
nomLeu3 Gibbon  Nomascus leucogenys     61853   Oct. 2012 (GGSC Nleu3.0/nomLeu3)
otoGar3 Bushbaby        Otolemur garnettii      30611   Mar. 2011 (Broad/otoGar3)
panPan3 Bonobo  Pan paniscus    9597    May 2020 (Mhudiblu_PPA_v0/panPan3)
panTro6 Chimp   Pan troglodytes 9598    Jan. 2018 (Clint_PTRv2/panTro6)
papAnu4 Baboon  Papio anubis    9555    Apr. 2017 (Panu_3.0/papAnu4)
papHam1 Baboon  Papio hamadryas 9562    Nov. 2008 (Baylor Pham_1.0/papHam1)
ponAbe3 Orangutan       Pongo pygmaeus abelii   9601    Jan. 2018 (Susie_PABv2/ponAbe3)
rheMac10        Rhesus  Macaca mulatta  9544    Feb. 2019 (Mmul_10/rheMac10)
rhiRox1 Golden snub-nosed monkey        Rhinopithecus roxellana 61622   Oct. 2014 (Rrox_v1/rhiRox1)
saiBol1 Squirrel monkey Saimiri boliviensis     39432   Oct. 2011 (Broad/saiBol1)
```

Run to see what tracks matching 'rmsk' and 'repeat' are there for each assembly

```
for db in $(cut -f1 ucsc_primate_assemblies.selected.tsv); do
  tracks=$(curl -s "https://api.genome.ucsc.edu/list/tracks?genome=${db}" | \
           jq -r --arg db "$db" '
             .[$db] | keys[] | select(test("rmsk|repeat"; "i"))' | \
           paste -sd, -)
  echo -e "${db}\t${tracks:-none}"
done
```

```
calJac4 nestedRepeats,rmsk,simpleRepeat
chlSab2 nestedRepeats,rmsk,simpleRepeat
gorGor6 nestedRepeats,rmsk,simpleRepeat
hg38    joinedRmsk,nestedRepeats,rmsk,simpleRepeat
hs1     simpleRepeat,t2tRepeatMasker
macFas5 rmsk,simpleRepeat
micMur2 nestedRepeats,rmsk,simpleRepeat
nasLar1 nestedRepeats,rmsk,simpleRepeat
nomLeu3 nestedRepeats,rmsk,simpleRepeat
otoGar3 nestedRepeats,rmsk,simpleRepeat
panPan3 nestedRepeats,rmsk,simpleRepeat
panTro6 nestedRepeats,rmsk,simpleRepeat
papAnu4 nestedRepeats,rmsk,simpleRepeat
papHam1 nestedRepeats,rmsk,simpleRepeat
ponAbe3 nestedRepeats,rmsk,simpleRepeat
rheMac10        nestedRepeats,rmsk,simpleRepeat
rhiRox1 nestedRepeats,rmsk,simpleRepeat
saiBol1 nestedRepeats,rmsk,simpleRepeat
```

Seems like every core primate assembly has rmsk except hs1 (T2T-CHM13). The T2T assembly uses a different track name (t2tRepeatMasker) because it's served as a hub rather than a traditional MySQL-backed assembly. So hs1 needs a different
download path than the others.

Download the pre-computed RepeatMasker output files for assemblies with rmsk track (calJac4 chlSab2 gorGor6 hg38 macFas5
micMur2 nasLar1 nomLeu3 otoGar3 panPan3 panTro6 papAnu4 papHam1 ponAbe3 rheMac10 rhiRox1 saiBol1).

```
for db in calJac4 chlSab2 gorGor6 hg38 macFas5 micMur2 nasLar1 nomLeu3 \
          otoGar3 panPan3 panTro6 papAnu4 papHam1 ponAbe3 rheMac10 \
          rhiRox1 saiBol1; do
  url="https://hgdownload.soe.ucsc.edu/goldenPath/${db}/bigZips/${db}.fa.out.gz"
  echo "Fetching ${db}..."
  wget -c -q "$url" -O "${db}.fa.out.gz" || echo "  FAILED: $db"
done
```

Then download the T2T-CHM13 RepeatMasker file.

```
wget https://hgdownload.soe.ucsc.edu/hubs/GCA/009/914/755/GCA_009914755.4/bbi/GCA_009914755.4_T2T-CHM13v2.0.t2tRepeatMasker/chm13v2.0_rmsk.out.gz
```

Extract Alus positions to bed files.

```
for f in *.fa.out.gz chm13v2.0_rmsk.out.gz; do
  db="${f%.fa.out.gz}"
  zcat "$f" | awk 'NR>3 && $11 ~ /^SINE\/Alu$/ {
    chrom=$5; start=$6-1; end=$7; strand=($9=="C"?"-":"+");
    name=$10;
    print chrom"\t"start"\t"end"\t"name"\t0\t"strand
  }' > "${db}.alu.bed"
  echo "${db}: $(wc -l < ${db}.alu.bed) Alu elements"
done
```

Download the fasta files for the assemblies.

```
mkdir -p genomes && cd genomes
for db in calJac4 chlSab2 gorGor6 hg38 macFas5 micMur2 nasLar1 nomLeu3 \
          otoGar3 panPan3 panTro6 papAnu4 papHam1 ponAbe3 rheMac10 \
          rhiRox1 saiBol1; do
  wget -c -q "https://hgdownload.soe.ucsc.edu/goldenPath/${db}/bigZips/${db}.fa.gz" \
       -O "${db}.fa.gz"
done

wget -c "https://hgdownload.soe.ucsc.edu/hubs/GCA/009/914/755/GCA_009914755.4/GCA_009914755.4.fa.gz" \
     -O hs1.fa.gz
```

Decompress for bedtools

```
for f in *.fa.gz; do
  gunzip "$f"
done
```

Running bedtools to get the Alu sequences to fasta files.

```
for f in beds/*.bed; do
  db=$(basename "$f" .alu.bed)
  echo "Extracting ${db}..."
  bedtools getfasta \
    -fi "../genome_fastas/${db}.fa" \
    -bed "$f" \
    -s \
    -name+ \
    -fo "${db}.alu.fa" \
    2> "${db}.alu.log"
done
```

## 4/22/2026

I did some sanity check and then sent the tree and alignment to Zihao.

## 4/23/2026

### Salicaceae Tree

I compared my tree with Zihao's tree and alignment. They are highly similar but some nodes still have pretty low
support. I think I will go ahead and use my tree and alignment to build a panMAT.

### Pimate Alu Tree

#### BLAST

Blast finished running for the `primate_alu.ad.curated.fa`. Now I have a `primate_alu.ad.curated.blast_hits.txt` to
parse.

(Mostly) Claude and I wrote a script `filter_blast_out.py` to filter the blast hits:

Post-processes BLAST outfmt 6 output to remove redundant overlapping HSPs within each (query, subject) group, keeping
the highest-bitscore hit when two hits overlap beyond a configurable threshold. Applies optional pre-filters for minimum
alignment length and percent identity, and can treat plus/minus strand hits as either independent or competing. Streams
the input file for bounded memory usage regardless of file size.

```bash
python3 filter_blast_out.py primate_alu.ad.curated.blast_hits.txt \
  --min-pident 80 --min-length 200 -o primate_alu.hq.tsv --stats
```

#### UCSC rmsk track

I want to get a quick statistics on the Alu sequences from the UCSC rmsk track.

Get a quick counts on the Alu families for each assembly.

```bash
 for file in ../beds/*bed; do
  prefix=$(basename $file .alu.bed)
  cut -f4 $file | sort | uniq -c | awk -v OFS='\t' '{print $2,$1}' > ${prefix}.alu_counts.tsv &
done; wait &
```

Then merge the counts into a single file.

```bash
python3 -c "
import glob, csv
from collections import defaultdict

files = sorted(f for f in glob.glob('*alu_counts.tsv') if not f.startswith('merged'))
samples = [f.replace('.alu_counts.tsv', '') for f in files]
counts = defaultdict(dict)
for f in files:
    for row in csv.reader(open(f), delimiter='\t'):
        counts[row[0]][f] = row[1]

families = sorted(counts)
print('\t'.join(['Alu', 'total'] + samples))
for fam in families:
    vals = [counts[fam].get(f, '0') for f in files]
    total = sum(int(v) for v in vals)
    print('\t'.join([fam, str(total)] + vals))
" > merged_alu_counts.tsv
```

Ok I think I will start with AluYe5... It's a young, hominoid-biased AluY subfamily with ~12,400 copies concentrated in
great apes. It's a good pilot because the size is computationally tractable and the sequences retain strong
phylogenetic signal without saturation.

## 4/24/2026

### Hominoid AluYe5 Tree

As said yesterday, I will start with AluYe5...

Data filtering:

- **I will not include non-hominoid AluYe5 sequences here.** Seems like a small number of AluYe5 sequences are also
present in non-hominoid assemblies. They are probably mis-annotated;  

- **I Will skip over the Gibbon assembly.** The Gibbon assembly only has AluJb, AluSx, AluY, FLAM_C, and FRAM annotated,
andquite a lot of them too. It's likely annotated with a coarser Alu library that only recognized the top-level
subfamilies (AluJ, AluS, AluY). 

- **I will skip over the T2T assembly.** hg38 is probably better annotated. And since other assemblies are non-T2T, I
will keep it consistent... *NOTE (4/27/2026) Oops, I actually forgot to exclude the T2T assembly and have already done 
a lot of the downstream stuff... Ig I'll just keep it for now.*

Get all the AluYe5 sequences from the filtered assemblies.

```bash
for assembly in hg38 hs1 panTro6 panPan3 gorGor6 ponAbe3 nomLeu3; do
  seqkit grep -rp 'AluYe5::' alu_fastas/${assembly}.alu.fa >> hominoid_aluye5.fa
done
seqkit seq -uv hominoid_aluye5.fa > tmp && mv tmp hominoid_aluye5.fa
seqkit rmdup -s -P hominoid_aluye5.fa > tmp && mv tmp hominoid_aluye5.fa
```

Seems like most sequences are near full length. Gonna filter out the sequences that are less than 250 bp.

```bash
seqkit seq -g -m 250 hominoid_aluye5.fa  > hominoid_aluye5.len_filtered.fa
```

Still have a lot of sequences left. Good.

```bash
$ seqkit stats -a hominoid_aluye5.fa hominoid_aluye5.len_filtered.fa 2> /dev/null 
file                             format  type  num_seqs    sum_len  min_len  avg_len  max_len   Q1   Q2   Q3  sum_gap  N50  N50_num  Q20(%)  Q30(%)  AvgQual  GC(%)
hominoid_aluye5.fa               FASTA   DNA     10,677  2,793,007       11    261.6      373  280  297  305        0  300       53       0       0        0  54.46
hominoid_aluye5.len_filtered.fa  FASTA   DNA      8,526  2,548,063      250    298.9      373  294  301  307        0  301       52       0       0        0  54.79
```

Gonna rename the headers so downstream tools don't get confused by special characters.

```bash
seqkit replace -p '.+' -r 'AluYe5_{nr}' hominoid_aluye5.len_filtered.fa > hominoid_aluye5.len_filtered.renamed.fa
```

Also need to trim off the polyA tail from the sequences. I'll just use mafft to align each sequence to the consensus to
find the polyA tails then trim them off using a simple python script `trim_polya.py`.

```bash
ls hominoid_aluye5_split/*.fa | parallel --jobs 64 '
  prefix=$(basename {} .fa)
  mafft --auto --thread 1 <(cat aluye5.fa {}) 2>/dev/null > mafft_out/${prefix}.mafft.aln
'

parallel -j 64 'f={}; p=$(basename "$f" .mafft.aln); python3 trim_polya.py "$f" -o "trimmed/${p}.polyATrimmed.fa"' ::: mafft_out/*
cat trimmed/AluYe5_*polyATrimmed.fa > hominoid_aluye5.trimmed.fa
```


I should also remove sequences that are too divergent from the consensus sequence... these could be misannotated
sequences from older Alus.

Use bwa mem to align the sequences to the consensus sequence then discard sequences that fail to align and secondary
aligments.

```bash
bwa index consensus/aluye5.polyATrimmed.fa
bwa mem -k15 -t 32 consensus/aluye5.polyATrimmed.fa hominoid_aluye5.trimmed.fa  | samtools view -h -F 2308 - > aluye5_to_consensus.mapped.sam
samtools fastq aluye5_to_consensus.mapped.sam | seqkit fq2fa > aluye5_to_consensus.mapped.fa
```

Since the sequences are quite short, I will use mafft to align the sequences then trimal. 

```bash
mafft --auto aluye5_to_consensus.mapped.fa --thread 32 > aluye5_to_consensus.mapped.auto.aln
mafft --retree 2 --maxiterate 1000 --thread 32 aluye5_to_consensus.mapped.fa > aluye5_to_consensus.mapped.retree.aln
trimal -in aluye5_to_consensus.mapped.auto.aln -out aluye5_to_consensus.mapped.auto.trimmed.aln -automated1
trimal -in aluye5_to_consensus.mapped.retree.aln -out aluye5_to_consensus.mapped.retree.trimmed.aln -automated1
```

Run iqtree on both the untrimmed and trimmed sequences just to see the difference.

```bash
iqtree -s aluye5_to_consensus.mapped.auto.aln -m MFP -B 1000 -T 32 --prefix iqtree_results/aluye5_to_consensus.mapped.auto
iqtree -s aluye5_to_consensus.mapped.auto.trimmed.aln -m MFP -B 1000 -T 32 --prefix iqtree_results/aluye5_to_consensus.mapped.auto.trimmed
```

> [!WARNING]
> Deduplicate exact substring? Tempted but unsure.

## 4/26/2026

### Hominoid AluYe5 Tree

iqtree is still running for both the untrimmed and trimmed sequences, and it seems like it's going to take a while. In
the meantime, I will also try DIPPER to build the tree. For completeness, I will also align the sequences with TWILIGHT.

Build a guide tree using DIPPER then run TWILIGHT and trim the alignment using trimal.

```bash
dipper_cpu -i r -I aluye5_to_consensus.mapped.fa -O dipper_results/aluye5_to_consensus.mapped.dipper.guide.nwk
twilight -t dipper_results/aluye5_to_consensus.mapped.dipper.guide.nwk -i aluye5_to_consensus.mapped.fa -o aluye5_to_consensus.mapped.twilight.aln -v -w --check -r 0.999  --cpu-only -C 32
trimal -in aluye5_to_consensus.mapped.twilight.aln -out aluye5_to_consensus.mapped.twilight.trimmed.aln -automated1
```

Convert bases to uppercase then run DIPPER on all the aligned and trimmed sequences (mafft-auto, mafft-auto-trimmed,
mafft-retree, mafft-retree-trimmed, twilight, twilight-trimmed).

```bash
seqkit seq -u -w 0 -o aluye5_to_consensus.mapped.auto.upper.aln aluye5_to_consensus.mapped.auto.aln
seqkit seq -u -w 0 -o aluye5_to_consensus.mapped.auto.trimmed.upper.aln aluye5_to_consensus.mapped.auto.trimmed.aln
seqkit seq -u -w 0 -o aluye5_to_consensus.mapped.retree.upper.aln aluye5_to_consensus.mapped.retree.aln
seqkit seq -u -w 0 -o aluye5_to_consensus.mapped.retree.trimmed.upper.aln aluye5_to_consensus.mapped.retree.trimmed.aln
seqkit seq -u -w 0 -o aluye5_to_consensus.mapped.twilight.upper.aln aluye5_to_consensus.mapped.twilight.aln
seqkit seq -u -w 0 -o aluye5_to_consensus.mapped.twilight.trimmed.upper.aln aluye5_to_consensus.mapped.twilight.trimmed.aln

for file in *upper.aln; do prefix=$(basename $file .upper.aln); mv $file ${prefix}.aln; done

dipper_cpu -i m -I aluye5_to_consensus.mapped.auto.aln -O dipper_results/aluye5_to_consensus.mapped.auto.nwk -d 4 --threads 32
dipper_cpu -i m -I aluye5_to_consensus.mapped.auto.trimmed.aln -O dipper_results/aluye5_to_consensus.mapped.auto.trimmed.nwk -d 4 --threads 32
dipper_cpu -i m -I aluye5_to_consensus.mapped.retree.aln -O dipper_results/aluye5_to_consensus.mapped.retree.nwk -d 4 --threads 32
dipper_cpu -i m -I aluye5_to_consensus.mapped.retree.trimmed.aln -O dipper_results/aluye5_to_consensus.mapped.retree.trimmed.nwk -d 4 --threads 32
dipper_cpu -i m -I aluye5_to_consensus.mapped.twilight.aln -O dipper_results/aluye5_to_consensus.mapped.twilight.nwk -d 4 --threads 32
dipper_cpu -i m -I aluye5_to_consensus.mapped.twilight.trimmed.aln -O dipper_results/aluye5_to_consensus.mapped.twilight.trimmed.nwk -d 4 --threads 32
```

Then build a PanMAT out of each tree and alignment. NOTE: For .trimmed.nwk tree, use the untrimmed alignment.

```bash 
/scratch1/alan/panmap/build/bin/panmanUtils -M aluye5_to_consensus.mapped.auto.aln -N dipper_results/aluye5_to_consensus.mapped.auto.nwk -o aluye5_to_consensus.mapped.auto.panman
/scratch1/alan/panmap/build/bin/panmanUtils -M aluye5_to_consensus.mapped.auto.aln -N dipper_results/aluye5_to_consensus.mapped.auto.trimmed.nwk -o aluye5_to_consensus.mapped.auto.trimmed.panman
/scratch1/alan/panmap/build/bin/panmanUtils -M aluye5_to_consensus.mapped.retree.aln -N dipper_results/aluye5_to_consensus.mapped.retree.nwk -o aluye5_to_consensus.mapped.retree.panman
/scratch1/alan/panmap/build/bin/panmanUtils -M aluye5_to_consensus.mapped.retree.aln -N dipper_results/aluye5_to_consensus.mapped.retree.trimmed.nwk -o aluye5_to_consensus.mapped.retree.trimmed.panman
/scratch1/alan/panmap/build/bin/panmanUtils -M aluye5_to_consensus.mapped.twilight.aln -N dipper_results/aluye5_to_consensus.mapped.twilight.nwk -o aluye5_to_consensus.mapped.twilight.panman
/scratch1/alan/panmap/build/bin/panmanUtils -M aluye5_to_consensus.mapped.twilight.aln -N dipper_results/aluye5_to_consensus.mapped.twilight.trimmed.nwk -o aluye5_to_consensus.mapped.twilight.trimmed.panman
```

Get some stats on the PanMATs

```bash
(for file in *panman; do   echo "=== FILE: $file ===";   /scratch1/alan/panmap/build/bin/panmanUtils -s "$file" 2>/dev/null; done) | awk '
BEGIN {
  OFS="\t"
  n_fields = split("file,nodes,samples,substitutions,insertions,deletions,inversions,max_depth,mean_depth,block_ins,block_del,block_inv,block_dup,block_trans", headers, ",")
  for (i = 1; i <= n_fields; i++) printf "%s%s", headers[i], (i==n_fields ? "\n" : OFS)
}
/^=== FILE: / { if (file != "") print_row(); reset(); file = $3; next }
/^Total Nodes in Tree:/         { nodes = $5 }
/^Total Samples in Tree:/       { samples = $5 }
/^Total Substitutions:/         { subs = $3 }
/^Total Insertions:/            { ins = $3 }
/^Total Deletions:/             { dels = $3 }
/^Total Inversions:/            { invs = $3 }
/^Max Tree Depth:/              { maxd = $4 }
/^Mean Tree Depth:/             { meand = $4 }
/^Total Block Insertions:/      { bins = $4 }
/^Total Block Deletions:/       { bdel = $4 }
/^Total Block Inversion:/       { binv = $4 }
/^Total Block Duplications:/    { bdup = $4 }
/^Total Block Translocation:/   { btrans = $4 }
END { if (file != "") print_row() }
function print_row() {
  print file, nodes, samples, subs, ins, dels, invs, maxd, meand, bins, bdel, binv, bdup, btrans
}
function reset() {
  nodes=samples=subs=ins=dels=invs=maxd=meand=bins=bdel=binv=bdup=btrans=""
}
' > all_stats.tsv
```

Look at some stats

```bash
$ cut -f 1-6,8,9 all_stats.tsv  | column -t
file                                                nodes  samples  substitutions  insertions  deletions  max_depth  mean_depth
aluye5_to_consensus.mapped.auto.panman              16533  8267     92806          11219       16065      85         52.0646
aluye5_to_consensus.mapped.auto.trimmed.panman      16533  8267     96264          10687       16153      69         42.2331
aluye5_to_consensus.mapped.retree.panman            16533  8267     92520          10676       15661      78         43.4753
aluye5_to_consensus.mapped.retree.trimmed.panman    16533  8267     92854          10509       15798      73         45.0177
aluye5_to_consensus.mapped.twilight.panman          16533  8267     92958          12471       18506      76         44.3179
aluye5_to_consensus.mapped.twilight.trimmed.panman  16533  8267     93247          12354       18775      82         49.7084
```

Plot the stats

```bash
python panmat_stats.py panmans/all_stats.tsv -o panmat_stats
```

## 4/27/2026

### Salicaceae chloroplast tree

I plan to first build an individual panMAT for each chloroplast region then concatenate them together to a single 
PanMAT.

First build a panMAT for each chloroplast region using the iqtree and alignment files. Actually need to first remove the
outgroup, trim all gap columns, then do upper.

```bash
mkdir -p /scratch1/alan/lab_notebook/panmama/salicaceae/data/salicaceae/salicaceae_panMAT && cd /scratch1/alan/lab_notebook/panmama/salicaceae/data/salicaceae/salicaceae_panMAT

cp ../salicaceae_regions/mummer_regions/alignments/LSC_oriented.merged.mafft-auto.aln .
cp ../salicaceae_regions/mummer_regions/alignments/SSC_oriented.merged.mafft-auto.aln .
cp ../salicaceae_regions/mummer_regions/alignments/IRa_oriented.merged.mafft-auto.aln .
cp ../salicaceae_regions/mummer_regions/alignments/IRb_oriented.merged.mafft-auto.aln .

seqkit grep -v -p 'OM177182.2' IRa_oriented.merged.mafft-auto.aln > tmp && mv tmp IRa_oriented.merged.mafft-auto.aln
seqkit grep -v -p 'OM177182.2' IRb_oriented.merged.mafft-auto.aln > tmp && mv tmp IRb_oriented.merged.mafft-auto.aln
seqkit grep -v -p 'OM177182.2' LSC_oriented.merged.mafft-auto.aln > tmp && mv tmp LSC_oriented.merged.mafft-auto.aln
seqkit grep -v -p 'OM177182.2' SSC_oriented.merged.mafft-auto.aln > tmp && mv tmp SSC_oriented.merged.mafft-auto.aln

trimal -in IRa_oriented.merged.mafft-auto.aln -out tmp -noallgaps && mv tmp IRa_oriented.merged.mafft-auto.aln
trimal -in IRb_oriented.merged.mafft-auto.aln -out tmp -noallgaps && mv tmp IRb_oriented.merged.mafft-auto.aln
trimal -in LSC_oriented.merged.mafft-auto.aln -out tmp -noallgaps && mv tmp LSC_oriented.merged.mafft-auto.aln
trimal -in SSC_oriented.merged.mafft-auto.aln -out tmp -noallgaps && mv tmp SSC_oriented.merged.mafft-auto.aln

seqkit seq -u -w 0 -o IRa_oriented.merged.mafft-auto.upper.aln IRa_oriented.merged.mafft-auto.aln
seqkit seq -u -w 0 -o LSC_oriented.merged.mafft-auto.upper.aln LSC_oriented.merged.mafft-auto.aln
seqkit seq -u -w 0 -o SSC_oriented.merged.mafft-auto.upper.aln SSC_oriented.merged.mafft-auto.aln
seqkit seq -u -w 0 -o IRb_oriented.merged.mafft-auto.upper.aln IRb_oriented.merged.mafft-auto.aln

/scratch1/alan/panmap/build/bin/panmanUtils -M alignment/SSC_oriented.merged.mafft-auto.upper.aln -N ../salicaceae_regions/mummer_regions/concatenated/iqtree_results/partition.txt.root_pruned.nwk -o SSC_oriented &
/scratch1/alan/panmap/build/bin/panmanUtils -M alignment/LSC_oriented.merged.mafft-auto.upper.aln -N ../salicaceae_regions/mummer_regions/concatenated/iqtree_results/partition.txt.root_pruned.nwk -o LSC_oriented &
/scratch1/alan/panmap/build/bin/panmanUtils -M alignment/IRa_oriented.merged.mafft-auto.upper.aln -N ../salicaceae_regions/mummer_regions/concatenated/iqtree_results/partition.txt.root_pruned.nwk -o IRa_oriented &
/scratch1/alan/panmap/build/bin/panmanUtils -M alignment/IRb_oriented.merged.mafft-auto.upper.aln -N ../salicaceae_regions/mummer_regions/concatenated/iqtree_results/partition.txt.root_pruned.nwk -o IRb_oriented &
wait
```

Look at some stats

```bash
$ cut -f 1-6,8,9 all_stats.tsv  | column -t
file                 nodes  samples  substitutions  insertions  deletions  max_depth  mean_depth
IRa_oriented.panman  243    122      2023           1906        1598       23         12.959
IRb_oriented.panman  243    122      1975           1903        1682       23         12.959
LSC_oriented.panman  243    122      25477          18636       13699      23         12.959
SSC_oriented.panman  243    122      6930           3667        2040       23         12.959
```

I added some code to concatenate the panMATs into a single panMAN and opened a PR to the main panman repo.

### AluYe5 tree

I'm goign to make a meta file mapping the renamed aluye5 sequences back to their original headers and origin.

Refer back to [4/24/2026](#4242026) for the original parsing and filtering of the AluYe5 sequences.

The original headers (e.g. `AluYe5::chr1:1280182-1280482(+)`) do not record which assembly each sequence came from.
Re-scan each per-assembly file and emit a TSV with sequence ID and source assembly.

```bash
cd /scratch1/alan/lab_notebook/tes/alu_ucsc/by_family/aluye5
for assembly in hg38 hs1 panTro6 panPan3 gorGor6 ponAbe3 nomLeu3; do
  seqkit grep -rp 'AluYe5::' ../../alu_fastas/${assembly}.alu.fa \
    | seqkit seq -n -i \
    | awk -v a="$assembly" 'BEGIN{OFS="\t"} {print $1, a}'
done > hominoid_aluye5.assembly_map.tsv
```

Merge in species-level metadata from `ucsc_primate_assemblies.selected.tsv`.

```bash
awk -F'\t' 'BEGIN{OFS="\t"}
  NR==FNR { meta[$1] = $2 OFS $3 OFS $4 OFS $5; next }
  { print $1, $2, meta[$2] }
' /scratch1/alan/lab_notebook/tes/alu_ucsc/ucsc_primate_assemblies.selected.tsv hominoid_aluye5.assembly_map.tsv \
  > hominoid_aluye5.assembly_map.full.tsv
```

Join the rename table to the metadata table so each new ID carries its full metadata.

```bash
awk -F'\t' 'BEGIN{OFS="\t"}
  NR==FNR { meta[$1] = $0; next }
  $2 in meta { print $1, meta[$2] }
' hominoid_aluye5.assembly_map.full.tsv rename_id_map.tsv \
  > renamed_meta.tsv
```

## 4/28/2026

### AluYe5 tree

Looking at the AluYe5 tree on taxonium, it seems like there's a major cluster of nodes that are almost exclusively of
AluYe5 copies from the T2T assembly (hs1), while the rest of the nodes are interspersed with copies from hg38 and other
primate assemblies...

Biologically relevent? Not sure. Too early to say.

Did I remove mostly hg38 sequences when I did deduplication? Gonna add it back and see what's up. Also maybe I shouldn't
have removed the duplicate sequences.

Recover the deduplicated sequences.

```bash
cd /scratch1/alan/lab_notebook/tes/alu_ucsc/by_family/aluye5
for assembly in hg38 hs1 panTro6 panPan3 gorGor6 ponAbe3 nomLeu3; do
  seqkit grep -rp 'AluYe5::' ../../alu_fastas/${assembly}.alu.fa >> hominoid_aluye5_withdup.fa
done
seqkit seq -guv -m 250 hominoid_aluye5_withdup.fa  > hominoid_aluye5_withdup.len_filtered.fa
comm -13 <(seqkit seq -n hominoid_aluye5.len_filtered.fa | sort) <(seqkit seq -n hominoid_aluye5_withdup.len_filtered.fa | sort) > deduped_ids.txt
seqkit grep -n -f deduped_ids.txt hominoid_aluye5_withdup.len_filtered.fa > deduped_seqs.fa

# reorganize the files
mkdir removed && mv deduped_* removed/ && cd removed
```

Then do the filtering again (refer to [4/24/2026](#4242026), followed exactly).

Rename the headers again, starting from the last sequence from the original file.

```bash
seqkit replace -p '.+' -r 'AluYe5_{nr}' deduped.mapped.fa   | awk '/^>/{n++; sub(/[0-9]+$/, n+8526)} 1' > deduped.mapped.renamed.fa
```

Gonna align the duplicated sequences to the existing alignment files I've made using mafft with `--add --keeplength`.
(There are way faster ways to do this but this avoids potential buggy scripts.)

```bash
wdir=/scratch1/alan/lab_notebook/tes/alu_ucsc/by_family/aluye5/alignment/dup_added
mkdir -p $wdir && cd $wdir
cp ../../filtered/aluye5_to_consensus.mapped.fa  ../../removed/deduped.mapped.renamed.fa .
python3 add_dups.py \            
  --rep-fa  aluye5_to_consensus.mapped.fa \
  --dups-fa deduped.mapped.renamed.fa \
  --aln     ../aluye5_to_consensus.mapped.auto.aln \
            ../aluye5_to_consensus.mapped.retree.aln \
            ../aluye5_to_consensus.mapped.twilight.aln
mv ../*expanded.aln .
```

I just realized I should have deduplicated the alignment files (as they were polyA-deduplicated after the initial dedup)
before trimming then apply the same trimming to the original alignment files as duplicated sequences can bias towards
duplicated sequences.

I just wrote a simple script to do that. Gonna run it on the dedup-added alignment files as well as the original
alignment files.

```bash
for file in $(find . -maxdepth 2 -name "*aln" ! -path "./backup_alns*"); do
  prefix=$(echo $file | sed 's/.aln$//g')
  bash ~/tools/misc/trimal.sh $file ${prefix}.trimmed.aln &
done
wait
```

Run DIPPER and panmanUtils again.

```bash
for file in $(find . -maxdepth 2 -name "*aln" ! -path "./backup_alns*"); do
  prefix=$(basename $file .aln)
  dipper_cpu -i m -I $file -O ../dipper_results/${prefix}.nwk -d 4 --threads 32
done

cd /scratch1/alan/lab_notebook/tes/alu_ucsc/by_family/aluye5/panmans
for file in $(find ../alignment/ -maxdepth 2 -name "*aln" ! -name "*trimmed.aln" ! -path "../alignment/backup_alns*"); do
  prefix=$(basename $file .aln)
  /scratch1/alan/panmap/build/bin/panmanUtils -M $file -N ../dipper_results/${prefix}.nwk --threads 32 -o ${prefix}
  /scratch1/alan/panmap/build/bin/panmanUtils -M $file -N ../dipper_results/${prefix}.trimmed.nwk --threads 32 -o ${prefix}.trimmed
done
```

Now let's look at the tree again... Still looks the same? Well, the hs1 assembly does have ~5.8k AluYe5 copies
while hg38 has ~1.3k and other assemblies have ~1k. So maybe those extra 4k copies for some reason cluster together?

I'm gonna roughly get the substree of the hs1 AluYe5 clusters.

```bash
wdir=/scratch1/alan/lab_notebook/tes/alu_ucsc/by_family/aluye5/hs1_cluster
mkdir -p $wdir && cd $wdir
nw_clade ../dipper_results/aluye5_to_consensus.mapped.retree.expanded.trimmed.nwk AluYe5_3102 AluYe5_4351 > hs1_cluster.nwk
```

Get AluYe5 copies from hs1 inside the cluster, AluYe5 copies from hs1 but outside of the cluster, and hg38 AluYe5 copies.

```bash
seqkit seq -g ../alignment/dup_added/aluye5_to_consensus.mapped.retree.expanded.aln > aluye5_polyAtrimmed.fa
# Get the sequences of the hs1 AluYe5 clusters.
nw_labels -I hs1_cluster.nwk > hs1_cluster_labels.txt
sed -e 's/^/^/g' -e 's/$/\\s/g' hs1_cluster_labels.txt | \
  grep -f - ../renamed_meta_withdup.tsv | \
  awk '$3 == "hs1"' | \
  cut -f 1 | \
  seqkit grep -n -f - aluye5_polyAtrimmed.fa > hs1_incluster.fa

# Get the sequences of the hs1 AluYe5 copies outside of the cluster.
nw_prune ../dipper_results/aluye5_to_consensus.mapped.retree.expanded.trimmed.nwk $(cat hs1_cluster_labels.txt) > hs1_cluster_comp.nwk
nw_labels -I hs1_cluster_comp.nwk > hs1_cluster_comp_labels.txt
sed -e 's/^/^/g' -e 's/$/\\s/g' hs1_cluster_comp_labels.txt | \
  grep -f - ../renamed_meta_withdup.tsv  | \
  awk '$3 == "hs1"' | \
  cut -f1  | \
  seqkit grep -n -f - aluye5_polyAtrimmed.fa > hs1_outcluster.fa

# Get the sequences of the hg38 AluYe5 copies.
nw_labels -I ../dipper_results/aluye5_to_consensus.mapped.retree.expanded.trimmed.nwk > all_labels.txt
comm -12 <(awk '$3 == "hg38"' ../renamed_meta_withdup.tsv  | cut -f 1 | sort) <(sort all_labels.txt) > hg38_labels.txt
seqkit grep -n -f hg38_labels.txt aluye5_polyAtrimmed.fa > hg38.fa
```

Take a look at if they have differnet identity-to-consensus distributions

```bash
for file in hg38.fa  hs1_incluster.fa  hs1_outcluster.fa; do
   prefix=$(basename $file .fa)
   bwa mem -k15 ../consensus/aluye5.polyATrimmed.fa $file > ${prefix}.bam
   samtools calmd -b ${prefix}.bam ../consensus/aluye5.polyATrimmed.fa | samtools sort - > ${prefix}.md.bam
   samtools index ${prefix}.md.bam
   python3 extract_identity.py ${prefix}.md.bam ${prefix}.stats.tsv
done
python3 plot_comparison.py --samples S1=hg38.stats.tsv S2=hs1_incluster.stats.tsv S3=hs1_outcluster.stats.tsv 
```

![](tes/alu_ucsc/by_family/aluye5/hs1_cluster/identity_comparison.png)

I looked at the identity distributions across the clade. The hs1 copies within the hs1-only cluster show much lower
identity to the consensus compared to hs1 copies in the mixed clade. Interestingly, both hg38 copies and hs1 copies
outside the main hs1 clade have very similar identity distributions.

**chr14 short arm**

![](tes/alu_ucsc/by_family/aluye5/hs1_cluster/chr14_short_arm.igv.png)

**chr15 short arm**

![](tes/alu_ucsc/by_family/aluye5/hs1_cluster/chr15_short_arm.igv.png)

Viewing the positions of the AluYe5 copies on the hs1 on IGV also doesn't show any obvious acrocentric/centromeric
enrichment. ERRRR.

How closely related are the hs1 copies within the hs1-only cluster?

```bash
seqkit seq -n  hs1_incluster.fa  | seqkit grep -n -f - ../alignment/dup_added/aluye5_to_consensus.mapped.retree.expanded.aln > hs1_incluster.aln
python3 majority_consensus.py hs1_incluster.noallgaps.aln  incluster_consensus > hs1_incluster.consensus.fa
# align to incluster consensus as alignment steps above
python3 extract_identity.py bams/hs1_incluster.to_incluster_consensus.md.bam  bams/hs1_incluster.to_incluster_consensus.stats.tsv
tail -n+2 bams/hs1_incluster.to_incluster_consensus.stats.tsv | awk '{sum += $10; sum2 += $11} END {print sum/NR, sum2/NR}'
```

```
avg_iden_clip_ignored    avg_iden_clip_penalized
0.892901                 0.840315
```

It actually has lower identity than I'd expect...

Is it mislabeled? Could it be some other Alu families? I will use famdb to get the hmm profile of the primate alus then
run hmmer to the score of all the incluster sequences against all the alu family profiles to see if aluye5 actually is
the best match

Get the hmm profile of the primate alus from famdb and index it

```bash
famdb.py -i famdb/ families -f hmm --include-class-in-name --class SINE/Alu -ad 9443 --curated > primate_alu.ad.curated.hmm
hmmpress primate_alu.ad.curated.hmm 
```

Run nhmmscan to score the incluster sequences against the primate alu hmm profile

```bash
nhmmscan --cpu 32 --tblout hs1_incluster.hmm_hits.tblout --noali ../../../../primate_alus/hmm/primate_alu.ad.curated.hmm hs1_incluster.fa
nhmmscan --cpu 32 --tblout hs1_outcluster.hmm_hits.tblout --noali ../../../../primate_alus/hmm/primate_alu.ad.curated.hmm hs1_outcluster.fa
nhmmscan --cpu 32 --tblout all_polyAtrimmed.hmm_hits.tblout --noali ../../../../primate_alus/hmm/primate_alu.ad.curated.hmm aluye5_polyAtrimmed.fa
```

Yep, almost all of the hs1 copies within the hs1-only cluster are mislabeled as aluye5, while the hs1 copies outside of
the cluster are not mislabeled.

```console
$ awk '!/^#/ && !seen[$3]++ {print $3"\t"$1"\t"$14"\t"$13}' hs1_incluster.hmm_hits.tblout  | cut -f2 | sort | uniq -c | grep AluYe5 && grep '^>' hs1_incluster.fa | wc -l
      7 AluYe5#SINE/Alu
3828

$ awk '!/^#/ && !seen[$3]++ {print $3"\t"$1"\t"$14"\t"$13}' hs1_outcluster.hmm_hits.tblout  | cut -f2 | sort | uniq -c | grep AluYe5 && grep '^>' hs1_outcluster.fa | wc -l
    810 AluYe5#SINE/Alu
863
```

#### AluYxx

I think I will just get all the AluYxx sequences of all the primate assemblies and get their profile scores
against the primate alu hmm profiles.

Following steps recorded in [4/24/2026](#4242026).

Get all the AluYxx sequences of all the primate assemblies with length > 250.

```bash
wdir=/scratch1/alan/lab_notebook/tes/alu_ucsc/by_family/aluyxx
mkdir -p $wdir && cd $wdir
for fa in ../../alu_fastas/*.fa; do
  assembly=$(basename $fa .alu.fa)
  seqkit grep -rp 'AluY[^:]*::' $fa | seqkit replace -p '^' -r "${assembly}::" >> primate_aluyxx.fa
done
seqkit seq -uvg -m 250 primate_aluyxx.fa > primate_aluyxx.len_filtered.fa
```

Check how many sequences we have

```console
$ seqkit stats primate_aluyxx.fa primate_aluyxx.len_filtered.fa 2> /dev/null 
file                            format  type   num_seqs      sum_len  min_len  avg_len  max_len
primate_aluyxx.fa               FASTA   DNA   2,614,999  691,075,124       11    264.3      816
primate_aluyxx.len_filtered.fa  FASTA   DNA   2,057,453  617,776,798      250    300.3      816
```

Run nhmmscan to score the primate alu yxx sequences against the primate alu hmm profile. Gonna do it on both the curated
and the full primate alu hmm profile from famdb.

```bash
nhmmscan --cpu 10 --noali -o /dev/null \
  --tblout >(awk '!/^#/ { if ($3 != p) {c=0; p=$3} if (c++ < 5) print }' > aluyxx.top5.curated.tblout) \
  ../../../primate_alus/hmm/primate_alu.ad.curated.hmm primate_aluyxx.len_filtered.fa

# actually not gonna do this because it takes too long
nhmmscan --cpu 10 --noali -o /dev/null \
  --tblout >(awk '!/^#/ { if ($3 != p) {c=0; p=$3} if (c++ < 5) print }' > aluyxx.top5.tblout) \
  ../../../primate_alus/hmm/primate_alu.ad.hmm primate_aluyxx.len_filtered.fa
```

It still will take too long. I split up the fasta into 50 fastas and made a slurm script to run it on phoenix in
parallel. See `bzhan146@emerald.prism:/private/home/bzhan146/scripts/tes/run_blast.sh`.

## 4/29/2026

### Salicaceae tree

Concatenate the individual salicaceae panMANs into a single panMAN.

```bash
wdir=/scratch1/alan/lab_notebook/panmama/salicaceae/data/salicaceae/salicaceae_panMAT/panmans

cd $wdir
/scratch1/alan/panmap/build/bin/panmanUtils -C panman_list.txt -O orientation_list.txt -o salicaceae_concatenated
/scratch1/alan/panmap/build/bin/panmap salicaceae_concatenated.panman --index-mgsr salicaceae_concatenated.idx  -k 15 -s 8 -l 1
```

Rerun Zihao and Bianca's data using the new tree.

```bash
wdir=/scratch1/alan/lab_notebook/panmama/salicaceae

cd $wdir
mkdir from_bianca from_zihao

/scratch1/alan/panmap/build/bin/panmap data/salicaceae/salicaceae_panMAT/panmans/salicaceae_concatenated.panman \
  -i data/salicaceae/salicaceae_panMAT/panmans/salicaceae_concatenated.idx \
  --meta --filter-and-assign \
  --batch-files-path data/batch.txt \
  --discard 0.6 --dust 5 -t 16 

/scratch1/alan/panmap/build/bin/panmap data/salicaceae/salicaceae_panMAT/panmans/salicaceae_concatenated.panman \
  data/from_zihao/KapK.cpDNA.dedup.fastq.gz \
  -i data/salicaceae/salicaceae_panMAT/panmans/salicaceae_concatenated.idx \
  --meta --filter-and-assign \
  --discard 0.6 --dust 5 -t 16 \
  --output results/from_zihao/KapK.cpDNA
```

Now I'm going to write a script that generates a Taxonium/Nextclade json file for viewing the results.

Need to first make a template json file.

```bash
wdir=/scratch1/alan/lab_notebook/panmama/salicaceae/results
cd $wdir

# make a short meta file with non-redundant information
cut -f 1,2,6,7 ../data/salicaceae/meta.tsv  > ../data/salicaceae/meta.short.tsv

# make a template json file using taxonium tools
newick_to_taxonium -i ../data/salicaceae/salicaceae_regions/mummer_regions/concatenated/iqtree_results/partition.txt.root_pruned.nwk \
  -m ../data/salicaceae/meta.short.tsv \
  -c orgName,family,genus \
  -t Salicaceae \
  --key_column accession \
  -o template.json

# Then make a template json file with panman nwk (identical topology as iqtree but have internal node labels anddistorted branch lengths)
newick_to_taxonium -i ../data/salicaceae/salicaceae_panMAT/panmans/salicaceae_concatenated.panman.nwk \
  -m ../data/salicaceae/meta.short.tsv \
  -c orgName,family,genus \
  -t Salicaceae \
  --key_column accession \
  -o template.panman.json

# Use iqtree's branch lengths with panman's internal node labels
python3 merge_taxonium.py --primary template.json --secondary template.panman.json -o template.merged.json
```
...to be continued.

### AluYxx (AlueYe5)tree

After transferring the split nhmmscan results from phoenix to silverbullet, I can now merge them and reorganize them
into a tabular format that are easier to parse for downstream analysis.

```bash
wdir=/scratch1/alan/lab_notebook/tes/alu_ucsc/by_family/aluyxx

cd $wdir
# Download the split nhmmscan results from phoenix.
tar -xvzf hmm_out_split.tar.gz
cat hmm_out_split/* | awk '{$1=$1; print}' OFS='\t' | cut -f 1-15 > hmm_out.tblout
# Keep the top 3 hits for each sequence and format the output so that each query is a single row.
awk -F'\t' -v OFS='\t' '
  function flush() {
    if (prev == "") return
    out = prev
    for (i = 1; i <= 3; i++) out = out "\t" (i <= n ? hits[i] : "")
    print out
  }
  $3 != prev {
    flush()
    prev = $3
    n = 0
  }
  { hits[++n] = $1 "\t" $13 "\t" $14 }
  END { flush() }
' hmm_out.tblout > hmm_out.top3.tsv
# Remove the split directory to clean up space
rm -r hmm_out_split/
```

As a sanity check, I'm going to confirm that the in-hs1-AluYe5-cluster sequences that are mislabeled align well to their top
Alu profile hits and align better than to AluYe5.

```bash
qry=$(grep '^hs1::AluYe5' hmm_out.top3.tsv  | awk '$2 != "AluYe5#SINE/Alu"' | head -n1 | cut -f 1)
seq=$(seqkit grep -n -p $qry primate_aluyxx.len_filtered.fa)
nhmmscan --cpu 32 ../../../primate_alus/hmm/primate_alu.ad.curated.hmm <(echo "$seq")  | grep -A 30 ">> AluSx1#SINE/Alu"
nhmmscan --cpu 32 ../../../primate_alus/hmm/primate_alu.ad.curated.hmm <(echo "$seq")  | grep -A 30 ">> AluYe5#SINE/Alu"
```
Yep they align well to their top hits and align better than to AluYe5. Not shown here cuz too long. 


Now get the AluYe5 related copies from hominoid assemblies:

1. Copies labeled as AluYe5 from UCSC genome browser, `aluye5_labeled`
2. Copies whose top hit is AluYe5, `aluye5_tophit`
3. Copies labeled as AluYe5 from UCSC genome browser && whose top hit is AluYe5, `aluye5_labeled_tophit`

```bash
echo "hg38 hs1 panTro6 panPan3 gorGor6 ponAbe3 nomLeu3" | tr ' ' '\n' | grep -f - hmm_out.top3.tsv | grep "::AluYe5::" > hmm_out.labeled_aluye5.tsv
echo "hg38 hs1 panTro6 panPan3 gorGor6 ponAbe3 nomLeu3" | tr ' ' '\n' | grep -f - hmm_out.top3.tsv | awk '$2 == "AluYe5#SINE/Alu"' > hmm_out.tophit_aluye5.tsv
comm -12 <(sort hmm_out.labeled_aluye5.tsv) <(sort hmm_out.tophit_aluye5.tsv) > hmm_out.labeled_tophit_aluye5.tsv
```

```console
$ wc -l hmm_out.labeled_tophit_aluye5.tsv hmm_out.labeled_aluye5.tsv hmm_out.tophit_aluye5.tsv
   4793 hmm_out.labeled_tophit_aluye5.tsv
   9014 hmm_out.labeled_aluye5.tsv
   9362 hmm_out.tophit_aluye5.tsv
```

Hmmm, only 4793 copies are both labeled as AluYe5 and have AluYe5 as their top hit...

```console
$ awk -F'\t' ' $2 == "AluYe5#SINE/Alu" { print $5 }' hmm_out.top3.tsv | sort | uniq -c | sort -rn | head
   8198 AluYe6#SINE/Alu
    780 AluYf1#SINE/Alu
    698 AluY#SINE/Alu
    322 AluYk3#SINE/Alu
     52 AluYm1#SINE/Alu
     37 PapAnu-1.6#SINE/Alu
     36 AluYh3#SINE/Alu
     31 AluYi6#SINE/Alu
     25 AluSc8#SINE/Alu
     12 AluSc5#SINE/Alu

$ awk -F'\t' ' $2 == "AluYe6#SINE/Alu" { print $5 }' hmm_out.top3.tsv | sort | uniq -c | sort -rn | head
   2912 AluYe5#SINE/Alu
    273 AluYf1#SINE/Alu
     66 AluYk3#SINE/Alu
     58 AluY#SINE/Alu
     28 PapAnu-1.6#SINE/Alu
     16 AluYc#SINE/Alu
     12 AluSc8#SINE/Alu
     10 AluYi6#SINE/Alu
      8 AluYm1#SINE/Alu
      8 AluSx3#SINE/Alu
```

Seems like AluYe5 and AluYe6 are very similar to each other, which make sense since AluYe6 has all 5 of AluYe5's
diagnostic mutations plus 1 additional one, so they are almost indistinguishable. 

**I will actually just built an AluYeX tree with all AluYeX sequences (which is just AluYe5 and AluYe6 for now).**

Get the AluYe5 and AluYe6 sequences from the hominoid assemblies. Then trim the poly-A tails as before (refer to
[4/24/2026](#4242026)). Intermediate files are removed to clear up space.

```bash
wdir=/scratch1/alan/lab_notebook/tes/alu_ucsc/by_family/aluyex
mkdir -p $wdir && cd $wdir

# Also made a hmm_out.top5.tsv for some more information
echo "hg38 hs1 panTro6 panPan3 gorGor6 ponAbe3 nomLeu3" | tr ' ' '\n' | sed 's/^/^/g' | grep -f - hmm_out.top5.tsv | \
  awk '$2 == "AluYe5#SINE/Alu" || $2 == "AluYe6#SINE/Alu"' > aluyex.tsv

# Get fasta
cut -f 1 aluyex.tsv  | seqkit grep -n -f - ../aluyxx/primate_aluyxx.len_filtered.fa > aluyex.fa

# Followed steps in [4/24/2026](#4242026).
seqkit stats aluyex.fa aluyex.polyA_trimmed.fa 2> /dev/null 
sed -e 's/:/_/g' -e 's/([+-])//'  aluyex.polyA_trimmed.fa > aluyex.polyA_trimmed.renamed.fa 
```

```console
file                     format  type  num_seqs    sum_len  min_len  avg_len  max_len
aluyex.fa                FASTA   DNA     12,546  3,774,059      250    300.8      391
aluyex.polyA_trimmed.fa  FASTA   DNA     12,546  3,497,145      210    278.7      362
```

Then align using mafft-auto, mafft-retree, and dipper-twilight.

```bash
mafft --auto aluyex.polyA_trimmed.renamed.fa --thread 32 > aluyex.polyA_trimmed.auto.aln
mafft --retree 2 --maxiterate 1000 --thread 32 aluyex.polyA_trimmed.renamed.fa > aluyex.polyA_trimmed.retree.aln
dipper_cpu -i r -I aluyex.polyA_trimmed.renamed.fa -O dipper_results/aluyex.polyA_trimmed.dipper.guide.nwk
twilight -t dipper_results/aluyex.polyA_trimmed.dipper.guide.nwk -i aluyex.polyA_trimmed.renamed.fa -o aluyex.polyA_trimmed.twilight.aln -v -w --check -r 0.999  --cpu-only -C 32

# Then trim
trimal -in aluyex.polyA_trimmed.auto.aln -out aluyex.polyA_trimmed.auto.trimmed.aln -automated1
trimal -in aluyex.polyA_trimmed.retree.aln -out aluyex.polyA_trimmed.retree.trimmed.aln -automated1
trimal -in aluyex.polyA_trimmed.twilight.aln -out aluyex.polyA_trimmed.twilight.trimmed.aln -automated1

# Run dipper
for file in *aln; do
  prefix=$(basename $file .aln)
  seqkit seq -uv -w 0 $file > tmptmptmp.fa && mv tmptmptmp.fa $file
  dipper_cpu -i m -I $file -O dipper_results/${prefix}.nwk -d 4 --threads 32
done
```

## 4/30/2026

### Salicaceae tree

Continuing from [4/29/2026](#4292026), gonna finish writing the script that generates a Taxonium/Nextclade json file for
viewing the results.

Run `gen_taxonium.py` or `plot_results.py` to generate the Taxonium/Nextclade json file for viewing the results.

```bash
cd /scratch1/alan/lab_notebook/panmama/salicaceae/results

# Zihao's data
python3 plot_results.py -t salicaceae_concatenated.panman.nwk \
  -l from_zihao/KapK.cpDNA.mgsr.assignedReadsLCANode.out \
  -n from_zihao/KapK.cpDNA.mgsr.assignedReads.out \
  -m ../data/salicaceae/meta.tsv \
  --color-node-labels genus \
  -g 0.5 \
  -o from_zihao/KapK.cpDNA

python3 plot_results.py -t salicaceae_concatenated.panman.nwk \
  -l from_zihao/KapK.cpDNA.mgsr.assignedReadsLCANode.out \
  -n from_zihao/KapK.cpDNA.mgsr.assignedReads.out \
  -m ../data/salicaceae/meta.tsv \
  --color-by lca_subtree_count \
  --size-by lca_count \
  --color-node-labels genus \
  -g 0.7 \
  -o from_zihao/KapK.cpDNA.lca_subtree

python3 gen_taxonium.py -t salicaceae_concatenated.panman.nwk \
  -l from_zihao/KapK.cpDNA.mgsr.assignedReadsLCANode.out \
  -n from_zihao/KapK.cpDNA.mgsr.assignedReads.out \
  -e template.merged.jsonl \
  -g 0.7 \
  > from_zihao/KapK.cpDNA.jsonl

# Bianca's data.. Not gonna do gen_taxonium for now.
for file in from_bianca/ERR10493*assignedReadsLCANode.out; do
  prefix=$(basename $file .assignedReadsLCANode.out)
  num_reads=$(awk '{sum += $2} END {print sum}' $file)
  python3 plot_results.py -t salicaceae_concatenated.panman.nwk \
  -l $file \
  -n from_bianca/${prefix}.assignedReads.out \
  -m ../data/salicaceae/meta.tsv \
  --color-node-labels genus \
  -g 0.7 \
  --color-by lca_subtree_count \
  --size-by lca_count \
  -o from_bianca/${prefix}.${num_reads}r.lca_subtree
done

for file in from_bianca/ERR10493*assignedReadsLCANode.out; do
  prefix=$(basename $file .assignedReadsLCANode.out)
  num_reads=$(awk '{sum += $2} END {print sum}' $file)
  python3 plot_results.py -t salicaceae_concatenated.panman.nwk \
  -l $file \
  -n from_bianca/${prefix}.assignedReads.out \
  -m ../data/salicaceae/meta.tsv \
  --color-node-labels genus \
  -g 0.7 \
  -o from_bianca/${prefix}.${num_reads}r
done
```

### AluYeX tree

Make a metadata file for this the AluYeX tree. Code written by Claude with my supervision.

```bash
python3 write_metadata.py \
  --fasta aluyex.polyA_trimmed.renamed.fa \
  --assemblies ../../ucsc_primate_assemblies.selected.tsv \
  --nhmm aluyex.tsv \
  --out metadata.tsv
```

The tree looks okish. I'm gonna try to download all the annotations for Primate Alus from Dfam.

```bash
wdir=/scratch1/alan/lab_notebook/tes/alu_dfam
mkdir -p $wdir && cd $wdir

# Get the primate alus using famdb.py with accession id
# Then get the accession ids
grep 'name=' primate_alu.ad.curated.acc.fa | cut -f 1 -d ' ' | sed -e 's/^>//g' -e 's/#SINE\/Alu//g' | cut -f 1 -d '.' > primate_alu.ad.curated.acc.txt
for acc in $(cat primate_alu.ad.curated.acc.txt); do
  bash download_annotations.sh $acc annotations
  sleep 5
done

# Download the assemblies
wget https://www.dfam.org/releases/ref-genomes/hg38/dfamseq
```

## 5/1/2026

### Alu tree

I will process the non-redundant human Alu copies from dfam for now.

```bash
wdir=/scratch1/alan/lab_notebook/tes/alu_dfam/hg38_annotations
cd $wdir && mkdir nr_fasta

# Get the fasta files
for file in nr/*; do prefix=$(basename $file .tsv); bedtools getfasta -fi ../ref_genomes/hg38.fa -bed <(awk 'BEGIN{OFS="\t"} !/^#/ {
  if ($10 < $11) { start = $10 - 1; end = $11 }
  else           { start = $11 - 1; end = $10 }
  print $1, start, end, $3"|"$2, $4, $9
}' $file) -s -name > nr_fasta/${prefix}.fa; done

# Get sequences > 2/3 of the length of the consensus sequence
seqkit fx2tab -n -l ../primate_alu.ad.curated.acc.fa  | grep 'name=' | awk '{print $1"\t"$2"\t"$NF}' | sed 's/name=//g'| \
  sort -k3,3 -gr > ../consensus_len.tsv

for file in nr_fasta/*hg38_nrph-true.fa; do
  acc=$(basename $file | cut -f 1 -d '_')
  len=$(grep $acc ../consensus_len.tsv | cut -f 3)
  min_len=$((len * 2 / 3))
  seqkit seq -uvg -w 0 -m $min_len $file > ${file%.fa}.${min_len}bp.fa
done


# remove reads with ambiguous bases and remove substring duplicates
cd nr_fasta
for file in *; do
  sga preprocess --pe-mode 0 $file > ${file%.fa}.preprocessed.fa 
  sga rmdup -t 4 ${file%.fa}.preprocessed.fa 
done

# remove poly-a tails
for file in *.rmdup.fa; do 
  acc=$(echo $file | cut -f 1 -d '_')
  consensus=../hard_modified_consensus/${acc}.fasta
  split_dir=$(basename $file .fa)_split
  merged_trimmed_file=${file%.fa}.trimmed.fa
  seqkit split -i --by-id-prefix '' -O $split_dir $file
  shopt -s extglob
  export consensus
  printf '%s\n' "$split_dir"/!(*trimmed.fa) | parallel -j 32 'f={}; p=${f%.fa}; cat "$consensus" "$f" | mafft --auto --thread 1 - 2>/dev/null | python3 ../../../scripts/trim_polya.py - -o "${p}.trimmed.fa" && rm "$f"'
  find "$split_dir" -maxdepth 1 -name '*.trimmed.fa' -print0 | xargs -0 cat > "$merged_trimmed_file"
  rm -r $split_dir
done

cat nr_fasta_atrimmed/DF00* > cleaned.merged.fa
```

!! NOTE: Consensus sequences for trimming poly-a tails were hard-modified to substitute non-A in the middle of the
poly-A tails for ease of code writing:

- `DF003893998.2#SINE/Alu name=ASR @Simiiformes`
- `DF000000073.4#SINE/Alu name=BC200 @Primates`

After moving things around the clearing up intermediate files, I now have the the final set of non-redundant human Alu 
copies with at least 2/3 of the consensus sequence length, no ambiguous bases, no substring duplicate, and poly-A
tails removed.

## 5/3/2026

We have quite a lot of sequences.

```console
$ seqkit stats cleaned.merged.name_clean.fa
file                          format  type   num_seqs      sum_len  min_len  avg_len  max_len
cleaned.merged.name_clean.fa  FASTA   DNA   1,013,046  268,789,363       34    265.3      484
```

Try to use dipper and twilight to buil a tree

```bash
# build a guide tree first using dipper
dipper_cpu -i r -I cleaned.merged.name_clean.fa -O cleaned.merged.name_clean.guide.nwk --threads 32 > dipper.log 2> dipper.err
# Then use twilight to align the sequences
```

While that's running, I'm gonna also try and build an aluyx tree, which has 126,235 sequences.

```bash
# build a guide tree first using dipper
dipper_cpu -i r -I cleaned.aluyx.name_clean.fa -O cleaned.aluyx.name_clean.guide.nwk --threads 32 > dipper.aluy.log 2> dipper.aluy.err
# Then use twilight to align the sequences
twilight -t cleaned.aluyx.name_clean.guide.nwk -i cleaned.aluyx.name_clean.fa -o cleaned.aluyx.name_clean.twilight.aln -v -w --check -r 0.999  --cpu-only -C 32
```

Maybe worth trying mafft (auto mode) to align the aluyx sequences as well. Gonna do it on phoenix.

```bash
input_fa=/private/groups/corbettlab/alan/lab_notebook/tes/alu_famdb/human_aluyx/cleaned.aluyx.name_clean.fa
sbatch -J aluyx_human_dbfam run_mafft_retree.sh $input_fa ${$input_fa%.fa}.aln 2 0
```

Trim the alignments
```bash
trimal -in cleaned.aluyx.name_clean.twilight.aln -out cleaned.aluyx.name_clean.twilight.trimmed.aln -automated1
trimal -in cleaned.aluyx.name_clean.mafftauto.aln -out cleaned.aluyx.name_clean.mafftauto.trimmed.aln -automated1
```

Then build a tree using dipper

```bash
# after alignment, run dipper 
sbatch -J dipper_mafft_aluyx_human run_dipper.sh \
  /private/groups/corbettlab/alan/lab_notebook/tes/alu_famdb/human_aluyx/cleaned.aluyx.name_clean.mafftauto.trimmed.aln \
  /private/groups/corbettlab/alan/lab_notebook/tes/alu_famdb/human_aluyx/cleaned.aluyx.name_clean.mafftauto.nwk

sbatch -J dipper_twilight_aluyx_human run_dipper.sh \
  /private/groups/corbettlab/alan/lab_notebook/tes/alu_famdb/human_aluyx/cleaned.aluyx.name_clean.twilight.trimmed.aln \
  /private/groups/corbettlab/alan/lab_notebook/tes/alu_famdb/human_aluyx/cleaned.aluyx.name_clean.twilight.nwk
```

The tree (from twilight.gappyout_trimmed.nwk) looks like this:

![](tes/alu_dfam/hg38_annotations/taxonium.twilight.gappyout_trimmed.5_4_2026.png)

It has some super long and weird branches. They might be misannotated or super divergent sequences. I'm gonna
investigate and toss out some bad quality sequences.

Make a meta data file first...
```bash
wdir=/scratch1/alan/lab_notebook/tes/alu_dfam/hg38_annotations
cd $wdir 
(for file in nr/DF00*AluY*; do prefix=$(basename $file .tsv); paste <(seqkit seq -n nr_fasta/${prefix}.fa | sed -e 's/(.)//g' -e 's/:/_/g') <(tail -n+2 $file) ; done) > aluyx.meta.tsv
((echo "sample  target  model_acc       model_name      bitscore        evalue  hmm_start       hmm_end hmm_len strand  aln_start aln_end   env_start  env_end    target_length" | awk -v OFS='\t' '{$1=$1; print}') && cat aluyx.meta.tsv) > tmp
mv tmp aluyx.meta.tsv
```

```bash
python3 plot_score_vs_distance.py cleaned.aluyx.name_clean.twilight.gappyout_trimmed.jsonl aluyx.meta.tsv 
python filter_alu.py aluyx.meta.tsv \
  -o aluyx.meta.filtered.tsv \
  --hard-cutoff 1.0 \
  --summary aluyx.meta.summary.tsv
```

In figures below, `x_coords` is the x-coordinate of the sequence in the Taxonium tree, which reflects its divergence
from the root. So it seems bit score is inversely correlated to divergence from the root, which is not surprising...

`Alignment length` is the length of the hmm alignment, which should, in large, have a positive, linearly correlation with
the bit score. I might be able to throw out sequences with low bit scores relative to their hmm alignment length.

![](tes/alu_dfam/hg38_annotations/dotplots.png)

## 5/4/2026 and 5/5/2026

And here is the density of bitscore / alignment length:

![](tes/alu_dfam/hg38_annotations/score_density_by_family.png)

And bitscore vs alignment length scatter plot (cut off bitscore / alignment length < 1.0):

![](tes/alu_dfam/hg38_annotations/bitscore_vs_alnlen.png)


Before I rerun the alignment and rebuild the tree all over again with the filtered set, I will first just prune the 
discarded sequences from the original tree and see what it looks like.

```bash
# discarded sequences...
comm -13 <(cut -f1 aluyx.meta.filtered.tsv | tail -n+2 | sort) <(cut -f1 aluyx.meta.tsv | tail -n+2 | sort) \
  > discarded.low_0nout8.txt 

# that are in the final tree
comm -12 <(seqkit seq -n  cleaned.aluyx.name_clean.fa | sort) discarded.low_0nout8.txt \
  > discarded.low_0nout8.intree.txt

# prune the tree
nw_prune -f cleaned.aluyx.name_clean.twilight.gappyout_trimmed.nwk discarded.low_bitscore.intree.txt \
  > cleaned.aluyx.name_clean.twilight.gappyout_trimmed.pruned.nwk &
nw_prune -f cleaned.aluyx.name_clean.mafftauto.gappyout_trimmed.nwk discarded.low_bitscore.intree.txt \
  > cleaned.aluyx.name_clean.mafftauto.gappyout_trimmed.pruned.nwk &
```

Now the tree (from mafftauto.gappyout_trimmed.pruned.nwk) looks like this:

![](tes/alu_dfam/hg38_annotations/taxonium.mafftauto.gappyout_trimmed.5_5_2026.png)

It still has a few long branches relative to the rest of the tree but I think it's good enough for now.

Now I will actually add all the consensus sequences as well as an outgroup to the filtered fasta, then realign and
rebuild the tree.

```bash
# clean up the old aln and tree files
mv cleaned.aluyx.name_clean.*aln cleaned.aluyx.name_clean.*jsonl cleaned.aluyx.name_clean.*nwk archive.aluyx.unfiltered_bitscore/
tar -czvf archive.aluyx.unfiltered_bitscore.tar.gz archive.aluyx.unfiltered_bitscore/

# grab the bitscore filtered fasta
seqkit grep -v -n -f discarded.low_0nout8.txt cleaned.aluyx.name_clean.fa > cleaned.aluyx.name_clean.bitscore_filtered.fa

# add the consensus sequences and an outgroup (an AluSc8 consensus). Poly-A trimmed of course.
cat hard_modified_consensus/* | seqkit grep -n -r -p '.*name=AluY.*' | seqkit seq -w 0 | \
  awk '/^>/{print; next} {sub(/A+$/, ""); print}' | sed -e 's/#SINE\/Alu//g' -e 's/ name=/_/g' | \
  seqkit seq -uvi -w 60 >> cleaned.aluyx.name_clean.bitscore_filtered.fa

seqkit seq -w 0 hard_modified_consensus/DF000000038.fasta | \
  awk '/^>/{print; next} {sub(/A+$/, ""); print}' | sed -e 's/#SINE\/Alu//g' -e 's/ name=/_/g' | \
  seqkit seq -uvi -w 60 >> cleaned.aluyx.name_clean.bitscore_filtered.fa


# On phoenix, run mafft --auto, mafft --retree 2, and dipper-twilight to get aln files
# On phoenix, run dipper on all aligned files
```

eghhhhhh... It seems like some of the annotations from dfam are also not correted..... which might explain super low
bitscore / alignment length ratios for some of the AluY subfamilies. See the score density of all families at
`tes/alu_dfam/hg38_annotations/all.score_density_by_family.png`...

I guess I will run my own `nhmmscan` on the annotated regions again.

```bash
# sync merged.nhmmerscan.tblout from phoenix to silverbullet
awk -v OFS='\t' '$1=$1 {print $0}' merged.nhmmerscan.tblout  > tmp && mv tmp merged.nhmmerscan.tblout
python3 filter_hmmscan.py merged.nhmmerscan.tblout merged.nhmmerscan
# sort by max-hit family name
sort -k16,16 merged.hmmerscan.hmm_max.tblout   > tmp && mv tmp merged.hmmerscan.hmm_max.tblout
```

## 5/6/2026

Reassign the copies to their max hit families and preprocess the fasta files for alignment and stuff.

```bash
wdir=/scratch1/alan/lab_notebook/tes/alu_dfam/hg38_annotations/hmm_max_hits
mkdir -p $wdir && cd $wdir
mv ../merged.nhmmerscan.* .

# rename the copies to match their max hit families
python3 scripts/rename_copy.py merged.nhmmerscan.hmm_max.tblout > merged.nhmmerscan.hmm_max.name_corrected.tblout

# split the file into multiple files by max hit family name
awk -F'\t' '{
  key = $1
  sub(/#.*/, "", key)
  print > (key ".tblout")
}' merged.nhmmerscan.hmm_max.name_corrected.tblout

# convert the hmmtblout files to bed files
for file in split_hmmout/*tblout; do
  python3 scripts/hmmtblout2bed.py $file > beds/$(basename $file .tblout).bed &
done
wait

# get fasta, this also flips the strand if the alignment is reversed
for file in beds/*bed; do
  prefix=$(basename $file .bed)
  bedtools getfasta -fi ../../ref_genomes/hg38.fa -bed $file -s -name > fastas/${prefix}.fa &
done
wait

# hmmalign copies to trim off poly-a tails
for file in fastas/*fa; do
  prefix=$(basename $file .fa)
  hmmalign --outformat afa -o ${prefix}.hmmaln \
    ../../../primate_alus/hmm/split/${prefix}.hmm \
    $file &
done
wait

# rerun nhmmscan (so reversed alignments are now in positive alignment coordinates)
for file in fastas/*fa; do
  prefix=$(basename $file .fa)
  q --cpu 2 --noali --tblout hmmscan_out/${prefix}.tblout ../../../primate_alus/hmm/split/${prefix}.hmm $file &
done
wait

# remove header lines and reformat spaces to tabs
for file in hmmscan_out/*; do
  sed -i '/^#/d' $file
  awk -v OFS='\t' '$1=$1 {print $0}' $file > tmp && mv tmp $file
done

# Keep only the max hit for each sequence
for file in hmmscan_out/*; do awk -F'\t' '
  $3 != prev {
    if (NR > 1) print best
    prev = $3
    best = $0
    max = $14
    next
  }
  $14 + 0 > max + 0 {
    best = $0
    max = $14
  }
  END { if (NR > 0) print best }
' $file > tmp && mv tmp $file
done

# get the position of the last non-poly-a base in consensus (1-based inclusive), i.e. lengh of the consensus sequence minus the length of the poly-a tail
(for file in hmmalign/*; do 
  prefix=$(basename "$file" .hmmaln)
  acc=$(awk -v p="$prefix" '$2 == p' ../../accession_to_name.tsv | cut -f1 -d '.')
  pos=$(seqkit seq -w 0 ../hard_modified_consensus/${acc}.fasta | sed 's/A*$//g' | seqkit stats | tail -n1 | awk '{print $5}')
  echo -e "${prefix}\t${acc}\t${pos}"
done) > nonpolya_last_base.tsv
```

## 5/8/2026

Now I will build the human AluYx panMAT map some simulated human reads onto it.

```bash
wdir=/scratch1/alan/lab_notebook/tes/alu_dfam/hg38_annotations
cd $wdir
mkdir -p panmat

# All the alignment and nwk files are copied from phoenix to alignment/ dir
# Build panman
for nwk in alignment/*nwk; do
  prefix=$(basename $nwk .nwk)
  alignment=alignment/${prefix}.upper.aln
  /scratch1/alan/panmap/build/bin/panmanUtils -M $alignment -N $nwk -o ${prefix} --threads 16 &
done
wait

# Then build index
for panman in panmat/*panman; do
  /scratch1/alan/panmap/build/bin/panmap $panman --index-mgsr ${panman%.panman}.idx &
done
wait
```

Simulate some hg38 reads using wgsim

```bash
# Simulate 100M paired reads of length 150bp -> ~10x coverage of the human genome. Random seed is 69, nice.
wgsim -1 150 -2 150 -N 100000000 -e 0.005 -r 0.001 -R 0.01 -X 0.1 -S 69 hg38.core.fa simulated.10x_R1.fastq simulated.10x_R2.fastq
```

Run panmap

```bash
/scratch1/alan/panmap/build/bin/panmap ../panmat/cleaned.aluyx.name_clean.bitscore_filtered.auto.panman \
  simulated.10x_R1.fastq simulated.10x_R2.fastq \
  -i ../panmat/cleaned.aluyx.name_clean.bitscore_filtered.auto.idx \
  --meta --filter-and-assign --discard 0.9  -t 16 \
  --output simulated.10x
```

## 5/10/2026

Run script to get the family distribution of the results

```bash
python3 ../get_family_distribution.py simulated.10x.assignedReadsLCANode.out  simulated.10x
```

Also it with the hs1 assembly

```bash
wdir=/scratch1/alan/lab_notebook/tes/alu_dfam/hg38_annotations/simulation/hs1; mkdir -p $wdir && cd $wdir

# Simulate 100M paired reads of length 150bp -> ~10x coverage of the human genome. Random seed is 69, nice.
wgsim -1 150 -2 150 -N 100000000 -e 0.005 -r 0.001 -R 0.01 -X 0.1 -S 69 hs1.fa simulated.10x_R1.fastq simulated.10x_R2.fastq

/scratch1/alan/panmap/build/bin/panmap ../panmat/cleaned.aluyx.name_clean.bitscore_filtered.auto.panman \
  simulated.10x_R1.fastq simulated.10x_R2.fastq \
  -i ../panmat/cleaned.aluyx.name_clean.bitscore_filtered.auto.idx \
  --meta --filter-and-assign --discard 0.9  -t 16 \
  --taxonomic-metadata ../hg38/meta.tsv --taxonomic-rank family \
  --output simulated.10x

python3 ../../get_family_distribution.py simulated.10x.mgsr.assignedReadsLCANode.out  simulated.10x
```

The family hit distribution for both the hg38 and hs1 assemblies are pretty similar to the true distributions in the hg38 genome.

```console
(echo -e "family\tperc(%)\thg38(%)\ths1(%)" && grep -v AluSc hg38/family_proportions.tsv | tail -n+2 |  while read family count perc; do hg38_perc=$(grep "$family\s" hg38/simulated.10x.family_distribution.tsv | cut -f 3); hs1_perc=$(grep "$family\s" hs1/simulated.10x.family_distribution.tsv | cut -f 3); echo -e "${family}\t${perc}\t${hg38_perc}\t${hs1_perc}"; done) | column -t

family   perc(%)    hg38(%)   hs1(%)
AluY     61.9227    63.27427  63.11049
AluYk3   14.2511    14.51000  14.41058
AluYm1   5.76422    5.65363   5.62706
AluYh3   3.49367    3.28037   3.31778
AluYf1   3.40791    3.45101   3.40657
AluYa5   2.65938    1.45155   1.43539
AluYc    1.87256    1.74491   1.72729
AluYb8   1.81594    1.00798   0.98499
AluYe5   0.980825   0.87819   0.86545
AluYk4   0.939194   0.92281   0.92472
AluYh7   0.888404   0.86860   0.88896
AluYi6   0.534541   0.46039   0.45302
AluYc3   0.397992   0.37334   0.37166
AluYg6   0.385502   0.28896   0.41987
AluYb9   0.203992   0.10433   0.11391
AluYe6   0.142378   0.13717   0.13962
AluYd8   0.137382   0.09318   0.09121
AluYk11  0.080764   0.06568   0.06650
AluYk12  0.0616138  0.02985   0.05585
AluYa8   0.0407983  0.02718   0.02270
AluYh9   0.0183176  0.01446   0.01587
```

## 5/11/2026

### Estimating AluY proportions of pantro6 using the hg38 AluY panmat

```bash
wdir=/scratch1/alan/lab_notebook/tes/alu_dfam/hg38_annotations/simulation/panTro6; mkdir -p $wdir && cd $wdir

# Simulate 100M paired reads of length 150bp -> ~10x coverage of the human genome. Random seed is 69, nice.
wgsim -1 150 -2 150 -N 100000000 -e 0.005 -r 0.001 -R 0.01 -X 0.1 -S 69 panTro6.core.fa simulated.10x_R1.fastq simulated.10x_R2.fastq

# Run panmap
/scratch1/alan/panmap/build/bin/panmap ../../panmat/cleaned.aluyx.name_clean.bitscore_filtered.auto.panman \
  simulated.10x_R1.fastq simulated.10x_R2.fastq \
  -i ../../panmat/cleaned.aluyx.name_clean.bitscore_filtered.auto.idx \
  --meta --filter-and-assign --discard 0.9  -t 16 \
  --taxonomic-metadata ../hg38/meta.tsv --taxonomic-rank family \
  --batch-size 100000 \
  --output simulated.10x

# Get the family distribution
python3 ../../get_family_distribution.py simulated.10x.mgsr.assignedReadsLCANode.out  simulated.10x

# Get the pantro AluY distributions from dfam annotations
zcat ../../../annotations/*DApanTro2_nrph-true.tsv.gz | cut -f3 | grep AluY | uniq -c | awk '{print $2"\t"$1}' > panTro6.aluy.counts.tsv
(echo -e "family\tcounts\tperc(%)" && awk 'BEGIN{FS=OFS="\t"} NR==FNR{sum+=$2; next} {printf "%s\t%d\t%.4f\n", $1, $2, ($2/sum)*100}' panTro6.aluy.counts.tsv panTro6.aluy.counts.tsv) > panTro6.aluy.proportions.tsv
```

Comparing the family proportions from the dfam annotations and the estimated proportions from the panmap results.

```console
(echo -e "family\tfamdb(%)\testimated(%)" && join -t $'\t' -a 1 -a 2 -e "0" -o '0,1.3,2.3'   <(tail -n +2 panTro6.aluy.proportions.tsv | sort -k1,1)   <(tail -n +2 simulated.10x.family_distribution.tsv | sort -k1,1)   | awk 'BEGIN{FS=OFS="\t"} {print}' | sort -k2,2 -k3,3 -gr) | column -t

family   famdb(%)  estimated(%)
AluY     61.7815   66.89192
AluYk3   12.9284   14.26758
AluYm1   6.1193    5.56688
AluYh3   5.1540    3.36265
AluYf1   3.6716    3.04981
AluYc    3.6282    1.45518
AluYh7   2.4261    0.79552
AluYk4   1.2719    0.87418
AluYe5   1.0675    0.79731
AluYc3   0.8213    0.46659
AluYh9   0.5666    0.01966
AluYe6   0.5636    0.17341
.        0         1.46770
AluYi6   0         0.38435
AluYa5   0         0.20022
AluYg6   0         0.06972
AluYb8   0         0.06614
AluYk11  0         0.04112
AluYk12  0         0.03754
AluYb9   0         0.00536
AluYd8   0         0.00358
AluYa8   0         0.00358
```

Some of the proportions are quite off... AluYi6 and AluYa5 were estimated to be present for a non-trivial proportion of
the reads but not in the dfam annotations. This might be due to incomplete representation of the AluYx copies or might 
simply be due to incomplete annotations of the AluYx copies in the dfam annotations. The same goes for cases like AluYh9
that has a very low estimated proportion but is present in the dfam annotations. Will investigate this further.


### HMM max hits processing

Pick up from the end of [5/6/2026](#562026)

Use the hmm alignment to trim off poly-As and positions outside of the hmm alignment.

```bash
# trim off poly-As and positions outside of the hmm alignment
for file in hmmalign/*; do 
  prefix=$(basename $file .hmmaln)
  python3 scripts/process_hmmalignment.py $file \
    inenv_polya_trimmed/${prefix}.fa \
    -p nonpolya_last_base.tsv \
    -a hmmscan_out/${prefix}.tblout &
done
wait

# only trim off poly-As
for file in hmmalign/*; do 
  prefix=$(basename $file .hmmaln)
  python3 scripts/process_hmmalignment.py $file \
    polya_trimmed/${prefix}.fa \
    -p nonpolya_last_base.tsv \
    -a hmmscan_out/${prefix}.tblout \
    -k &
done
wait
```

## 5/12/2026 and 5/13/2026

### Build a primate AluYx panMAT

```bash
wdir=/scratch1/alan/lab_notebook/tes/alu_dfam/primate_annotations/fastas; mkdir -p $wdir && cd $wdir

# Get the fasta files for the primate AluYx copies
for ref in ../dfam_references/*fa; do
  assembly=$(basename "$ref" .fa)
  mkdir -p "${assembly}"
  
  for annotation_file in $(find ../../ -name "*${assembly}*_nrph-true.tsv.gz"); do
    model_name=$(basename "$annotation_file" | cut -f1 -d '_')
    model_name=$(grep "^${model_name}" ../../accession_to_name.tsv | cut -f2)
    bedtools getfasta -fi "$ref" -bed <(zcat "$annotation_file" | awk 'BEGIN{OFS="\t"} !/^#/ {
      start = ($10 < $11) ? $10 - 1 : $11 - 1
      end = ($10 < $11) ? $11 : $10
      print $1, start, end, $3"|"$2, $4, $9
    }') -s -name > "${assembly}/${model_name}.fa" &
  done
  wait
done

# Merge the fasta files for each assembly
for dir in *; do
  cat ${dir}/*.fa > ${dir}.fa
done

# Add assembly name to the fasta headers
for dir in */; do
  assembly="${dir%/}"
  for file in "$dir"*.fa; do
    [ -f "$file" ] || continue
    sed -i "s/^>/>${assembly}|/" "$file"
  done
done

for f in *.fa; do
  assembly="${f%.fa}"
  sed -i "s/^>/>${assembly}|/" "$f"
done

# Group fastas by family. This is the unfiltered fastas grabbed directly from fastas dir.
wdir=/scratch1/alan/lab_notebook/tes/alu_dfam/primate_annotations/fastas_by_family; mkdir -p $wdir && cd $wdir
for family in $(l -1 ../fastas/*/* | awk -F'/' '{print $NF}' | sed 's/.fa//g' | sort -u); do
  cat ../fastas/*/${family}.fa > ${family}.fa
done

# Run hmmalign on unfiltered, family grouped fastas on phoenix for poly-A trimming and maybe some other filtering. Used run_hmmalign_2026-05-12.sh on phoenix.
wdir=/scratch1/alan/lab_notebook/tes/alu_dfam/primate_annotations
seqkit seq -w 0 ../../primate_alus/primate_alu.ad.curated.fa | sed 's/[Aa]*$//g' | seqkit fx2tab  -i -n -l | grep -v short | sed 's/#SINE\/Alu//g' > nonpolya_last_base.tsv

find hmmalign_out_by_family -mindepth 2 -type f -name "*hmmaln"| \
  parallel -j 8 'python3 ../../scripts/process_hmmalignment.py {} {.}.fa -p nonpolya_last_base.tsv -k'

mkdir fastas_by_family_polya_trimmed
for dir in hmmalign_out_by_family/*; do
  family=$(basename $dir _hmmalign_out)
  cat $dir/*fa > fastas_by_family_polya_trimmed/${family}.fa
done

# clean up intermediate files
rm -r hmmalign_out_by_family/*_hmmalign_out

# Make assembly-specific hmm profile for each assembly
wdir=/scratch1/alan/lab_notebook/tes/alu_dfam/primate_annotations/hmm_models; mkdir -p $wdir && cd $wdir
for dir in ../fastas/*; do
  [ -f "$dir" ] && continue
  assembly=$(basename "$dir")
  families=$(ls -1 $dir/*.fa | sed 's/\.fa//g' | awk -F '/' '{print $NF}')
  for family in $families; do
    cat ../../../primate_alus/hmm/split/${family}.hmm >> ${assembly}.hmm
  done
done

# Run nhmmscan on phoenix to get the max hmm hits using assembly specific hmm profiles. Used run_nhmmerscan_2026-05-12.sh on phoenix.
# Reformat and reassign the copies to their max hit families
wdir=/scratch1/alan/lab_notebook/tes/alu_dfam/primate_annotations/nhmmscan_out; cd $wdir
for file in *tblout; do
  awk -v OFS='\t' '$1=$1 {print $0}' $file > ${file%.tblout}.tmp && mv ${file%.tblout}.tmp $file
  awk -v OFS='\t' '$16 == "7SL" {$16="7SLRNA"} 1' $file > ${file%.tblout}.tmp && mv ${file%.tblout}.tmp $file
  python3 ../../../scripts/filter_hmmscan.py $file ${file%.tblout}
  python3 ../../../scripts/rename_copy.py ${file%.tblout}.hmm_max.tblout > ${file%.tblout}.hmm_max.renamed.tblout
done
```


Work on AluYx specific panMAT

```bash
wdir=/scratch1/alan/lab_notebook/tes/alu_dfam/primate_annotations/aluyx; mkdir -p $wdir && cd $wdir
cat ../fastas_by_family/AluY*fa > aluyx.unfiltered.fa
cat ../fastas_by_family_polya_trimmed/AluY*fa > aluyx.polya_trimmed.fa

# Check out the length distribution
seqkit fx2tab -n -l -i aluyx.unfiltered.fa  | cut -f2 | sort -g | uniq -c > len_stats.unfiltered.tsv &
seqkit fx2tab -n -l -i aluyx.polya_trimmed.fa  | cut -f2 | sort -g | uniq -c > len_stats.polya_trimmed.tsv &

python3 ../../../scripts/plot_len_distribution.py len_stats.unfiltered.tsv len_dist.unfiltered.png --bin-width 5 &
python3 ../../../scripts/plot_len_distribution.py len_stats.polya_trimmed.tsv len_dist.polya_trimmed.png --bin-width 5 &

# Length filter reads. Remove reads with ambiguous bases and remove substring duplicates
sga preprocess --pe-mode 0 aluyx.polya_trimmed.fa -m 100 > aluyx.polya_trimmed.preprocessed.fa
sga index -t 16 aluyx.polya_trimmed.preprocessed.fa
sga rmdup -t 16 aluyx.polya_trimmed.preprocessed.fa

# Rename headers to be tree friendly
seqkit seq -i -w 0 aluyx.polya_trimmed.preprocessed.rmdup.fa | sed 's/(.)//g' | tr ':' '_' | sed 's/#SINE\/Alu//g' > aluyx.polya_trimmed.preprocessed.rmdup.renamed_headers.fa
```

I also ran the same commands on AluSx and AluJx copies.

## 5/14/2026

### AluYx PanMAT related

After the job failed initially I tried different other options and better GPUs but nothing worked.

I talked with Sumit and changed the sketch-seed setting to 100 from the default (1000). The default number of seeds is 
too high for the relatively short Alu sequences, making most of the sequences look very similar to each other, resulting
in the clusters bigger than the backbone.

If that doesn't work I will try to increase the backbone size, which Sumit told me how to do.

Yay the tree built after ~8 hours! Wait.. TWILIGHT failed.. It seemed like the final tree built is missing exactly one
sample?

I'm going to rerun dipper with both the docker pulled version and the latest github version.

I also found that DIPPER fails when the job gets assigned to phoenix-02 but works on phoenix-05. I forced the job to be
assigned to phoenix-05 for now. Will report this issue to Sumit and maybe the cluster admins.

It seemed like the error arose actually because TWILIGHT wasn't parsing the floats correctly when reading in the tree.
I removed the branch length information from the nwk and it worked.

### Side project: make a better verterbate mito tree

I learned that, similar to chloroplast genomes, mitochondrial genomes should also be built by splitting the genome into
regions then aligning each region individually followed by concatenation and tree inference. I'm trying out this tool
called `mitos` to annotate the mitochondrial genomes and extract the regions.

```bash
wdir=/scratch1/alan/lab_notebook/vertebrate_new_panman; mkdir -p $wdir && cd $wdir
ls v_mtdna.dedupped.split/*.fa | parallel -j 32 --bar '
  name=$(basename {} .fa)
  mkdir -p annotations/$name
  runmitos \
    -i {} \
    -c 2 \
    -o annotations/$name \
    -r refseq89m \
    -R . \
    --noplots \
    > logs/$name.stdout 2> logs/$name.stderr
'
```

I will be aligning the 13 protein coding genes (PCGs) and the 2 rRNAs (12S and 16S) and using their alignments to build
the tree. Gonna exclude the tRNAs.

For future me the 13 PCGs are:

```console
nad1 — ND1, ~950 nt
nad2 — ND2, ~1040 nt
nad3 — ND3, ~345 nt
nad4 — ND4, ~1380 nt
nad4l — ND4L, ~295 nt
nad5 — ND5, ~1810 nt
nad6 — ND6, ~525 nt
cox1 — COI, ~1540 nt
cox2 — COII, ~680 nt
cox3 — COIII, ~785 nt
atp6 — ATP6, ~680 nt
atp8 — ATP8, ~200 nt
cob — CYTB, ~1140 nt
```

12 of the 13 are on the heavy strand. nad6 is the odd one out, on the light strand, and has weird base composition
because of it — might need to give it special treatment (separate partition, or RY-coding, or just drop it) if I see
weirdness in the tree later.
And the 2 rRNAs:

```console
rrnS — 12S, ~950 nt
rrnL — 16S, ~1550 nt
```

That gives me 15 partitions total, ~14kb of aligned sequence per taxon once everything's concatenated.

**For alignment:**

PCGs: use MACSE v2, codon-aware. Vertebrate mitochondrial genetic code = translation table 2:

```
AUA = Met (not Ile)
UGA = Trp (not stop)
AGA/AGG = stop (not Arg)
MACSE handles this if I pass --gc_def 2 (or whatever the flag is — check docs).
```

rRNAs: MAFFT --auto is fine, or L-INS-i for better accuracy. Structure-aware would be better but probably overkill.
Trim with trimAl -automated1 or -gappyout before concatenation.

### aeDNA related 

Gonna run some more samples Bianca and Zihao sent me.

1. cpDNA of Populus genus, Populus_besthit_diff_1.fastq. *Already ran earlier*.
2. FLB sediment Salix data, competitively mapped, Salix reads extracted
3. FLB sediment data, not mapped, just straight up clean fqs
4. FLB bone data, not mapped
5. KapK competitively mapped Populus reads, deduplicated, potential NUPTs and MTPTs removed.


First run the flb sediment salix reads

```bash
wdir=/scratch1/alan/lab_notebook/panmama/salicaceae/data/flb_sediment_salix
cd $wdir

# convert bam to fastqs
for file in *bam; do
  samtools fastq $file >> merged.primary.fq
  samtools fastq $file -F 0 >> merged.all.fq
done

# Run panmap
wdir=/scratch1/alan/lab_notebook/panmama/salicaceae/results/flb_sediment_salix; mkdir -p $wdir && cd $wdir

/scratch1/alan/panmap/build/bin/panmap ../../data/salicaceae/salicaceae_panMAT/panmans/salicaceae_concatenated.panman \
  ../../data/flb_sediment_salix/merged.all.fq \
  -i ../../data/salicaceae/salicaceae_panMAT/panmans/salicaceae_concatenated.idx \
  --meta --filter-and-assign \
  --discard 0.5  -t 8 \
  --output merged.all

# Plot results (all and primary, plain and lca_subtree)
python3 ../../plot_results.py -t ../../data/salicaceae/salicaceae_panMAT/panmans/salicaceae_concatenated.panman.nwk -l merged.all.mgsr.assignedReadsLCANode.out -n merged.all.mgsr.assignedReads.out -m ../../data/salicaceae/meta.tsv --color-node-labels genus -g 0.7 -o merged.all
python3 ../../plot_results.py -t ../../data/salicaceae/salicaceae_panMAT/panmans/salicaceae_concatenated.panman.nwk -l merged.all.mgsr.assignedReadsLCANode.out -n merged.all.mgsr.assignedReads.out -m ../../data/salicaceae/meta.tsv --color-node-labels genus --color-by lca_subtree_count --size-by lca_count -g 0.7 -o merged.all.lca_subtree
```

FLB sediment unmapped reads

```bash
wdir=/scratch1/alan/lab_notebook/panmama/salicaceae/results/flb_sediment_unmapped; mkdir -p $wdir && cd $wdir

data_path=/storage2/alan/flb_sediment_unmapped
/scratch1/alan/panmap/build/bin/panmap ../../data/salicaceae/salicaceae_panMAT/panmans/salicaceae_concatenated.panman \
  --batch-files-path ${data_path}/batch_files.tsv \
  -i ../../data/salicaceae/salicaceae_panMAT/panmans/salicaceae_concatenated.idx \
  --meta --filter-and-assign \
  --discard 0.6  -t 16 --batch-size 100000

python3 ../gather_assignedReads.py -i *clean.mgsr.assignedReadsLCANode.out -o assignedReadsLCANode.merged.out
python3 ../gather_assignedReads.py -i *clean.mgsr.assignedReads.out -o assignedReads.merged.out

python3 ../../plot_results.py -t ../../data/salicaceae/salicaceae_panMAT/panmans/salicaceae_concatenated.panman.nwk -l assignedReadsLCANode.merged.out -n assignedReads.merged.out -m ../../data/salicaceae/meta.tsv --color-node-labels genus -g 0.7 -o merged.all
python3 ../../plot_results.py -t ../../data/salicaceae/salicaceae_panMAT/panmans/salicaceae_concatenated.panman.nwk -l assignedReadsLCANode.merged.out -n assignedReads.merged.out -m ../../data/salicaceae/meta.tsv --color-node-labels genus --color-by lca_subtree_count --size-by lca_count -g 0.7 -o merged.all.lca_subtree

```

FLB bone data

```bash
wdir=/scratch1/alan/lab_notebook/panmama/salicaceae/results/flb_bone; mkdir -p $wdir && cd $wdir

data_path=/storage2/alan/flb_bone
/scratch1/alan/panmap/build/bin/panmap ../../data/v_mtdna/v_mtdna.new.panman \
  --batch-files-path ${data_path}/batch_files.tsv \
  -i ../../data/v_mtdna/v_mtdna.pmai \
  --meta --filter-and-assign \
  --taxonomic-metadata ../../data/v_mtdna/v_mtdna.meta.tsv \
  --taxonomic-rank Family \
  --discard 0.6  -t 16 
```

KapK competitively mapped Populus reads

```bash
wdir=/scratch1/alan/lab_notebook/panmama/salicaceae/results/kapk_populus_dedup_remap; mkdir -p $wdir && cd $wdir

for idx_param in k11_s6_l1 k13_s7_l1 k15_s8_l1; do
  outprefix=B.KapK.Populus.cpDNA.dedup.remap.${idx_param}
  /scratch1/alan/panmap/build/bin/panmap ../../data/salicaceae/salicaceae_panMAT/panmans/salicaceae_concatenated.panman \
    ../../data/kapk_populus_dedup_remap/B.KapK.Populus.cpDNA.dedup.remap.fastq.gz \
    -i ../../data/salicaceae/salicaceae_panMAT/panmans/salicaceae_concatenated.${idx_param}.idx \
    --meta --filter-and-assign \
    --discard 0.5  -t 8 \
    --output ${outprefix}

done
```

## 5/15/2026

### aeDNA related

Let's first take a look at the FLB bones.

```bash
wdir=/scratch1/alan/lab_notebook/panmama/salicaceae/results/flb_bone; mkdir -p $wdir && cd $wdir

mkdir -p analysis/taxon_reads
for file in *fastq; do
  prefix=$(basename $file .mgsr.assignedReads.fastq)
  python3 ../../../sediment_arctic_4mil/collect_reads.py $prefix all analysis/taxon_reads/${prefix} &
done
wait

mkdir -p analysis/aligned_family
for file in *fastq; do
  prefix=$(basename $file .mgsr.assignedReads.fastq)
  python3 ../../../sediment_arctic_4mil/align_family.py \
    ../../../sediment_arctic_4mil/data/fastas \
    analysis/taxon_reads/${prefix} \
    ../../../sediment_arctic_4mil/data/v_mtdna.meta.tsv \
    --out_dir analysis/aligned_family/${prefix} &
done
wait

mkdir -p analysis/damage
for sample in analysis/aligned_family/*; do
  outdir=analysis/damage/$(basename ${sample})
  mkdir -p $outdir
  for bam in ${sample}/*.merged.bam; do
    prefix=$(basename $bam .merged.bam)
    out_prefix=${outdir}/${prefix}
    python3 ../../../sediment_arctic_4mil/plot_damage.py compute --bam ${bam} --out-subs ${out_prefix}.subs.tsv --out-plot ${out_prefix}.plot.png --tax-name ${prefix} && \
    python3 ../../../sediment_arctic_4mil/plot_damage.py score --in-subs ${out_prefix}.subs.tsv --out-tsv ${out_prefix}.score.tsv --tax-name ${prefix} &
  done
  wait
done
```

The FLB bone reads don't seem to be qc'd. Gonna qc it now.

Now it's qc'd. Looking at the results, it doesn't seem like Panmap identified any significant hits.

### Side project: make a better verterbate mito tree

Just got all the annotations. Get some stats

```bash
wdir=/scratch1/alan/lab_notebook/vertebrate_new_panman; mkdir -p $wdir && cd $wdir

bash summarize_mitos.sh annotations/ regions.tsv
```

## 5/18/2026

Maybe the nhmmscan and dfam annotation discrepancies was because dfam annotations were scanned across the entire genome,
which might have genomic contexts that could make it different from nhmmscan on only the sequence itself. Gonna try
to give it ~200 bp of flanks

```bash
wdir=/scratch1/alan/lab_notebook/tes/alu_dfam/primate_annotations_with_flanks/fastas; mkdir -p $wdir && cd $wdir

for ref in ../../primate_annotations/dfam_references/*fa; do
  assembly=$(basename "$ref" .fa)
  mkdir -p "${assembly}"

  genome_file="../../primate_annotations/dfam_references/${assembly}.genomes.txt"
  for annotation_file in $(find ../../annotations -maxdepth 1 -name "*${assembly}*_nrph-true.tsv.gz"); do
    model_name=$(basename "$annotation_file" | cut -f1 -d '_')
    model_name=$(grep "^${model_name}" ../../accession_to_name.tsv | cut -f2)
    bedtools getfasta -fi "$ref" -bed <(zcat "$annotation_file" | awk 'BEGIN{OFS="\t"} !/^#/ {
      start = ($10 < $11) ? $10 - 1 : $11 - 1
      end = ($10 < $11) ? $11 : $10
      print $1, start, end, $3"|"$2, $4, $9
    }' | bedtools slop -i - -g "$genome_file" -b 200 -s) -s -name > "${assembly}/${model_name}.fa" &
  done
  wait
done
```

Then follow steps from [5/12/2026 and 5/13/2026](#5122026-and-5132026)

## 5/19/2026

I'm going to run dfamscan on whole genomes to see if I get the same results as dfam annotations. This is done on
Phoenix.

Meanwhile I will work on my lab meeting.

### Content for lab meeting

Make a stacked bar plot of the Alu copies from dfam annotations.

```bash
wdir=/scratch1/alan/lab_notebook/tes/alu_dfam/primate_annotations; mkdir -p $wdir && cd $wdir

python3 ../../scripts/plot_alu_counts.py fastas/all.stats.tsv -o alu_totals.png
```

Make a bar plot comparing the counts of UCSC annotated AluYe5 copies between assemblies.

```bash
wdir=/scratch1/alan/lab_notebook/tes/alu_ucsc/by_family/aluye5; mkdir -p $wdir && cd $wdir
python3 plot_counts.py counts_by_assembly.tsv -o aluye5_counts.png -t "AluYe5 counts by assembly"
```

Make two pie charts comparing the nhmmscan results between hs1-incluster and hs1-outcluster aluye5 copies.

```bash
wdir=/scratch1/alan/lab_notebook/tes/alu_ucsc/by_family/aluye5/hs1_cluster; mkdir -p $wdir && cd $wdir
python3 plot_alu_pies.py hs1_incluster.hmm_tophit_counts.tsv hs1_outcluster.hmm_tophit_counts.tsv   --label1 "In-cluster" --label2 "Out-of-cluster" -o alu_family_pies.png
```

Simualted run results

```
family   hg38_famdb(%)  hg38_famdb_filtered(%)  estm_hg38(%)  estm_hs1(%) | family   pantro(%)  estimated(%)
AluY     55.3862        61.9227                 63.27427      63.11049    |  AluY     61.7815    66.89192
AluYk3   13.2246        14.2511                 14.51000      14.41058    |  AluYk3   12.9284    14.26758
AluYm1   5.5811         5.76422                 5.65363       5.62706     |  AluYm1   6.1193     5.56688
AluYh3   4.7205         3.49367                 3.28037       3.31778     |  AluYh3   5.1540     3.36265
AluYf1   3.6093         3.40791                 3.45101       3.40657     |  AluYf1   3.6716     3.04981
AluYa5   3.4467         2.65938                 1.45155       1.43539     |  AluYc    3.6282     1.45518
AluYc    3.1104         1.87256                 1.74491       1.72729     |  AluYh7   2.4261     0.79552
AluYb8   2.0589         1.81594                 1.00798       0.98499     |  AluYk4   1.2719     0.87418
AluYe5   1.1516         0.980825                0.87819       0.86545     |  AluYe5   1.0675     0.79731
AluYk4   1.0569         0.939194                0.92281       0.92472     |  AluYc3   0.8213     0.46659
AluYh7   2.0974         0.888404                0.86860       0.88896     |  AluYh9   0.5666     0.01966
AluYi6   0.7652         0.534541                0.46039       0.45302     |  AluYe6   0.5636     0.17341
AluYc3   0.6479         0.397992                0.37334       0.37166     |  AluYi6   0          0.38435
AluYg6   0.6843         0.385502                0.28896       0.41987     |  AluYa5   0          0.20022
AluYb9   0.3857         0.203992                0.10433       0.11391     |  AluYg6   0          0.06972
AluYe6   0.4372         0.142378                0.13717       0.13962     |  AluYb8   0          0.06614
AluYd8   0.3438         0.137382                0.09318       0.09121     |  AluYk11  0          0.04112
AluYk11  0.2443         0.080764                0.06568       0.06650     |  AluYk12  0          0.03754
AluYk12  0.2128         0.0616138               0.02985       0.05585     |  AluYb9   0          0.00536
AluYa8   0.4269         0.0407983               0.02718       0.02270     |  AluYd8   0          0.00358
AluYh9   0.4084         0.0183176               0.01446       0.01587     |  AluYa8   0          0.00358
Ambigu   .              .                       1.36214       1.55051     |  Ambigu   0          1.46770
```

## 5/20/2026

I gave lab meeting and sent over the Salicaeceae and the bone DNA placements to Zihao and Bianca.

## 5/21/2026

Yay, Zihao will include the KapK Salicaceae placements in his paper, and Bianca will include me on the ellesmere paper
author list (not my results but she appreciated my work)!

## 5/22/2026

I think I will just throw away all the copies that are assigned different families by hmmerscan across the whole genome
and local sequence only. i.e. only keep the ones that could be confidently assigned to one family regardless of its
flanking region (keep when the top hit and the annotation are in the same family but not necessarily the same subfamily).
I think we'd lose about ~10% of the copies but that should be fine.

```bash
wdir=/scratch1/alan/lab_notebook/tes/alu_dfam/primate_annotations/nhmmscan_out; cd $wdir
for file in $(find . -maxdepth 1 -name "*tblout" ! -name "*hmm_max*tblout"); do
  file=$(basename $file)
  awk '{
    f1 = $1; sub(/#.*/, "", f1)
    f3 = $3; sub(/\|.*/, "", f3)
    if (f1 == f3) print
  }' $file > ${file%.tblout}.matched.tblout &
done
wait &

# then reprocess the AluY copies first
wdir=/scratch1/alan/lab_notebook/tes/alu_dfam/primate_annotations/aluyx; cd $wdir
cat ../fastas_by_family_polya_trimmed/AluY*fa > aluyx.polya_trimmed.fa
sed -i 's|#SINE/Alu||g' aluyx.polya_trimmed.fa

for file in ../nhmmscan_out/*.hmm_max.tblout; do
  assembly=$(basename $file .hmm_max.tblout)
  awk '$1 ~ /AluY/ && $3 ~ /AluY/' $file | cut -f3 | sed "s/^/${assembly}\|/g" > ${assembly}.aluy.matched.ids.txt &
done
wait &

cat *.aluy.matched.ids.txt | sed 's|#SINE/Alu||g' > aluyx.matched.ids.txt && rm *.aluy.matched.ids.txt
seqkit grep -f aluyx.matched.ids.txt aluyx.polya_trimmed.fa  | seqkit seq -w 0 > aluyx.polya_trimmed.matched.fa

# len filter
seqkit seq -m 250 -M 350 -w 0 aluyx.polya_trimmed.matched.fa  > aluyx.polya_trimmed.matched.len_filtered.fa

# # remove exact substring duplicates
# sga preprocess --pe-mode 0 aluyx.polya_trimmed.matched.fa > aluyx.polya_trimmed.matched.preprocessed.fa
# sga index -t 16 aluyx.polya_trimmed.matched.preprocessed.fa
# sga rmdup -t 16 aluyx.polya_trimmed.matched.preprocessed.fa


# cluster the matched AluY copies to subsample the sequences
cd-hit-est -i aluyx.polya_trimmed.matched.fa -o aluy_rep.fasta -c 0.95 -n 10 -M 16000 -T 16 -d 0
cd-hit-est -i aluyx.polya_trimmed.matched.len_filtered.fa -o aluy_rep.len_filtered.fasta -c 0.95 -n 10 -M 16000 -T 16 -d 0
cd-hit-est -i aluyx.polya_trimmed.matched.len_filtered.fa -o aluy_rep.len_filtered.99.fasta -c 0.99 -n 10 -M 128000 -T 32 -d 0
```

While the clustering is running, I'm gonna run DIPPER on the unclsutered, len filtered AluY copies 

```bash
for m in 3 1; do
  for s in 100 200; do
    for k in 15 13; do
      sbatch -J dipper_m${m}_s${s}_k${k} run_dipper_gpu_from_fa.sh \
        /private/groups/corbettlab/alan/lab_notebook/tes/alu_famdb/primate_alus/aluyx/aluyx.polya_trimmed.matched.len_filtered.renamed.fa \
        /private/groups/corbettlab/alan/lab_notebook/tes/alu_famdb/primate_alus/aluyx/aluyx.polya_trimmed.matched.len_filtered.renamed.m${m}_s${s}_k${k}.nwk \
        ${m} ${s} ${k}
    done
  done
done
```

## 5/26/2026

All of the dipper runs failed but the `m=3, s=100, k=15` setting. Going to trim it off now.

I also don't think I will use the clustering results for anything. Setting the sequence identity threshold low (0.95)
loses signals and setting it high (0.99) create too many clusters and still risks losing some signals. THe tree itself
automatically clusters the sequences and compresses similar sequences using "evolutionary compression" anyway...

Today I will start working on the Github tutorial for peqg.

## 5/27/2026

I just finished the Github tutorial for peqg.

## 5/28/2026

Finished trimming the first aluyx alignment. Gonna run dipper on it now.

```bash
for m in 3 1; do
  sbatch -J dipper_m${m} run_dipper_gpu_from_aln.sh \
    /private/groups/corbettlab/alan/lab_notebook/tes/alu_famdb/primate_alus/aluyx/aluyx.polya_trimmed.matched.len_filtered.renamed.trimmed.aln \
    /private/groups/corbettlab/alan/lab_notebook/tes/alu_famdb/primate_alus/aluyx/aluyx.polya_trimmed.matched.len_filtered.renamed.trimmed.m${m}.nwk \
    ${m}
done
```

And build a PanMAN from the first alignment to see how it looks

```bash

/scratch1/alan/panmap/build/bin/panmanUtils \
  -M aluyx.polya_trimmed.matched.len_filtered.renamed.aln \
  -N aluyx.polya_trimmed.matched.len_filtered.renamed.m3_s100_k15.nwk \
  --threads 32 \
  -o aluyx
```

## Land plant-wide chloroplast PanMAN

### 5/29/2026

I just had my committee meeting. It went well. The committee members were happy with my progress. They warn me that
modern sequencing capability allows us to fairly easily and cheaply assemble genomes to discover TEs and curate
reference pangenome pangenomes could also be used by mapping reads directly to the pangenome to identify TE copies.

**They recommended me to fully consider building a land plant-wide chloroplast PanMAN**

Searching all the complete land plant chloroplast genomes in RefSeq. In NCBI (nucleotide), search use this query:

```text
Embryophyta[ORGN] AND chloroplast[filter] AND ("complete genome"[Title] OR "complete plastid genome"[Title] OR "complete chloroplast genome"[Title]) AND refseq[filter] 
```

There are 12449 results as of 11:36 AM on 5/29/2026.

Get the list of accession numbers and metadata.

```bash
# Get accession numbers, taxids, names, etc.
esearch -db nuccore -query 'Embryophyta[ORGN] AND chloroplast[filter] AND ("complete genome"[Title] OR "complete plastid genome"[Title] OR "complete chloroplast genome"[Title]) AND refseq[filter]' \
  | esummary \
  | xtract -pattern DocumentSummary \
      -element AccessionVersion,TaxId,Organism,Slen,CreateDate,UpdateDate,Title \
  > plastome_metadata.tsv

# Get taxonomy
cut -f2 plastome_metadata.tsv | sort -u > taxids.txt

epost -db taxonomy -input taxids.txt \
  | efetch -format xml \
  | xtract -pattern Taxon -tab "\t" \
      -KING "(-)" -PHYL "(-)" -CLSS "(-)" -ORDR "(-)" -FMLY "(-)" -GNUS "(-)" \
      -block "*/Taxon" \
        -if Rank -equals "kingdom" -KING ScientificName \
        -if Rank -equals "phylum"  -PHYL ScientificName \
        -if Rank -equals "class"   -CLSS ScientificName \
        -if Rank -equals "order"   -ORDR ScientificName \
        -if Rank -equals "family"  -FMLY ScientificName \
        -if Rank -equals "genus"   -GNUS ScientificName \
      -group Taxon \
        -element TaxId,ScientificName,"&KING","&PHYL","&CLSS","&ORDR","&FMLY","&GNUS" \
  > taxonomy.tsv

python3 merge_meta.py plastome_metadata.tsv taxonomy.tsv > plastome_metadata_with_taxonomy.tsv
```

Special cases are handled as follows:

| Scientific name type | Example | Parsing strategy |
| --------------------- | ------- | ----------- |
| Genus hybrid cultivar | Aeonium hybrid cultivar | set species name to `sp.` and hybrid type to `hybrid_cultivar` |
| \[Genus\] species ... | \[Polygonum\] chinense var. procumbens | set genus to `Genus` without the square brackets |
| aff. samples | Euphorbia aff. neococcinea Luke s.n. | set species name to `sp.`|
| cf. samples | Urera cf. cordifolia Ur15 | set species name to `sp.`|
| subsp. samples | Brassica rapa subsp. rapa | set subspecies name |
| var. samples | Capsicum baccatum var. pendulum | set varietas name |
| sp. samples | Tetrastigma sp. Wen 12461 | set species name to `sp.` and uniden_id to identifiers after `sp.`, e.g. `Wen 12461` |

We have chloroplast genomes from a total of 16 classes.

```console
$ cut -f 6 plastome_metadata_with_taxonomy.tsv | sort -u
Andreaeopsida
Anthocerotopsida
Bryopsida
Cycadopsida
Gnetopsida
Haplomitriopsida
Jungermanniopsida
Leiosporocerotopsida
Lycopodiopsida
Magnoliopsida
Marchantiopsida
Pinopsida
Polypodiopsida
Polytrichopsida
Sphagnopsida
Tetraphidopsida
```

Seems like we are missing:

- Ginkgoopsida (Ginkos): the missing gymnosperm class
  - **Found many complete Ginko chloroplast genomes but none are in refseq, using `AB684440.1`**
- Andreaeobryopsida and Oedipodiopsida (two moss classes): very tiny classes with very few species
  - **Andreaeobryopsida has a single species: Andreaeobryum *macrosporum*. No complete chloroplast genome is available.**
  - **Oedipodiopsida has a single species: Oedipodium *griffithianum*. Using `OZ374950.1`.**
- Takakiopsida (two species in genus Takakia): enigmatic moss lineage.
  - Two species in Akakipsida: Takakia *ceratophylla* and Takakia *lepidozioides*
  - **Found complete chloroplast genomes for *lepidozioides*: `AP014702.1`**

### 6/1/2026

I just filled in some assemblies from missing classes. All the ones added are not in RefSeq but are complete chloroplast 
genomes.

```bash
cat plastome_metadata.tsv addons/*/*meta.tsv  > plastome_metadata.with_addon.tsv
cat taxonomy.tsv addons/*/*taxonomy.tsv  > taxonomy.with_addon.tsv
python3 merge_meta.py plastome_metadata.with_addon.tsv taxonomy.with_addon.tsv > plastome_metadata_with_taxonomy.tsv
```

Now download all the assemblies. Just gonna download the entire refseq plastid archive then add the missing assemblies.

```bash
for i in 1 2 3; do
  wget https://ftp.ncbi.nlm.nih.gov/refseq/release/plastid/plastid.${i}.1.genomic.fna.gz
done

for file in *genomic.fna.gz; do
  cut -f1 plastome_metadata.with_addon.tsv | seqkit grep -f - $file >> plastomes.fa
done

cat addons/*/*fa >> plastomes.fa 
rm *genomic.fna.gz
```

Looks like we got everything.

```console
$ seqkit stats plastomes.fa
file          format  type  num_seqs        sum_len  min_len    avg_len  max_len
plastomes.fa  FASTA   DNA     12,452  1,910,342,521   15,553  153,416.5  345,184

$ wc -l plastome_metadata.with_addon.tsv
12452 plastome_metadata.with_addon.tsv
```

Plot genome length distribution.

```bash
python3 scripts/plot_genome_len_dist.py misc_data/genome_len.tsv 
```

![genome length distribution](embryophyta_cp/figures/genome_length_distribution.png)

345,184 is very long for a chloroplast genome. And it belongs to `NC_066227.1`.

```console
$ seqkit fx2tab plastomes.fa -nl | awk '{print $1,$NF}' |  sort -k2,2 -gr | head -n 5
NC_066227.1 345184
NC_031206.1 242575
NC_084419.1 234657
NC_088443.1 232302
NC_056348.1 232020
```

According to this [paper](https://www.biorxiv.org/content/10.1101/2025.10.06.680833v3), it's misassembled. The two other
assemblies on [NCBI](https://www.ncbi.nlm.nih.gov/nuccore/?term=Magnolia+patungensis%5BORGN%5D+AND+chloroplast%5Bfilter%5D+AND+(%22complete+genome%22%5BTitle%5D+OR+%22complete+plastid+genome%22%5BTitle%5D+OR+%22complete+chloroplast+genome%22%5BTitle%5D)) that are in RefSeq have 160,120 bp and 160,139 bp.

**Perhaps, I need to do some more qc.**

Nonetheless, I'm gonna self-align the chloroplast genomes using nucmer to identify the quadripartite regions.

```bash
parallel -j 64 'prefix=$(basename {1} .fa); bash ~/tools/misc/mummer_align_self_cpg.sh {1} {1} regions/${prefix}/${prefix}' ::: fastas/*.fa
```

I will also get all the Embryophyta chloroplast genomes from NCBI without the refseq filter to check for weird refseq
genomes. Follow commands from 5/29/2026 without the `refseq[filter]`.

### 6/4/2026

I'm gonna follow [paper](https://www.biorxiv.org/content/10.1101/2025.10.06.680833v3) to do some quality control.

First let's take a look at the 16S rRNA and 23s rRNA genes, or RF00177 (bacterial SSU/16S rRNA) and RF02541 (bacterial
LSU/23S rRNA)

```bash
wdir=/scratch1/alan/lab_notebook/embryophyta_cp/rrnas; mkdir -p $wdir; cd $wdir

mkdir rfam_models
wget https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/fasta_files/RF00177.fa.gz -O rfam_models/RF00177.fa.gz
wget https://rfam.org/family/RF00177/cm -O rfam_models/RF00177.cm
wget https://rfam.org/family/RF02541/cm -O rfam_models/RF02541.cm

cat rfam_models/*cm > rfam_models/plastid_rrna.cm
cmpress rfam_models/plastid_rrna.cm

mkdir rrna_hits
parallel -j 64 '
  chunk={1};
  base=$(basename $chunk .fa);
  cmsearch \
    --rfam \
    --cpu 1 \
    -E 0.001 \
    --tblout rrna_hits/${base}.tbl \
    --noali \
    rfam_models/plastid_rrna.cm \
    $chunk > rrna_hits/${base}.out
' ::: ../fastas/*.fa
```

Consolidate the results and create a summary.

```bash
head -n2 rrna_hits/AB684440.1.tbl > rrna_hits_all.tbl
for f in rrna_hits/*.tbl; do
  grep -v '^#' $f >> rrna_hits_all.tbl
done

awk '!/^#/ {print $1"\t"$3"\t"$16}' rrna_hits_all.tbl \
  | sort -k1,1 -k2,2 \
  | awk '
      {
        key=$1"\t"$2;
        count[key]++;
        if (best[key]=="" || $3+0 < best[key]+0) best[key]=$3;
      }
      END {
        for (k in count) print k"\t"count[k]"\t"best[k];
      }' \
  | sort -k1,1 -k2,2 > rrna_per_genome.tsv

cut -f1 ../plastome_metadata.tsv | sort -u > all_accessions.txt

awk '$2=="SSU_rRNA_bacteria"' rrna_per_genome.tsv | cut -f1 | sort -u > has_16s.txt
awk '$2=="LSU_rRNA_bacteria"' rrna_per_genome.tsv | cut -f1 | sort -u > has_23s.txt

comm -23 all_accessions.txt has_16s.txt > missing_16s.txt
comm -23 all_accessions.txt has_23s.txt > missing_23s.txt
```

Only one genome is missing the 16S and 23S rRNA genes. `NC_037503.1` is *Asarum minus* in the Asaraceae family. There
are 5 other Asarum species in our dataset, all of which have the 16S and 23S rRNA genes. `NC_037503.1` is also
significantly shorter than the other Asarum species (~15 kb vs ~190 kb), an obvious sign of misassembly/miss-annotation.

<details>

<summary>See details</summary>

```console
$ tail -n+1 missing_*
==> missing_16s.txt <==
NC_037503.1

==> missing_23s.txt <==
NC_037503.1

$ awk '$9 == "Asarum"' ../plastome_metadata.with_taxonomy.tsv 
NC_086579.1     3113364 190179  Viridiplantae   Streptophyta    Magnoliopsida   Piperales       Asaraceae       Asarum  Asarum chungbuensis     None    Asarum  chungbuensis    None    None    None    None    None    None    None    None    Asarum chungbuensis chloroplast, complete genome
NC_077489.1     447319  192892  Viridiplantae   Streptophyta    Magnoliopsida   Piperales       Asaraceae       Asarum  Asarum splendens        None    Asarum  splendens       None    None    None    None    None    None    None    None    Asarum splendens chloroplast, complete genome
NC_058740.1     366670  193105  Viridiplantae   Streptophyta    Magnoliopsida   Piperales       Asaraceae       Asarum  Asarum sieboldii f. maculatum   None    Asarum  sieboldii       None    None    None    None    None    None    None    None    Asarum maculatum voucher NIBRVP0000640516, 2006.4.1. chloroplast, complete genome
NC_058739.1     366669  193163  Viridiplantae   Streptophyta    Magnoliopsida   Piperales       Asaraceae       Asarum  Asarum misandrum        None    Asarum  misandrum       None    None    None    None    None    None    None    None    Asarum misandrum voucher NIBRVP0000640514, 2007.4.1. chloroplast, complete genome
NC_037503.1     76132   15553   Viridiplantae   Streptophyta    Magnoliopsida   Piperales       Asaraceae       Asarum  Asarum minus    None    Asarum  minus   None    None    None    None    None    None    None    None    Asarum minus chloroplast, complete genome
NC_037190.1     76098   193356  Viridiplantae   Streptophyta    Magnoliopsida   Piperales       Asaraceae       Asarum  Asarum sieboldii        None    Asarum  sieboldii       None    None    None    None    None    None    None    None    Asarum sieboldii voucher NIBR-VP0000640510 chloroplast, complete genome

$ awk '$9 == "Asarum"' ../plastome_metadata.with_taxonomy.tsv | cut -f1 | grep -f - rrna_per_genome.tsv
NC_037190.1     LSU_rRNA_bacteria       2       0
NC_037190.1     SSU_rRNA_bacteria       2       0
NC_058739.1     LSU_rRNA_bacteria       2       0
NC_058739.1     SSU_rRNA_bacteria       2       0
NC_058740.1     LSU_rRNA_bacteria       2       0
NC_058740.1     SSU_rRNA_bacteria       2       0
NC_077489.1     LSU_rRNA_bacteria       2       0
NC_077489.1     SSU_rRNA_bacteria       2       0
NC_086579.1     LSU_rRNA_bacteria       2       0
NC_086579.1     SSU_rRNA_bacteria       2       0

$ awk '$9 == "Asarum"' ../plastome_metadata.with_taxonomy.tsv | cut -f1 | grep -f - ../misc_data/genome_len.tsv 
NC_037190.1     193356
NC_058739.1     193163
NC_058740.1     193105
NC_077489.1     192892
NC_086579.1     190179
NC_037503.1     15553
```

</details>
<br/>

Good idea to check if genomes within the same family have consistent genome sizes. Plot the genome length distribution
by family, skipping families with less than 2 genomes and genome length range less than 10,000 bp.

```bash
python3 scripts/plot_len_by_taxon.py plastome_metadata.with_taxonomy.tsv -o figures/len_dist_by_taxon.png --box-width 0.3 --xtick-length 6.0
```

![genome length distribution by family](embryophyta_cp/figures/len_dist_by_taxon.png)

## Population Populus PanMAT

### 6/1/2026

Working on the population Populus PanMAT that Zihao sent me. This to do population levle analysis for genomes of 
P. *trichocarpa* and P. *balsamifera*.

```bash
wdir=/scratch1/alan/lab_notebook/panmama/salicaceae/data/For_Alan_population_level_panmat/cp; cd $wdir
```

According to Zihao, one IR was removed and the regions were aligned and concatenated using 100 Ns.

Gonna remove the N's using a custom script (made with Claude with close supervision).

```bash
python3 split_regions.py PopPt_31_conca_raw_Ns.fasta 100
mkdir -p split
mv region_1.fasta split/LSC.fasta
mv region_2.fasta split/IR.fasta
mv region_3.fasta split/SSC.fasta
```

```console
$ seqkit seq -w 0 PopPt_31_conca_raw_Ns.fasta | tr -d '-' | seqkit fx2tab -nl
KJ664926        129809
MW376847        128879
NC_037417       151181
SRR11622733     129521
SRR12235434     129424
SRR13324505     129463
SRR13324525     129291
SRR13324536     129161
SRR1569613      129546
SRR1759777      129429
SRR25364191     129408
SRR25364192     129379
SRR25364205     129630
SRR25364217     129314
SRR25364243     129430
SRR25364246     129364
SRR25364261     129410
SRR25364275     129491
SRR25364283     129288
SRR25364330     129320
SRR25364397     129444
SRR25364471     129591
SRR25364489     129550
SRR25364510     129333
SRR25364569     129428
SRR25364585     129226
SRR25364663     129339
SRR25364685     129395
SRR25364729     129478
SRR3045866      129030
SRR35248031     129750
```

One assemnly `NC_037417` seems weirdly long. after removing gaps. 

```console
$ seqkit grep -p NC_037417 PopPt_31_conca_raw_Ns.fasta | seqkit seq -w 0 | tr -d '-' | sed  's/NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN//g' | seqkit stats
file  format  type  num_seqs  sum_len  min_len  avg_len  max_len
-     FASTA   DNA          1  149,981  149,981  149,981  149,981
```

The full genome length is 158,551 bp but its sequence length after removing an IR region is 149,981 bp? `NC_037417` is
Populus *pruinosa* and likely an outgroup. For some reason the outgroup is not split up correctly...

Gonna ask Zihao to check if the outgroup is split up and concatednated correctly.

### 6/2/2026

Zihao just sent me the fixed version. Gonna check it now.

`NC_037417` looks good now.

```console
$ seqkit seq -w 0 PopPt_31_conca_raw_Ns_fixed.fasta | tr -d '-' | seqkit fx2tab -nl | head
KJ664926        129809
MW376847        128879
NC_037417       130626
SRR11622733     129521
SRR12235434     129424
SRR13324505     129463
SRR13324525     129291
SRR13324536     129161
SRR1569613      129546
SRR1759777      129429
```

Now remove N's and split the regions.

```bash
wdir=/scratch1/alan/lab_notebook/panmama/salicaceae/data/For_Alan_population_level_panmat/cp; cd $wdir

python3 split_regions.py PopPt_31_conca_raw_Ns.fasta 100
mkdir -p split
mv region_1.fasta split/LSC.fasta
mv region_2.fasta split/IR.fasta
mv region_3.fasta split/SSC.fasta
```

Regions also have expected lengths with no weirdness.

```console
$ (echo -e 'file\tnum_seq\tmin_len\tmax_len\tavg_len' && for file in split/*; do stats=$(seqkit seq -w 0 $file | tr -d '-' | seqkit stats | tail -n 1 | awk -v OFS='\t' '{print $4,$6,$8,$7}'); echo -e "${file}\t${stats}"  ; done) | column -t
file             num_seqs  min_len  max_len  avg_len
split/IR.fasta   31        27,034   28,125   27,654.5
split/LSC.fasta  31        84,641   85,761   85,015.2
split/SSC.fasta  31        16,318   16,659   16,573.7
```

I just confirmed that the regions should be concatenated in the order of LSC, IR, SSC. (NOT SSC, IR, LSC).

<details>
<summary>See alignment</summary>

`NC_037417.partial.fa` is concated from LSC, IR, SSC. It's obviously the correct one compared to
`NC_037417.partial2.fa`. The mismatch was the artifact of me arbitrarily substituing an `M` to an `A` in the original 
fasta file for `cpstools`.

```console
$ seqkit concat split/LSC.fasta split/IR.fasta split/SSC.fasta | seqkit grep -p NC_037417 | seqkit seq -w 0  | tr -d '-' > NC_037417.partial.fa
$ minimap2 -x asm5 -a --MD NC_037417.fasta  NC_037417.partial.fa  2> /dev/null | cut -f 10 --complement | grep -v '^@'  | column -t 
NC_037417  0  NC_037417.1  1  60  130426M  *  0  0  *  NM:i:1  ms:i:130424  AS:i:130424  nn:i:1  tp:A:P  cm:i:13022  s1:i:130378  s2:i:28109  de:f:0.0000  MD:Z:8987A121438  rl:i:0
```

`NC_037417.partial2.fa` is concated from SSC, IR, LSC.

```console
$ seqkit concat split/SSC.fasta split/IR.fasta split/LSC.fasta| seqkit grep -p NC_037417 | seqkit seq -w 0  | tr -d '-' > NC_037417.partial2.fa
$ minimap2 -x asm5 -a --MD NC_037417.fasta  NC_037417.partial2.fa  2> /dev/null | cut -f 10 --complement | grep -v '^@'  | column -t 
NC_037417  0     NC_037417.1  1       60  44665S85761M        *  0  0  *  NM:i:1  ms:i:85759  AS:i:85759  nn:i:1  tp:A:P  cm:i:8577  s1:i:85720  s2:i:0      de:f:0.0000  SA:Z:NC_037417.1,85760,+,16538S28129M85759S,1,0;NC_037417.1,113887,+,16540M113886S,60,0;  MD:Z:8987A76773  rl:i:0
NC_037417  2048  NC_037417.1  85760   1   16538H28129M85759H  *  0  0  *  NM:i:0  ms:i:28129  AS:i:28129  nn:i:0  tp:A:P  cm:i:2782  s1:i:28114  s2:i:28109  de:f:0       SA:Z:NC_037417.1,1,+,44665S85761M,60,1;NC_037417.1,113887,+,16540M113886S,60,0;           MD:Z:28129       rl:i:0
NC_037417  272   NC_037417.1  130427  0   85761S28125M16540S  *  0  0  *  NM:i:0  ms:i:28125  AS:i:28125  nn:i:0  tp:A:S  cm:i:2781  s1:i:28109  de:f:0      MD:Z:28125   rl:i:0                                                                                                     
NC_037417  2048  NC_037417.1  113887  60  16540M113886H       *  0  0  *  NM:i:0  ms:i:16540  AS:i:16540  nn:i:0  tp:A:P  cm:i:1658  s1:i:16519  s2:i:0      de:f:0       SA:Z:NC_037417.1,1,+,44665S85761M,60,1;NC_037417.1,85760,+,16538S28129M85759S,1,0;        MD:Z:16540       rl:i:0
```

</details>

Remove the 100-N linker from MSA

```bash
seqkit seq -w 0 PopPt_31_conca_raw_Ns_fixed.fasta \
  | sed 's/NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN//g' \
  > PopPt_31_conca_raw_Ns_removed.fasta
```

<details>
<summary>Sanity check</summary>

```console
$ seqkit concat split/LSC.fasta split/IR.fasta split/SSC.fasta > PopPt_31_regions_concat.fasta
$ seqkit sum PopPt_31_conca_raw_Ns_removed.fasta PopPt_31_regions_concat.fasta
seqkit.v0.1_DLS_k0_5aaecdb09efb77c02506dc2f83b39afa     PopPt_31_conca_raw_Ns_removed.fasta
seqkit.v0.1_DLS_k0_5aaecdb09efb77c02506dc2f83b39afa     PopPt_31_regions_concat.fasta
```

</details>

Now build a PanMAN and run panmap on the data

```bash
/scratch1/alan/panmap/build/bin/panmanUtils -M PopPt_31_conca_raw_Ns_removed.fasta -N PopPt_31_ultrametric.nwk -o PopPt_31

/scratch1/alan/panmap/build/bin/panmap PopPt_31.panman --index-mgsr PopPt_31.panman.idx -k 15 -s 8 -l 1

mkdir results
for read in aeDNA/Clade1_besthit_diff_1.ChrPt.remap.dedup.fastq.gz aeDNA/Populus_besthit_diff_1.ChrPt.remap.fastq.gz; do
  /scratch1/alan/panmap/build/bin/panmap PopPt_31.panman \
    $read \
    -i PopPt_31.panman.idx \
    --meta --filter-and-assign \
    --discard 0.5  -t 8 \
    --output results/$(basename $read .fastq.gz)
done

for file in results/*mgsr.assignedReadsLCANode.out; do
  prefix=$(basename $file .assignedReadsLCANode.out)
  python3 plot_results.py -t PopPt_31.panman.nwk \
    -l $file \
    -n results/${prefix}.assignedReads.out \
    -m pop.meta.tsv \
    --color-node-labels common_name \
    -g 0.7 \
    --color-by lca_subtree_count \
    --size-by lca_count \
    -o results/${prefix}.subtree &

  python3 plot_results.py -t PopPt_31.panman.nwk \
    -l $file \
    -n results/${prefix}.assignedReads.out \
    -m pop.meta.tsv \
    --color-node-labels common_name \
    -g 0.7 \
    -o results/${prefix}.lca &
  
  python3 plot_results.py -t PopPt_31.panman.nwk \
    -l $file \
    -n results/${prefix}.assignedReads.out \
    -m pop.meta.tsv \
    --color-node-labels common_name \
    -g 0.7 \
    --color-by fractionated_read_count \
    --size-by lca_count \
    -o results/${prefix}.fractionated &
  wait
done

```

### 6/3/2026

Ooops, there's also mito data to run for Zihao in the same data directory that he sent me.

```bash
wdir=/scratch1/alan/lab_notebook/panmama/salicaceae/data/For_Alan_population_level_panmat/mt
cd $wdir
```

The fasta file has the expected number of 100-N linked bases for all the entries.

```console
$seqkit seq -w 0 PopMt_32_conca.raw.Ns.fasta | grep -on 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN' | cut -f1 -d ':' | sort | uniq -c | awk '{print $1}' | sort | uniq -c 
     32 593
```

Remove the 100-N linked bases. Then build a PanMAN and run panmap on the data

```bash
seqkit seq -w 0 PopMt_32_conca.raw.Ns.fasta \
  | sed 's/NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN//g' \
  > PopMt_32_conca.raw.Ns_removed.fasta

/scratch1/alan/panmap/build/bin/panmanUtils -M PopMt_32_conca.raw.Ns_removed.fasta -N PopMt_32_taxon_ultrametric.nwk -o PopMt_32

/scratch1/alan/panmap/build/bin/panmap PopMt_32.panman --index-mgsr PopMt_32.panman.idx -k 15 -s 8 -l 1

mkdir results
for read in aeDNA/Clade1_besthit_diff_1.ChrMt.remap.dedup.fastq.gz  aeDNA/Populus_besthit_diff_1.ChrMt.remap.dedup.fastq.gz; do
  /scratch1/alan/panmap/build/bin/panmap PopMt_32.panman \
    $read \
    -i PopMt_32.panman.idx \
    --meta --filter-and-assign \
    --discard 0.5  -t 8 \
    --output results/$(basename $read .fastq.gz)
done

for file in results/*mgsr.assignedReadsLCANode.out; do
  prefix=$(basename $file .assignedReadsLCANode.out)

  python3 plot_results.py -t PopMt_32.panman.nwk \
    -l $file \
    -n results/${prefix}.assignedReads.out \
    -g 0.7 \
    --color-by lca_subtree_count \
    --size-by lca_count \
    -o results/${prefix}.subtree &

  python3 plot_results.py -t PopMt_32.panman.nwk \
    -l $file \
    -n results/${prefix}.assignedReads.out \
    -g 0.7 \
    -o results/${prefix}.lca &

  python3 plot_results.py -t PopMt_32.panman.nwk \
    -l $file \
    -n results/${prefix}.assignedReads.out \
    -g 0.7 \
    --color-by fractionated_read_count \
    --size-by lca_count \
    -o results/${prefix}.fractionated &

  wait
done 

```

### 6/4/2026

It seems like Zihao might have wanted the read counts to be split between tied genomes. I will implement this and rerun
the samples... Actually, I can just use the output file to compute it.. Added in the code blocks above.


