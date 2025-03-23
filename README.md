# Lab notebook

This notebook will track the progress of my work in the lab.

## 1/30/2025

### panMANIA

~~To build the panMAN using Dockerfile, see [panmania/build_panman/](panmania/build_panman/) directory.~~
To build the panMAN see notes on [2/2/2025](#222025).

PanMAN built using Dockerfile is broken. Currently using Summit's docker image.

## 1/31/2025

### panMANIA

Today I want to start by examining the alignments of sequences on a panMAN. I will write scripts that would compare the global alignments of sequences from the panMAN to the alignments of sequences using mafft. I expect to see that the mafft alignments should be better than from the panMAN because panMAN alignments are from MSA of all the sequences while mafft only compares two sequences at a time. Big differences, however, indicate particular problematic alignments, likely on the block level.

Inside evaluate_alignments, `parse_newick_to_pairs.py` is a script that would parse the newick tree and output the parent-child pairs in depth first order. The sequences of each parent-child pair will be compared.

Following commands will generate the pairs, run the alignment comparison, and plot the alignment differences.

```
python parse_newick_to_pairs.py data/hiv/hiv_panman.newick > data/hiv/hiv20000_pairs.tsv

sbatch evaluate_alignments/evaluate_alignments_hiv.sh

python evaluate_alignments/plot_alignment_diff.py evaluate_alignments/out/hiv20000_alignment_differences.tsv evaluate_alignments/data/hiv/hiv20000_pairs.tsv hiv evaluate_alignments/out/hiv20000_alignment_differences.png
```

For results, see [panmania/evaluate_alignments](panmania/evaluate_alignments).

## 2/1/2025

### panMANIA

I realized that there are newer versions of panmans availble. I'm going to build the latest version of panman and look at the alignments comparison again for newer panmans.

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

Edit the CMakeLists.txt by adding the following lines between `TARGET_LINK_LIBRARIES(...)` and `target_include_directories(...)`
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

I have started writing the specific aim section of my thesis proposal. For aim 1, other than figures for accuracy assessment, I think I can also have figure that shows the correlation between the pseudo-chaining score and the alignment score.

### panMAMA

When I'm in the progress of formalizing the method for calling consensus, I noticed a pretty significant bug in `heuristic_v6.py` that gives inconsistent `position_info_by_haplotype` for the same haplotype. `position_info_by_haplotype` is a dictionary that stores other haplotypes' position info at a specific aligned position.

Will need to fix this before I can continue with the formalization.

### panMANIA

The panMAN gives `Exceeded message traversal limit` error for SARS and HIV trees. The reason is that the default message traversal limit is too small.

To fix this, replace `uint64_t traversalLimitInWords = 8 * 1024 * 1024 * 128` in `/home/capnproto-c++-1.0.2/src/capnp/message.h` with `uint64_t traversalLimitInWords = 8 * 1024 * 1024 * 256`.

Then run command below to rebuild capnproto.
```
make clean && ./configure && make -j && make install
```

\
As a continuation of `2/1/2025`, I need to assess the alignment quality of the new panMANs. Writing and split the fasta files for nodes on the panMANs are currently taking too long for `SARS_20000` and `HIV20000`.

## 2/4/2025

### To-do

- [ ] Complete the alignment quality assessment for new panMANs.
- [ ] Make a figure that illustrates the effectiveness of pseudo-chaining scores at discerning between sequences.
- [x] ~~Make a figure that shows the correlation between the pseudo-chaining score and the alignment score.~~ *need to find a better way to illustrate this*
- [x] ~~Fix `position_info_by_haplotype` bug in `heuristic_v6.py`.~~ *Was redefining `pos` in debug print statements.*

### panMAMA

To make the pseudo-chaining score vs sequence similarity figure (refer to `2/3/2025`), I'm going to write a script that would give all 150-bp kmers from a random node (preferably containing no ambiguous bases, `OM857280.1`) and generate 50 mutated kmers from each kmer, each containg 1 - 50 mutations. Then I will calculate align them using minimap2 and pseudo-chaining and compare their scores. Scripts and data are placed in `/private/groups/corbettlab/alan/lab_notebook/panmama/pseudo_chaining-vs-sequence_similarity`.

The final figure is saved in `panmama/pseudo_chaining-vs-sequence_similarity/pseudo_chaining_vs_seq_similarity.png`. It doesn't look very good... A better to illustrate the effectiveness of pseudochaining is to show how it's able to discern between haplotypes, since we only care about if it's sensitive enough to tell different sequences apart.

## 2/5/2025

### To-do

- [ ] Finish writing the first draft of the specific aim section of my thesis proposal.
- [ ] Really do finish formalizing the method for calling consensus.

### Logistics and general work

#### Meeting with Russ
I talked a little about panMANIA then the consensus calling step with Russ. Russ gave a really good sugesstion on assuming there no back mutations. See notes below in panMAMA section.

### panMANIA

I just plotted the alignment differences for the new SARS and RSV panMANs (HIV is skipped because it's taking too long. Seems like the aligned sequences in the new HIV panMAN are very long). The new SARS panMAN looks pretty identical to the old SARs panMAN but **the new RSV panMAN looks **worse** than the old one**. See [panmania/evaluate_alignments](panmania/evaluate_alignments) for more details.

Maybe it is worth to also fix the alignments in panMAN for Aim 2.

Notes below need to be beautified and executed. This is a place holder so I don't forget

*I have an idea to identify problematic block states. Will try to implement and see how it goes.*
- *Record positions of mismatches in the panMAN alignment. Merge positions into ranges. Identify nearby ranges with identical lengths of mismatches.*
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
- [x] ~~Finish writing the first draft of the specific aim section of my thesis proposal.~~ *See [thesis_proposal/](thesis_proposal/) directory for the draft*

### Logistics and general work

I finished the first draft of my specific aims section for my thesis proposal. Aim 3 still needs a lot of work. See [thesis_proposal/](thesis_proposal/) section for the first draft.

### panMANIA

With the limited time I have today, I will write up the script for identifying problematic block states that cause large chunks of misalignment. See [panmania/evaluate_alignments/](panmania/evaluate_alignments/) for details

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

Finally finished formalizing consensus calling. See [panmama/consensus_calling/](panmama/consensus_calling/) for sample scripts.

Pseudocode:
```
for line in vcf_all:
  possible_alleles = all alleles with depth > abundance - tolerance

  split reads into groups by what other haplotypes they are assigned to

  if an allele-determining group (reads assigned to only the current haplotype) exists:
    assign the allele aligned in the allele-determining group
  else
    iterate through the groups of reads and calculate the likelihoods of the current haplotype having each possible allele
    
    split the groups of reads futher into groups by the number of haplotypes they are assigned to.
    
    within each group split by the number of haplotypes they are assigned to, select the subgroup with the highest likelihood ratio. This is the group's representative likelihood ratio

    accumulate the likelihood ratio for each possible allele across the groups

    if highest accumulated likelihood ratio - second highest accumluated likelihood ratio > threshold:
      assign allele with the highest accumulate likelihood ratio
    else:
      assign reference allele
```

### panMANIA


I fininally finished quantifying the prevalence of problematic block states (see [panmania/evaluate_alignments/](panmania/evaluate_alignments/#quantify-the-prevalence-of-problematic-block-states)). It looks like it's definitely worth fixing the misaligned blocks or representing them differently (a different layer of block coordinate?). Will talk to Russ about it.

Fixing this should also substantially improve the runtime of `panmap` and `panMAMA`.

## 2/11/2025

### To-do

- [ ] Characterize errors in consensus calling in panMAMA
- [x] ~~Develop a formal pipeline for evaluating panMAN internal alignments~~

### panMAMA

Today I will manually look through the output of consensus calling and characterize the errors. See [panmama/consensus_calling/](panmama/consensus_calling/) for detail.

Idea: perhaps I can also calculate the average sequence similarity of the reference and reads supporting each allele. 

I think I'll also assess and plot the threshold for assigning alleles.

### panMANIA

Steps for developing a formal pipeline for testing panMAN internal alignments.
  
1. Add a function in panMAN to output the ranges of all blocks

2. Add a function in panMAN to randomly selectly N pairs of nodes and output their aligned sequences and unaligned sequences

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

2. Partially incorporated read to reference sequence similairty into consensus calling. Currently support raw counts but can take position into account as well (counting reads mapped to the same position only once).

3. Still need to implement the 2nd round of ambiguity resolution.

4. Need to implement a way to handle indels.

## 2/19/2025

### To-do

- [ ] Finalize panMAMA consensus calling

### panMAMA

For the past few days, I've been improving panMAMA consensus calling. I think I've improved it enough to start finalizing the method.

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

- [x] ~~**Make a figure of read assignments for Russ**~~ *See [panmama/illustrating_figures_andor_diagrams/](panmama/illustrating_figures_andor_diagrams/)*

### panMAMA

I'm going to do some final optimization of panMAMA runtime. This will be done on silverbullet for efficiency.

Keeping the command and options here so I don't have to re-write them.

```
panmap example.panman example_R1.fastq example_R2.fastq --place-per-read --redo-read-threshold 0 --em-filter-round 2 --remove-threshold 0.01 --rounds-remove 5 --preem-filter-method mbc --save-kminmer-binary-coverage
```

After meeting with Richard Durbin, we got some good ideas to further improve accuracy and speed. We can compress homopolymers when sketching syncmers, and this could reduce syncmer errors due to homopolymer related sequencing errors. Additionally, we can implement a sorted and merged syncmer/kminmer index as direct hash lookup, while has a O(1) lookup time, is much slower than linear lookup of adjacent hash values.

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

I'm going to graph the distribution of the types of kminmer changes. Specifically, the number of kminmer changes for a read and the type of kminmer changes (Whether the change had added a kminmer to the chain, removed a kminmer from the chain, or didn't change the chain). See [panmama/run_time_optimization/](panmama/run_time_optimization/) for details.

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

I will install and test out a profiler (linaro map) today to see where panMAMA is spending most of its time. See [panmama/run_time_optimization/](panmama/run_time_optimization/) for more details.


## 3/9/2025

### To-do

- [x] Optimize panMAMA runtime and memory usage (high priority)

- [ ] Modify panmap build and src files to read new panMANs (high priority)

- [ ] Evaluate panMAMA using real data (high priority)

- [ ] Finalize panMAMA consensus calling (see notes on [2/19/2025](#2192025)) (high priority)

- [ ] Implement homopolymer compressed syncmers (low priority)

- [ ] Implement sorted and merged syncmer/kminmer index if possible (low priority)

### panMAMA

I think I've finished optimizing panMAMA runtime and memory usage for the time being. Now, panMAMA scores 1 million reads against 20,000-sample SARS tree in 5 minutes using 9GB of memory. EM step takes ~14 minutes.

Now to evaluate panMAMA using real data. See [panmama/real_data/](panmama/real_data/) for details.

## 3/18/2025

### To-do

- [ ] Evaluate panMAMA using real data (high priority)

- [ ] Implement coverage filter for panMAMA (high priority)

- [ ] Modify panmap build and src files to read new panMANs (high priority)

### panMAMA

I have been trying different ways to analyze the RSV data that Marc sent me. They are looking a little better now that I ran reads generated from RSVA priemrs and RSVB primers separately. Assuming Marc's assignment of mixture/isolated is correct, the false positives from panMAMA are because of low depth over all or high depth of likely contaminants at one or two loci. This is actually very useful information, as it's probably a good idea to consider coverage in the tool during demixing (set a filter to ignore haplotypes with coverage less than some threshold, [at least for amplicon samples]).

See my [notes](https://docs.google.com/spreadsheets/d/14gnGfWoBKk4Jzjnv00v8XbslElWe5N_IhrAk8pv9s_k/edit?usp=sharing) for details.

There are some cases I'm not sure when the samples would have 40-50% coverage with low depth in both RSV A and RSV B, I have these highlighted in bright yellow in the notes column. I also found two instances of "false negatives" where I had assigned them to be unmixed while you assigned them to be mixed:

```
RSV00174
RSV00220
```

After manual inspection of these two samples on the IGV, I actually found them to be unmixed (their alignments on IGV look like the false positives where they have high depth at only one or two loci). They are also highlighted in bright yellow.

It seems like a good idea for me to also incorporate kminmer overlap into prefiltering the nodes before EM (at least for amplicon samples). Most of the false positives are also cases where contaminant read stacks map to one or two loci on one of the RSV types.

Will stare at the google sheet for a while to see if I can find any other patterns.

## 3/23/2025

### To-do

- [ ] Evaluate panMAMA using real data (high priority)

- [ ] Implement coverage filter for panMAMA (high priority)

- [ ] Modify panmap build and src files to read new panMANs (high priority)

### panMAMA

I discussed with Marc our reasonings for determining whether the RSV samples are mixed or unmixed. My logic is that if there are truly both RSVA and RSVB present in sample, then the reads generated from RSVA primers should map well to RSVA ref and cover most of the genome, and the reads generated from RSVB primers should map well to RSVB ref and cover most of the genome. On the other hand, samples where reads only map well to a few regions on the genome likely contain contaminants that get amplified in the PCR process. Marc, on the other hand, offers an alternative interpretation that we should believe what the data is telling us:

>The target genomes in any of the samples (not just the mixed/co-infected ones) may not always be full-length (due to sample degradation over time(?)), and may not always be present in equimolar amounts (along the genome for one subtype). For me [Marc], I find it much harder to understand how RSVA or RSVB contaminants could be introduced in an essentially stochastic way among 192 samples all collected on different dates in a variety of locations and processed on different dates.  You might want to discuss these alternatives (and others) with Russ.

I will discuss this with Russ.

### panMANIA

Thanks to Faith Okamoto's amazing rotation project, I have some really good preliminary data and code to work with. I copied over Faith's google docs to my google drive in case they some how disappear from the internet.

I will generate some additional figures as preliminary data for my advancement.

See [panmania/imputation/](panmania/imputation/) for details.


