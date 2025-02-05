# Lab notebook

This notebook will track the progress of my work in the lab.

## 1/30/2025

### panMANIA

~~To build the panMAN using Dockerfile, see `build_panman` directory.~~
To build the panMAN see notes on `2/2/2025`.

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

For results, see `panmania/evaluate_alignments/README.md`.

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

### panMAMA 

### To-do

- [ ] Make a figure that shows the correlation between the pseudo-chaining score and the alignment score.
- [ ] Complete the alignment quality assessment for new panMANs. 
- [x] ~~Fix `position_info_by_haplotype` bug in `heuristic_v6.py`.~~ *Was redefining `pos` in debug print statements.*

### panMAMA

To make the pseudo-chaining score vs sequence similarity figure (refer to `2/3/2025`), I'm going to write a script that would give all 150-bp kmers from a random node (preferably containing no ambiguous bases, `OM857280.1`) and generate 50 mutated kmers from each kmer, each containg 1 - 50 mutations. Then I will calculate align them using minimap2 and pseudo-chaining and compare their scores. Scripts and data are placed in `/private/groups/corbettlab/alan/lab_notebook/panmama/pseudo_chaining-vs-sequence_similarity`.

For detail see `panmama/pseudo_chaining-vs-sequence_similarity/README.md`.

\
sequence similarity will be calculated using sam file generated from bwa mem.
```
python3 gen_seq.py OM8572801.fasta > simulated_reads.fastq

bwa mem OM8572801.fasta simulated_reads.fastq > simulated_reads.sam

samtools sort simulated_reads.sam > simulated_reads.sorted.bam

samtools calmd simulated_reads.sorted.bam OM8572801.fasta > simulated_reads.sorted.md.bam

samtools view simulated_reads.sorted.md.bam > simulated_reads.sorted.md.sam
```

\
pseudo-chaining score will be generated by adding code below to `pmi.cpp` after DFS scoring step.
```c++
  // print scores for OM8572801.1
  const auto& nodeScores = allScores.at("OM857280.1");
  std::ofstream nodeScoreOut("OM857280.1.scores.tsv");
  for (size_t i = 0; i < nodeScores.size(); ++i) {
    for (size_t readIdx : readSeedmersDuplicatesIndex[i]) {
      nodeScoreOut << readNames[readIdx] << "\t" << reads[i].seedmersList.size() << "\t" << nodeScores[i].first << "\n";
    }
  }
  nodeScoreOut.close();
```

After rebuilding panmap, run:
```
sbatch get_chain_scores.sh
```
Output is saved in `OM857280.1.scores.tsv`

\
The figure is saved in `panmama/pseudo_chaining-vs-sequence_similarity/pseudo_chaining_vs_seq_similarity.png`. It doesn't look very good... A better to illustrate the effectiveness of pseudochaining is to show how it's able to discern between haplotypes, since we only care about if it's sensitive enough to tell different sequences apart.


