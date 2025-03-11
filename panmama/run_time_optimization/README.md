# Runtime optimization

[Distribution of kminmer changes](#distribution-of-kminmer-changes)

[Profiling](#profiling)

## Distribution of kminmer changes

I'm going to graph the distribution of the types of kminmer changes. Specifically, the number of kminmer changes for a read and the type of kminmer changes (Whether the change had added a kminmer to the chain, removed a kminmer from the chain, or didn't change the chain).

I modified `pmi.cpp` to output some information ([affected_kminmers/sars20000_changes.txt](affected_kminmers/sars20000_changes.txt)) on the affected kminmers for each read at each node. I then wrote a script to parse this information and graph the distribution of the types of kminmer changes.

Added code below right before and after the main loop in `pmi.cpp` that scores the pseudochains.
```cpp
debugOut << "> " << node->identifier << "_before\n";
for (const auto& [readIndex, affectedSeedmerIndices] : readToAffectedSeedmerIndexVec) {
  debugOut << "@" << readIndex << "," << affectedSeedmerIndices.size() << "\n";
  for (const auto& match : reads[readIndex].matches) {
    debugOut << match.first << "," << match.second << " ";
  }
  debugOut << "\n";
  for (const auto& affectedSeedmerIndex : affectedSeedmerIndices) {
    debugOut << affectedSeedmerIndex << " ";
  }
  debugOut << "\n";
}
```

```cpp
debugOut << "> " << node->identifier << "_after\n";
for (const auto& [readIndex, affectedSeedmerIndices] : readToAffectedSeedmerIndexVec) {
  debugOut << "@" << readIndex << "," << affectedSeedmerIndices.size() << "\n";
  for (const auto& match : reads[readIndex].matches) {
    debugOut << match.first << "," << match.second << " ";
  }
  debugOut << "\n";
  for (const auto& affectedSeedmerIndex : affectedSeedmerIndices) {
    debugOut << affectedSeedmerIndex << " ";
  }
  debugOut << "\n";
}
std::cout << std::endl;
```

The read files used to generate the plots are 

- `affected_kminmers/reads/hiv20000_5hap-a_12000_rep1_R*.fastq`
- `affected_kminmers/reads/rsv4000_5hap-a_10000_rep1_R*.fastq`
- `affected_kminmers/reads/sars20000_5hap-a_10000_rep1_R*.fastq`



Blue bars - kminmers added to chain

Gold bars - kminmers removed from chain

Purple bars - kminmers unchanged on the chain (same as the previous node, not added or removed)

Column 5 means that out of all the reads that had 5 kminmers affected by the change in reference, about 1 million of them were added to the chain, about 7 million of them were removed from the chain, and a small fraction of them were unchanged.

SARS
![affected_kminmers/sars_kminmer_changes.png](affected_kminmers/sars_kminmer_changes.png)

RSV
![affected_kminmers/rsv_kminmer_changes.png](affected_kminmers/rsv_kminmer_changes.png)

HIV
![affected_kminmers/hiv_kminmer_changes.png](affected_kminmers/hiv_kminmer_changes.png)

Graphs below shows the kminmer changes by reads. Column 5 shows that out all of the reads that had 5 kminmer changes, about 0.2 million of them had all 5 kminmers added to the chain, about 1.4 million of them had all 5 kminmers removed from the chain, and a fraction of them caused zero or mixed changes to the chain.

SARS
![affected_kminmers/sars_kminmer_changes.png](affected_kminmers/sars_kminmer_changes_by_reads.png)

RSV
![affected_kminmers/rsv_kminmer_changes.png](affected_kminmers/rsv_kminmer_changes_by_reads.png)

HIV
![affected_kminmers/hiv_kminmer_changes.png](affected_kminmers/hiv_kminmer_changes_by_reads.png)

## Profiling

I will install and test out a profiler (linaro map) today to see where panMAMA is spending most of its time.

The profile is pretty easy to use. After following the installation instructions, run the command below to profile panMAMA.

```
map --profile panmap example.panman example_R1.fastq example_R2.fastq --place-per-read --redo-read-threshold 0 --em-filter-round 2 --remove-threshold 0.01 --rounds-remove 5 --preem-filter-method mbc --save-kminmer-binary-coverage
```

This outputs a `panmap_1p_1n_YYYY-MM-DD_HH-MM.map` file. Load this file into the linaro forge client and view the profile. 


**Ideas to improve runtime or memory usage**

| Costly operation | Time or memory costly | Idea to improve |
| --- | --- | --- |
|allScores[node->identifier] = allScores[node->parent->identifier]| memory and time | Change the way to store allScores. Instead of using an `unordered_map`, use a single object to store the scores and updates to the scores during a tree traversal. This should reduce both runtime (no longer need to copy the scores of aprents to children) and memory usage (no longer need to store the scores of all nodes). |
|hashToPositionsMap.find()|time|| 
|getMax()|time|Can likely make a condensed prob matrix by removing reads with identical scores and index them. Need to sketch it out to see if it's feasible.|