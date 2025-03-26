# PanMANIA - Imputation

This is a sub-notebook for the annotation part of the panMANIA project. 

## 3/26/2025

I will start by exploring some of the notable existing annotation liftover tools

### UCSC genome browser liftOver

Installed from the [UCSC genome browser website](https://genome-store.ucsc.edu).

Here is the [example](https://genome.ucsc.edu/FAQ/FAQdownloads.html#liftOver) of running the liftOver commandline tool.

I didn't find any liftover chain files for the SARS-CoV-2 genome, so I couldn't test it. This actually proves the point of the `LiftOff` tool, which is that it can liftover from any genome to any genome.

### LiftOff

Liftoff paper summarized using `Claude 3.5`:

> Liftoff is a genome annotation transfer tool that improves upon traditional liftover methods by focusing on gene integrity rather than just coordinate conversion. It aligns entire gene sequences using Minimap2, then employs a graph-based algorithm to find mappings that preserve gene structure while maximizing sequence identity. Unlike tools like UCSC liftOver that process features independently, Liftoff maintains hierarchical relationships between exons, transcripts, and genes, preventing biological disruption. It also intelligently resolves paralogous genes and can identify additional gene copies not present in the reference. Testing showed Liftoff outperformed existing tools when transferring annotations between human genome versions and from human to chimpanzee.

Installed from conda:

```
conda install -c bioconda liftoff
```

Tried running the tool on the SARS-CoV-2 genome:

```
liftoff -g wuhCor1.gtf Germany_IMS_10293_CVDP_0A817728_2DA8_4F59_8351_DD8F8246B908_2021_OV210285_1_2021_11_07.unaligned.fasta wuhCor1.fa
```

and got an error

> GFF does not contain any gene features. Use -f to provide a list of other feature types to lift over.

Maybe it's failing because it's a GTF file and not a GFF3 file, although the github page says it supports GTF.

Let me try to convert the GTF file to GFF3.

Installing a gtf to gff3 converter:

```
conda install -c bioconda agat
```

and runnning:

```
agat_convert_sp_gxf2gxf.pl -g wuhCor1.gtf -o wuhCor1.gff3

liftoff -g wuhCor1.gff3 Germany_IMS_10293_CVDP_0A817728_2DA8_4F59_8351_DD8F8246B908_2021_OV210285_1_2021_11_07.unaligned.fasta wuhCor1.fa > Germany_IMS_10293_CVDP_0A817728_2DA8_4F59_8351_DD8F8246B908_2021_OV210285_1_2021_11_07.unaligned.gff3
```

It worked!


