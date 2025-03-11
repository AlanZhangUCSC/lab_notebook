# Real data evaluation

## SARS data

I'm still waiting on Alex to send me the new CMakeLists.txt and other src files changed to support the new panMANs. Will work on RSV data for now.

## RSV data

Marc sent me 192 RSV samples. Since the reads are generated using amplicon sequencing, I need to trim the adapters first.

Since I'm working in the silverbullet server, all the paths are in `alan@silverbullet.ucsc.edu`.

First I'm going to manually check the primers. 

`RSVA.primer.bed` and `RSVB.primer.bed` contains the where the primers are located in the reference genome.

`$head -n 10 RSVA.primer.bed`

```
RS20000581	44	66	RSVA_1_LEFT	1	+
RS20000581	434	464	RSVA_1_RIGHT	1	-
RS20000581	359	385	RSVA_2_LEFT	2	+
RS20000581	749	773	RSVA_2_RIGHT	2	-
RS20000581	669	699	RSVA_3_LEFT	1	+
RS20000581	1057	1083	RSVA_3_RIGHT	1	-
RS20000581	990	1016	RSVA_4_LEFT	2	+
RS20000581	1366	1389	RSVA_4_RIGHT	2	-
RS20000581	1302	1330	RSVA_5_LEFT	1	+
RS20000581	1678	1706	RSVA_5_RIGHT	1	-

```

`RSVA.primer.tsv` and `RSVB.primer.tsv` contains the primer sequences.

`$head -n 10 RSVA.primer.tsv`

```
name	pool	seq	size	%gc	tm (use 65)
RSVA_1_LEFT	1	TTGGTTAGAGATGGGCAGCAAC	22	50	61.06
RSVA_1_RIGHT	1	AGGTCAAATCCAAGTAATTCAGATAATTGA	30	30	60.35
RSVA_2_LEFT	2	TGGCCTAATAGATGACAATTGTGAAA	26	34.62	59.55
RSVA_2_RIGHT	2	TGACCAGAAATGTAAATGTGGCCT	24	41.67	61.08
RSVA_3_LEFT	1	TCTAACCAGAGACATCATAACACATAAATT	30	30	59.95
RSVA_3_RIGHT	1	CCCATCTTTCATCTTATGTCTCTCCT	26	42.31	60.4
RSVA_4_LEFT	2	ACAACTTTATGCATAATCACACTCCA	26	34.62	59.78
RSVA_4_RIGHT	2	TGTTACATCCACTCCATTTGCCT	23	43.48	60.5
RSVA_5_LEFT	1	GCTATGTCTAGATTAGGAAGAGAAGACA	28	39.29	60.34
```

Let's inspect a random pair of samples `RSV00059_A_S83_L002_R1_001.fastq` and `RSV00059_A_S83_L002_R2_001.fastq` in `/scratch1/alan/rsv_reads/RSV00059/A`

Check if the primer sequence is in the reads. Using `RSVA_1_LEFT` primer as an example.

`$grep -n 'TTGGTTAGAGATGGGCAGCAAC' RSV00059_A_S83_L002_R1_001.fastq | head -n 10 | cut -c 1-50`

```
18:TTGGTTAGAGATGGGCAGCAACACTCTTCTAAAGGGCCAGGGAGTTT
38:TTGGTTAGAGATGGGCAGCAACACAAAGCAATGGCCGCCTCACGCAC
110:TTGGTTAGAGATGGGCAGCAACCCTTCAGCGTAAGAGGTCAGGGTG
210:TTGGTTAGAGATGGGCAGCAACGCGCTCGGCTCGCACAGCAGTTCC
458:TTGGTTAGAGATGGGCAGCAACTGCATCACATGCACACATAGGCAT
466:TTGGTTAGAGATGGGCAGCAACTGGGGGGATTGGTAGGGAGTGACT
522:TTGGTTAGAGATGGGCAGCAACTGGGGGGATTGGTAGGGAGTGACT
526:TTGGTTAGAGATGGGCAGCAACTGGGGGGATTGGTAGGGAGTGACT
686:TTGGTTAGAGATGGGCAGCAACTCCTTGCTGGTGGCCGTGGCGGCA
742:TTGGTTAGAGATGGGCAGCAACCATGTCGCTGTCCTGGAAGGGTGA
```

The primer is in the reads. It also appears that the primers are also anchored to the 5' end of the reads. After trying different primers, they all appear to be anchored to the 5' end of the reads...

I need split the primer.tsv into two files, one for the left primers and one for the right primers.

```
awk 'BEGIN {
  # Clear the files first
  system("rm -f RSVA_forward.primers.fa RSVA_reverse.primers.fa")
}
NR>1 {
  if ($1 ~ /LEFT/) {
    print ">"substr($1, 1, length($1)-5) >> "RSVA_forward.primers.fa"
    print $3 >> "RSVA_forward.primers.fa"
  } else if ($1 ~ /RIGHT/) {
    print ">"substr($1, 1, length($1)-6) >> "RSVA_reverse.primers.fa"
    print $3 >> "RSVA_reverse.primers.fa"
  }
}' /scratch1/alan/rsv_reads/other_files_you_might_need/RSVA.primer.tsv &

awk 'BEGIN {
  # Clear the files first
  system("rm -f RSVB_forward.primers.fa RSVB_reverse.primers.fa")
}
NR>1 {
  if ($1 ~ /LEFT/) {
    print ">"substr($1, 1, length($1)-5) >> "RSVB_forward.primers.fa"
    print $3 >> "RSVB_forward.primers.fa"
  } else if ($1 ~ /RIGHT/) {
    print ">"substr($1, 1, length($1)-6) >> "RSVB_reverse.primers.fa"
    print $3 >> "RSVB_reverse.primers.fa"
  }
}' /scratch1/alan/rsv_reads/other_files_you_might_need/RSVB.primer.tsv &
```

~~Now I have the primer files split into forward and reverse, `RSVA_forward.primers.fa` and `RSVA_reverse.primers.fa`; `RSVB_forward.primers.fa` and `RSVB_reverse.primers.fa`.~~

~~After closer inspection of the `RSV00059` sample, I found that most paired reads actually do NOT have the same pair of primers. I'm not quite sure what do about this so I will trim the primers in two different ways and see how they look. I will firs only trim the primers off reads that have the same pair of primers and call it `trimmed_strict`. Then I will trim the primers off all reads and call it `trimmed_loose`.~~

~~Then I will run panmap for the each trimmed sample. Since the RSV samples I got are not very similar to any of the samples in th panMANs, I will investigate whether panMAMA can discern RSV A and RSV B from the samples. To do that I will:~~

~~1. Run panmap for the trimmed samples~~

~~2. For each haplotype estimated to be present, use MAFFT to align it against both RSVA reference and RSVB reference~~

~~3. Assign the haplotype to be RSV A or RSV B based on to which reference it has the least distance~~

~~4. Accumulate the total proportions of RSV A and RSV B for each sample~~

~~See `run_panmap_parallel.sh` and `run_panmap.sh` in [/scripts](/scripts) for the scripts I used to do the things above.~~

I will actually use `ivar` to trim the primers. Following the pipeline that Marc sent me, I will first align the reads generated from the RSV A primers and RSV B primers to their respective reference genomes. Then I will run `ivar` to trim the primers using the primer bed files. Then I will convert the trimmed bam files to fastq files and merge the reads from the RSV A primers and RSV B primers.

`panMAMA` didn't work well. Assuming that most of these samples only contain one strain, panMAMA estimated many strains to be present. I think this is happening because none of the samples on the RSV panMAN is particularly close to the new samples Marc provided. There are ~100 differences between the closest sample on the panMAN and the new samples.

Nevertheless, I will evaluate whether panMAMA can correctly identify whether the sample contains RSVA, RSVB, or a mixture of both. To do that, I accumulated the proportions of the strains estimated to be present to either RSVA or RSVB depending on which reference genome the haplotype is closer to.

See  `run_panmap_ivar_trim.sh` in [/scripts](/scripts) for the script I used to do the things above.

For samples estimated to contain a mixture of both RSVA and RSVB, I wouldn't trust the estiimated proportions as I'm merging the reads from the two separate amplicon sequencing runs using different primers. They merely indicate that it's possible that the sample contains both RSVA and RSVB.

| RSV sample | RSV A | RSV B | Num reads after trimming|
| --- | --- | --- | --- |
| RSV00059 | 0.000 | 1.000 | 438946 |
| RSV00060 | 0.000 | 1.000 | 18684 |
| RSV00061 | 1.000 | 0.000 | 4740409 |
| RSV00062 | 0.771 | 0.229 | 25434 |
| RSV00063 | 0.000 | 1.000 | 13903347 |
| RSV00064 | 0.000 | 1.000 | 26460687 |
| RSV00065 | 0.000 | 1.000 | 122104 |
| RSV00066 | 0.000 | 1.000 | 457125 |
| RSV00067 | 0.000 | 1.000 | 1726 |
| RSV00068 | 0.000 | 1.000 | 72365 |
| RSV00069 | 0.000 | 1.000 | 8377485 |
| RSV00070 | 0.000 | 1.000 | 19099705 |
| RSV00071 | 0.000 | 1.000 | 1125467 |
| RSV00072 | 0.000 | 1.000 | 2788055 |
| RSV00073 | 0.000 | 1.000 | 7705223 |
| RSV00074 | 1.000 | 0.000 | 417496 |
| RSV00075 | 0.000 | 1.000 | 1113120 |
| RSV00076 | 0.000 | 1.000 | 17037316 |
| RSV00077 | 0.000 | 1.000 | 4 |
| RSV00078 | 1.000 | 0.000 | 823269 |
| RSV00079 | 0.000 | 1.000 | 1319508 |
| RSV00080 | 0.000 | 1.000 | 7417643 |
| RSV00081 | 0.000 | 1.000 | 3645831 |
| RSV00082 | 0.000 | 1.000 | 211655 |
| RSV00083 | 0.000 | 1.000 | 27274491 |
| RSV00084 | 0.000 | 1.000 | 50089 |
| RSV00085 | 1.000 | 0.000 | 39725 |
| RSV00086 | 1.000 | 0.000 | 20296068 |
| RSV00087 | 0.000 | 1.000 | 60835 |
| RSV00088 | 1.000 | 0.000 | 30294538 |
| RSV00089 | 1.000 | 0.000 | 145218 |
| RSV00090 | 0.000 | 1.000 | 538 |
| RSV00091 | 0.000 | 1.000 | 1929 |
| RSV00092 | 1.000 | 0.000 | 2944020 |
| RSV00093 | 0.000 | 1.000 | 225 |
| RSV00094 | 1.000 | 0.000 | 4945288 |
| RSV00095 | 0.000 | 1.000 | 3101985 |
| RSV00096 | 1.000 | 0.000 | 4601163 |
| RSV00097 | 1.000 | 0.000 | 18281 |
| RSV00098 | 0.000 | 1.000 | 111 |
| RSV00099 | 0.000 | 1.000 | 4389005 |
| RSV00100 | 1.000 | 0.000 | 37835029 |
| RSV00101 | 0.000 | 1.000 | 1602929 |
| RSV00102 | 0.000 | 1.000 | 962809 |
| RSV00103 | 0.000 | 1.000 | 10327818 |
| RSV00104 | 1.000 | 0.000 | 1253502 |
| RSV00105 | 1.000 | 0.000 | 850385 |
| RSV00106 | 1.000 | 0.000 | 9280933 |
| RSV00107 | 1.000 | 0.000 | 200 |
| RSV00108 | 1.000 | 0.000 | 29992 |
| RSV00109 | 0.000 | 1.000 | 3903991 |
| RSV00111 | 1.000 | 0.000 | 4677150 |
| RSV00112 | 0.412 | 0.588 | 17 |
| RSV00113 | 1.000 | 0.000 | 52803 |
| RSV00114 | 0.000 | 1.000 | 6609969 |
| RSV00115 | 0.000 | 1.000 | 2263551 |
| RSV00116 | 1.000 | 0.000 | 9694111 |
| RSV00117 | 0.000 | 1.000 | 1017 |
| RSV00118 | 0.000 | 1.000 | 4173374 |
| RSV00119 | 0.000 | 1.000 | 6391312 |
| RSV00120 | 0.167 | 0.833 | 13600 |
| RSV00121 | 0.000 | 1.000 | 7671 |
| RSV00122 | 1.000 | 0.000 | 35213901 |
| RSV00123 | 1.000 | 0.000 | 14133365 |
| RSV00124 | 1.000 | 0.000 | 14582704 |
| RSV00125 | 0.000 | 1.000 | 18060 |
| RSV00126 | 1.000 | 0.000 | 1908841 |
| RSV00127 | 1.000 | 0.000 | 5085 |
| RSV00128 | 0.000 | 1.000 | 425 |
| RSV00129 | 0.000 | 1.000 | 2228660 |
| RSV00130 | 0.056 | 0.944 | 124 |
| RSV00131 | 1.000 | 0.000 | 3400767 |
| RSV00132 | 1.000 | 0.000 | 46653 |
| RSV00133 | 0.000 | 1.000 | 3681804 |
| RSV00134 | 0.000 | 1.000 | 56602 |
| RSV00135 | 1.000 | 0.000 | 10880533 |
| RSV00136 | 0.000 | 1.000 | 21082 |
| RSV00137 | 0.000 | 1.000 | 1426166 |
| RSV00138 | 0.000 | 1.000 | 1063429 |
| RSV00139 | 1.000 | 0.000 | 2912915 |
| RSV00140 | 1.000 | 0.000 | 480211 |
| RSV00141 | 0.000 | 1.000 | 35716188 |
| RSV00142 | 1.000 | 0.000 | 39646791 |
| RSV00143 | 0.000 | 1.000 | 179171 |
| RSV00144 | 1.000 | 0.000 | 110223 |
| RSV00145 | 0.000 | 1.000 | 22738773 |
| RSV00146 | 0.000 | 1.000 | 667829 |
| RSV00147 | 0.081 | 0.919 | 37 |
| RSV00148 | 1.000 | 0.000 | 12151044 |
| RSV00149 | 0.000 | 1.000 | 8760907 |
| RSV00150 | 0.000 | 1.000 | 168 |
| RSV00151 | 0.000 | 1.000 | 1010829 |
| RSV00152 | 0.000 | 1.000 | 63200 |
| RSV00153 | 0.000 | 1.000 | 4391093 |
| RSV00154 | 1.000 | 0.000 | 145854 |
| RSV00159 | 0.222 | 0.778 | 9 |
| RSV00164 | 0.000 | 1.000 | 498796 |
| RSV00165 | 0.000 | 1.000 | 768372 |
| RSV00166 | 1.000 | 0.000 | 36341 |
| RSV00167 | 0.000 | 1.000 | 1815 |
| RSV00168 | 0.000 | 1.000 | 826257 |
| RSV00169 | 0.000 | 1.000 | 651327 |
| RSV00170 | 0.000 | 1.000 | 212363 |
| RSV00171 | 0.000 | 1.000 | 53969 |
| RSV00172 | 0.000 | 1.000 | 5621 |
| RSV00173 | 0.000 | 1.000 | 8962 |
| RSV00174 | 0.000 | 1.000 | 963797 |
| RSV00175 | 0.000 | 1.000 | 31740 |
| RSV00176 | 0.000 | 1.000 | 121626 |
| RSV00177 | 1.000 | 0.000 | 93871 |
| RSV00178 | 0.000 | 1.000 | 898732 |
| RSV00179 | 0.000 | 1.000 | 520236 |
| RSV00180 | 0.000 | 1.000 | 4388 |
| RSV00181 | 0.000 | 1.000 | 246730 |
| RSV00182 | 1.000 | 0.000 | 495 |
| RSV00183 | 1.000 | 0.000 | 2716 |
| RSV00184 | 0.000 | 1.000 | 899549 |
| RSV00185 | 1.000 | 0.000 | 1 |
| RSV00186 | 1.000 | 0.000 | 563693 |
| RSV00187 | 0.000 | 1.000 | 151269 |
| RSV00188 | 1.000 | 0.000 | 6792 |
| RSV00189 | 1.000 | 0.000 | 458369 |
| RSV00190 | 0.000 | 1.000 | 342242 |
| RSV00191 | 0.000 | 1.000 | 909361 |
| RSV00192 | 1.000 | 0.000 | 1253371 |
| RSV00193 | 1.000 | 0.000 | 586680 |
| RSV00194 | 1.000 | 0.000 | 127 |
| RSV00195 | 0.000 | 1.000 | 850760 |
| RSV00196 | 0.000 | 1.000 | 443491 |
| RSV00197 | 1.000 | 0.000 | 9111 |
| RSV00198 | 1.000 | 0.000 | 2223018 |
| RSV00199 | 0.000 | 1.000 | 2795 |
| RSV00200 | 0.000 | 1.000 | 92099 |
| RSV00201 | 0.000 | 1.000 | 50232 |
| RSV00202 | 0.000 | 1.000 | 255 |
| RSV00203 | 1.000 | 0.000 | 969737 |
| RSV00204 | 0.000 | 1.000 | 765421 |
| RSV00205 | 0.000 | 1.000 | 78034 |
| RSV00206 | 1.000 | 0.000 | 135850 |
| RSV00207 | 0.000 | 1.000 | 195 |
| RSV00208 | 0.000 | 1.000 | 88891 |
| RSV00209 | 0.000 | 1.000 | 1005854 |
| RSV00210 | 1.000 | 0.000 | 420249 |
| RSV00211 | 0.000 | 1.000 | 906043 |
| RSV00212 | 0.000 | 1.000 | 272675 |
| RSV00213 | 1.000 | 0.000 | 974363 |
| RSV00214 | 1.000 | 0.000 | 1570513 |
| RSV00215 | 0.000 | 1.000 | 622367 |
| RSV00216 | 1.000 | 0.000 | 285671 |
| RSV00217 | 0.000 | 1.000 | 576 |
| RSV00218 | 1.000 | 0.000 | 5429 |
| RSV00219 | 1.000 | 0.000 | 941102 |
| RSV00220 | 0.000 | 1.000 | 331259 |
| RSV00221 | 1.000 | 0.000 | 300096 |
| RSV00222 | 0.000 | 1.000 | 7443 |
| RSV00223 | 1.000 | 0.000 | 242428 |
| RSV00224 | 0.000 | 1.000 | 625889 |
| RSV00225 | 1.000 | 0.000 | 635745 |
| RSV00226 | 1.000 | 0.000 | 2029957 |
| RSV00227 | 0.000 | 1.000 | 655519 |
| RSV00228 | 0.000 | 1.000 | 286214 |
| RSV00229 | 1.000 | 0.000 | 1486842 |
| RSV00230 | 0.000 | 1.000 | 168826 |
| RSV00231 | 0.000 | 1.000 | 70052 |
| RSV00232 | 0.000 | 1.000 | 756285 |
| RSV00233 | 0.000 | 1.000 | 913834 |
| RSV00234 | 0.000 | 1.000 | 1018849 |
| RSV00235 | 1.000 | 0.000 | 338148 |
| RSV00236 | 1.000 | 0.000 | 787122 |
| RSV00237 | 1.000 | 0.000 | 16949 |
| RSV00238 | 0.000 | 1.000 | 76715 |
| RSV00239 | 0.000 | 1.000 | 1283341 |
| RSV00240 | 1.000 | 0.000 | 4775 |
| RSV00241 | 1.000 | 0.000 | 1645723 |
| RSV00242 | 0.000 | 1.000 | 1239849 |
| RSV00243 | 0.000 | 1.000 | 601460 |
| RSV00244 | 0.000 | 1.000 | 1214799 |
| RSV00245 | 1.000 | 0.000 | 346870 |
| RSV00246 | 0.000 | 1.000 | 1137802 |
| RSV00247 | 0.000 | 1.000 | 722866 |
| RSV00248 | 0.000 | 1.000 | 399359 |
| RSV00249 | 0.000 | 1.000 | 900240 |
| RSV00250 | 1.000 | 0.000 | 92330 |
| RSV00251 | 0.000 | 1.000 | 157 |
| RSV00252 | 1.000 | 0.000 | 164313 |
| RSV00253 | 0.000 | 1.000 | 363414 |
| RSV00254 | 0.000 | 1.000 | 323371 |
| RSV00255 | 1.000 | 0.000 | 1455189 |
| RSV00256 | 1.000 | 0.000 | 82921 |
| RSV00257 | 0.000 | 1.000 | 476717 |
| RSV00258 | 1.000 | 0.000 | 6 |
| RSV00259 | 1.000 | 0.000 | 163406 |


