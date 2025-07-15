After meeting with Liz 7/3:

- Basically the bigwig plots for R2 are *close* but not exactly what we're trying to do. 

"Raw PRO-seq reads were trimmed using cutadapt (-a TGGAATTCTCGG -A
GATCGTCGGACT -m 26 -q 10) to remove Illumina primers, and UMItools extract (-p NNNNNN --bc-pattern2=NNNNNN) to extract the 6nt Unique Molecular Identifiers (UMIs). Trimmed reads were aligned to the reference sequence using Bowtie2 (--very- sensitive). Samtools sort was used to sort and compress the Bowtie2 output. bedtools bam2bed was used to create bed files for the reference sequence. The bed files were read into Python (ver.3.10.12) for further processing. Since PRO-seq’s R1 contains the RNA 3’ end, reads aligning to the “-” strand were the reads aligned to the reference sequence. The first base of each R1 read aligned to the “-” strand were counted and plotted to generate RNA Pol II traces. The traces were binned every 10nt to reduce noise and were plotted with 95% confidence intervals.
RNA Pol II pausing indices (PI) were also calculated using the 3’ end data, where PI is the ratio of total 3’ reads in [TSS-50 – TSS+50] and in [TSS+50 – TSS+250]. Average read is defined by total counts in a window, divided by window size and was calculated for the following windows for the reporter gene: TSS–TSS+150, TSS+150– TES-150, TES-150–TES, TES–TES+150, and TES+150–TES+500. Normalized average reads were calculated by normalizing average reads in each window by average reads in [TSS+150 – TES-150] window."

- essentially we're trying to get the *polymerase traces* not the whole reads
- difference between sense vs antisense strands, R1 vs R2 and 5' vs 3' is very important. 
- We want to take the R1 reads which contain the RNA 3' end and use the reads aligning to the "-" strand. Then take the first base of each R1 read and count and plot them to generate these RNA Pol II traces. 
- The "pause indices" can be calculated by total counts in a window 

- Our RNA expression background across samples is roughly equal. I could run DeSeq2 Analysis and normalize with something like Salmon or to the Drosphilia spike in. Use feature count and estimate transcript abundance. However this is low priority so I'm shelving it for now. 