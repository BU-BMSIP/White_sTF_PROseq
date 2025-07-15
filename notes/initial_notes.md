Measuring kinetics of a synthetic transcription factor SynZiFTER 2.0 (vs 1.0)


Cell Bio background:

Macrophages are immune cells, polarize to 2 subtypes
M1 are pro-inflammation and show anti-tumor activity
M2 treat anti-inflammatory diseases
Cancer hijacks, and there are issues with current treatments of antibodies, AAV gene therapy, peptides and small molecule
3 tools developed for inducible polarization:
engineered STAT effector proteins that can drive macrophage polarization to M1 or M2 state
synZiFTRs 2.0, an optimized human gene regulation toolkit built upon our previous work, enabling drug-inducible and enhanced transgene expression with minimal basal activity
monocyte-targeting nanoparticles for cell-specific delivery of synZiFTR induction agents.

Synthetic Biology background / engineering

On/off switch enables more effective cell therapies
CAR-M is like CAR-T, but aimed at solid tumors and still pre-clinical
Necessary to have certain features
No “leakiness” 
Low toxicity
Strong activation
Easily tunable 

Experimental tool PROSeq

PROseq:
Cells are permeabilized, followed by a nuclear run-on assay
Engaged RNA polymerase complexes incorporate a biotinylated nucleotide into the nascent RNA chain 
Labeled RNAs are enriched using streptavidin beads and sequencing of the final libraries reveals RNA 3’ ends at single nucleotide resolution.

PRO-Seq illuminates process of elongation: RNAP2 activities such as pausing, pause release, RNA synthesis and termination are rate-related and time-sensitive
PRO-Seq is particularly suitable to gain insights in TIMs’ effects on not only stimulating transcriptional activity, but also mediating the kinetics of different stages of elongation
Quantifying the average counts of 3’ end reads in different regions on the gene body, such as promoter proximal, gene body, distal and post termination
Normalize the average counts of 3’ end reads in each region with that of the gene body. 
TIMs should be able to facilitate efficient elongation by increasing RNAP2 activity level by modulating RNAP2 dynamics such as post-termination clearance


Data Analysis Particulars:

Quantifying the average counts of 3’ end reads in different regions on the gene body, such as promoter proximal, gene body, distal and post termination
Normalize the average counts of 3’ end reads in each region with that of the gene body. 
TIMs should be able to facilitate efficient elongation by increasing RNAP2 activity level by modulating RNAP2 dynamics such as post-termination clearance

*Note this is actually wrong*, what we want is to measure the exact polII trace, ie single base pair resolution of R1 - strand reads mapped to plasmid. 

Pipeline bullshit:

Liz built a pipeline out of bash spit and grit
Danko lab has a pipeline but only for general PROseq: https://github.com/Danko-Lab/proseq2.0
