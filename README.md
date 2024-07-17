# pipeline_miRNA_seq

This pipeline performs miRNAseq mapping with STAR and counts mapped reads to
STAR reference genome using STAR --quantmode TranscriptomeSAM GeneCounts as well 
counts mapped reads using featureCounts if a more specific GTF annotation file wants 
to be  used.


## Requirements

On top of the default CGAT setup, the pipeline requires the following.

Software:
    python (v3.8.12 with pysam v0.17.0 when built)  
    STAR (v2.7)
    featurecounts (v2.0.1)
    multiqc (v1.11)


## Configuration

The pipeline requires a configured :file: pipeline.yml file.

Make a directory with your project name. Configure the pipeline with 

python [path_to_repo]/pipeline_motif_analysis.py config

A pipeline.log and pipeline.yml file(s) will be added to your new directory.
Modify the pipeline.yml according to your project.


## Inputs

- Trimmed unpaired fasta files ending in .fastq.gz
- a STAR genome directory, pat hto be specified in pipeline.yml
- a GTF adapted for featureCounts for miRNAs


## Run pipeline

Run the pipeline with python [path_to_repo]/pipeline_miRNA_seq.py make full -v5.


## Ouputs

- star.dir directory with STAR outputs,
- counts.dir directory with featureCounts outputs
- multiqc outputs for both STAR and featureCounts
