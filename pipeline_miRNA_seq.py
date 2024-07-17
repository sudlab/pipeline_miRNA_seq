import sys
import os
import sqlite3

from cgatcore import pipeline as P
import cgat.GTF as GTF
import cgatcore.iotools as IOTools
from ruffus import *

#Load config file options
PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "../pipeline.yml",
     "pipeline.yml"])

#First roung mapping
PARAMS = P.get_parameters("pipeline.yml")

@subdivide("*.fastq.gz",
           regex("(.+).fastq.gz"),
           [r"star.dir/\1.Aligned.sortedByCoord.out.bam",
           r"star.dir/\1.Unmapped.out.mate1"])
def star_mapping(infile, outfiles):
    job_memory = "8G"
    job_threads = 4
    mismax = PARAMS["star_max_mismatch"]
    multi = PARAMS["star_max_multi"]
    cover = PARAMS["star_mis_cover"]
    genomeDir = PARAMS["star_genome"]
    length = PARAMS["star_min_length_match"]
    outname = P.snip(outfiles[0], "Aligned.sortedByCoord.out.bam")
    '''STAR mapping'''
    statement="""
    STAR --runThreadN 4
    --genomeDir %(genomeDir)s
    --readFilesIn %(infile)s
    --readFilesCommand zcat
    --outSAMtype BAM SortedByCoordinate
    --quantMode TranscriptomeSAM GeneCounts
    --outReadsUnmapped Fastx
    --outFileNamePrefix %(outname)s
    --outFilterMultimapScoreRange 0
    --outFilterMultimapNmax %(multi)s
    --outFilterScoreMinOverLread 0
    --outFilterMatchNminOverLread 0
    --outFilterMatchNmin %(length)s
    --outFilterMismatchNmax %(mismax)s
    --alignSJDBoverhangMin 1000
    --alignEndsType EndToEnd
    --alignIntronMax 1
    --outMultimapperOrder Random
    --outFilterMismatchNoverLmax %(cover)s
    """
    P.run(statement)


@collate(star_mapping,
         regex("(star.dir)/(.+).Aligned.sortedByCoord.out.bam"),
         r"\1/multiqc_report.html")
def multiQCstar(infile, outfile):
    """Run multiqc on star mapping"""
    job_memory = "2G"
    statement = '''
    export LANG=en_GB.UTF-8 && 
    export LC_ALL=en_GB.UTF-8 && 
	multiqc -f -n %(outfile)s --no-data-dir star.dir/ 
    '''
    P.run(statement)


@transform(star_mapping,
           regex("(.+).dir/(.+).Aligned.sortedByCoord.out.bam"),
           r"counts.dir/\2_counts.tsv")
def featureCounts(infile, outfile):
    """running featureCounts on mapped reads"""
    job_memory = "4G"
    job_threads = 4
    anno = PARAMS["featurecounts_gtf"]
    more_params = PARAMS["featurecounts_options"]
    #theinfile = str(infile)+""
    '''featureCounts'''
    statement="""
    featureCounts -a %(anno)s
    -T %(job_threads)s
    -t miRNA
    %(more_params)s
    -o %(outfile)s %(infile)s
    """
    P.run(statement)


@collate(featureCounts,
         regex("(counts.dir)/(.+).tsv"),
         r"\1/multiqc_report.html")
def multiQCfc(infile, outfile):
    """Run multiqc on star mapping"""
    job_memory = "2G"
    statement = '''
    export LANG=en_GB.UTF-8 && 
    export LC_ALL=en_GB.UTF-8 && 
	multiqc -f -n %(outfile)s --no-data-dir counts.dir/ 
    '''
    P.run(statement)


@follows(featureCounts, multiQCfc, multiQCstar)
def full():
    ''' Later alligator '''
    pass

P.main()
