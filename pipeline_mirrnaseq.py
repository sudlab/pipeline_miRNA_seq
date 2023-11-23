import sys
import os
import sqlite3

from cgatcore import pipeline as P
import cgat.GTF as GTF
import cgatcore.iotools as IOTools
from ruffus import *


#First roung mapping
PARAMS = P.get_parameters("pipeline.yml")

@subdivide("*.fastq.gz",
         regex("(.+).fastq.gz"),
         [r"star_first.dir/\1.Aligned.sortedByCoord.out.bam", r"star_first.dir/\1.Unmapped.out.mate1"])
def star_mapping(infile, outfiles):
    job_memory = "8G"
    job_threads = 4
    mismax = PARAMS["star1_max_mismatch"]
    multi = PARAMS["star1_max_multi"]
    cover = PARAMS["star1_mis_cover"]
    genomeDir = PARAMS["star1_genome"]
    length = PARAMS["star1_min_length_match"]
    outname = P.snip(outfiles[0], "Aligned.sortedByCoord.out.bam")
    '''STAR mapping out of cgat'''
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

#Second round mapping
@subdivide(star_mapping,
         regex("star_first.dir/(.+).out.mate1"),
         [r"star_second.dir/\1.Aligned.sortedByCoord.out.bam", r"star_second.dir/\1.Unmapped.out.mate1"])
def star_mapping2(infile, outfiles):
    job_memory = "8G"
    job_threads = 4
    mismax = PARAMS["star2_max_mismatch"]
    multi = PARAMS["star2_max_multi"]
    cover = PARAMS["star2_mis_cover"]
    genomeDir = PARAMS["star2_genome"]
    length = PARAMS["star2_min_length_match"]
    outname = P.snip(outfiles[0], "Aligned.sortedByCoord.out.bam")
    '''STAR mapping out of cgat'''
    statement="""
    STAR --runThreadN 4
    --genomeDir /shared/sudlab1/General/mirror/genomes/STAR/hg38_noalt_junc93_49.dir
    --readFilesIn %(infile)s
    --outSAMtype BAM SortedByCoordinate
    --sjdbGTFfile /shared/sudlab1/General/annotations/hg38_noalt_ensembl93/ensembl.dir/geneset_all.gtf
    --quantMode TranscriptomeSAM GeneCounts
    --outFileNamePrefix %(outname)s
    --outFilterMultimapNmax %(multi)s
    --outFilterMatchNmin %(length)s
    --outFilterMismatchNmax %(mismax)s
    --alignEndsType EndToEnd 
    --alignIntronMax 1
    --outFilterMismatchNoverLmax %(cover)s
    """
    P.run(statement)


#Stats
@collate([star_mapping, star_mapping2],
         regex("(.+).dir/(.+).Aligned.sortedByCoord.out.bam"),
         r"\1.dir/resume_star_stats.txt")
def computeStarstats(infile, outfile):
    job_memory = "2G"
    script_path = os.path.join(os.path.join((os.path.dirname(__file__)),
                               "Rscripts",
                               "computeStarStats.R"))
    out_path = P.snip(outfile,"/resume_star_stats.txt")
    statement = '''
    Rscript %(script_path)s -o %(out_path)s
    '''
    P.run(statement)

#Counts
@transform([star_mapping],
           regex("star_(.+).dir/(.+).Aligned.sortedByCoord.out.bam"),
           r"counts_\1.dir/\2_countsUnique.txt")
def featureCountsUnique(infile, outfile):
    job_memory = "4G"
    job_threads = 4
    anno = PARAMS["gtf_gff"]
    #theinfile = str(infile)+""
    '''featureCounts'''
    statement="""
    featureCounts -a %(anno)s
    -T 4
    -O  --fracOverlap 0.8
    -o %(outfile)s %(infile)s
    """
    P.run(statement)


@transform([star_mapping],
           regex("star_(.+).dir/(.+).Aligned.sortedByCoord.out.bam"),
           r"counts_\1.dir/\2_countsPrimary.txt")
def featureCountsPrimary(infile, outfile):
    job_memory = "4G"
    job_threads = 4
    anno = PARAMS["gtf_gff"]
    #theinfile = str(infile)+""
    '''featureCounts'''
    statement="""
    featureCounts -a %(anno)s
    -t miRNA
    -g Alias
    -T 4
    --primary
    -O  --fracOverlap 0.8
    -o %(outfile)s %(infile)s
    """
    P.run(statement)

@collate([featureCountsUnique,featureCountsPrimary],
         regex("counts(.+).dir/.+_counts(.*).txt"),
         r"counts\1.dir/\2_resume_counts.txt")
def resumeCounts(infile, outfile):
    job_memory = "4G"
    script_path = os.path.join(os.path.dirname(__file__),
                               "Rscripts",
                               "computeReadCounts.R")
    statement = '''
    Rscript %(script_path)s -o %(outfile)s
    '''
    P.run(statement)

@follows(resumeCounts,computeStarstats)
def full():
    ''' Later alligator '''
    pass

P.main()
