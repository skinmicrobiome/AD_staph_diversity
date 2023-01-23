# vim::w set ft=python:
import os
from os.path import join
import sys
import glob
import pandas as pd
import csv


#metagenomes directory
metagenomesdir="data/00_preprocessing/processed/metagenomes/AD/" 

#specify files to work on
df_run=pd.read_csv("metadata_all_metagenomes.tsv",sep="\t")
metagenomes=df_run["Sample_ID"]

#specify output directory
readmapping_dir="data/metagenome_analyses/readmapping_SMGC_plusviral/" 


#specify database
SMGC_dir="data/databases/SMGC_plusviral/"
MAG_ISOLATES,=glob_wildcards(SMGC_dir+"fastas/{bins}.fa")



rule all_readmapping: 
      input: expand(readmapping_dir+"map2ref/{run}_ref-db_total.tab", run=metagenomes),
             readmapping_dir+"map2ref/bwa_counts-unique.csv"



rule rename_fasta:
     input: SMGC_dir+"fastas/{bins}.fa"
     output: SMGC_dir+"fastas_renamed/{bins}_ren.fasta"
     params:
          name="{bins}",
     shell:"""
          scripts/rename_multifasta_prefix.py -f {input} -p {params.name} > {output}
           """


rule cat_index_ref:
     input:expand(SMGC_dir+"fastas_renamed/{bins}_ren.fasta", bins=MAG_ISOLATES),
     output:join(SMGC_dir+"fastas_renamed/bwa-ref_name_vf/ref-db.fasta")
     params:
          indir=SMGC_dir+"fastas_renamed/bwa-ref_name_vf/"
     conda: "envs/readmapping.yaml"
     shell: """
            rm -rf {params.indir}
            mkdir -p {params.indir}
            cat {input}>{output}
            bwa index {output}
            """


rule map_reads:
     input:
         fwd = ancient(metagenomesdir+"{run}_1.fastq.gz"),
         rev = ancient(metagenomesdir+"{run}_2.fastq.gz"),
         bwa_ref_fasta=ancient(SMGC_dir+"fastas_renamed/bwa-ref_name_vf/ref-db.fasta")
     output: 
         readmapping_dir+"map2ref/{run}_ref-db_total.tab",
         readmapping_dir+"map2ref/{run}_ref-db_unique.tab"
     params:
        out_prefix=readmapping_dir+"map2ref/{run}"
     threads: 15 
     conda: "envs/readmapping.yaml"
     shell: """
            scripts/map2ref.sh -t {threads} -i {input.fwd} -n {input.rev} -r {input.bwa_ref_fasta} -o {params.out_prefix} -c contigs
            """


rule map_reads_single:
     input:
         fwd = metagenomesdir+"{run}.fastq.gz",
         bwa_ref_fasta=ancient(SMGC_dir+"fastas_renamed/bwa-ref_name_vf/ref-db.fasta")
     output:
         readmapping_dir+"map2ref/{run}_ref-db_total.tab2",
         readmapping_dir+"map2ref/{run}_ref-db_unique.tab2"
     params:
        out_prefix=readmapping_dir+"map2ref/{run}"
     priority: 100
     threads: 10
     conda: "envs/readmapping.yaml"
     shell: """
            scripts/map2ref.sh -t {threads} -i {input.fwd} -r {input.bwa_ref_fasta} -o {params.out_prefix} -c contigs
            """


rule summarize_read_mapping:
     input:
          expand(readmapping_dir+"map2ref/{run}_ref-db_total.tab", run=metagenomes),
          expand(readmapping_dir+"map2ref/{run}_ref-db_unique.tab", run=metagenomes),
     output: 
          readmapping_dir+"map2ref/bwa_counts-total.csv",
          readmapping_dir+"map2ref/bwa_counts-unique.csv",
     params:
        indir=readmapping_dir+"map2ref/"
     shell:"""
     scripts/parse_bwa.R -t {params.indir} -u {params.indir} -c 30 -m loose
     mv bwa* {params.indir}
"""



