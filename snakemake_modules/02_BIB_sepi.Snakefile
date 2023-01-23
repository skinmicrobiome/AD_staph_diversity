# -*- coding: utf-8
import pandas as pd
import csv


#only use nih genomes for strain abundance

#output dir
BIB_v2_dir="data/metagenome_analyses/strain_abundance/sepidermidis_manual_NIH_subset/"

#metagenomes
df_run=pd.read_csv("metadata_all_metagenomes.tsv",sep="\t")
df_run=df_run[df_run["Study"]=="To be announced"]
metagenomes=df_run["Sample_ID"]
df_run=pd.read_csv("metadata_files_final/BIB_sepi.csv")
genomes=df_run["genome"]
metagenomesdir="data/00_preprocessing/processed/metagenomes/AD/"





rule all:
   input: BIB_v2_dir+"abundance_all_vf.tsv",expand(BIB_v2_dir+"abundance/{run}/{run}_abundance_vf.tsv",run=metagenomes)




rule copy_genomes:
   input: "data/genomes/renamed/{id}.fasta"
   output: BIB_v2_dir+"genomes/{id}.fasta"
   group: "fast"
   shell:"""
   scp {input} {output}
"""


rule merge_reads:
    input:
      fwd=metagenomesdir+"{sample}_1.fastq.gz",
      rev=metagenomesdir+"{sample}_2.fastq.gz",
    output:"/data/sahebkashafs2/AD_analyses_project_vf/data/00_preprocessing/processed/cat_metagenomes/{sample}.fastq.gz"
    shell:"""
   cat {input.fwd} {input.rev} > {output}
"""


rule progressivemauve_v2:
   input: expand(BIB_v2_dir+"genomes/{id}.fasta", id=genomes) 
   output: BIB_v2_dir+"full_alignment.xmfa2" 
   params:
     outdir=BIB_v2_dir
   shell:"""
cd {params.outdir}
progressiveMauve --output full_alignment.xmfa genomes/*fasta
"""


rule LCB_v2:
   input: BIB_v2_dir+"full_alignment.xmfa"
   output: BIB_v2_dir+"core_alignment_gapless.fasta2" 
   params:
     out= BIB_v2_dir
   shell:"""
module load bowtie
cd {params.out}
stripSubsetLCBs full_alignment.xmfa full_alignment.xmfa.bbcols core_alignment.xmfa 500 4
perl /data/sahebkashafs2/AD_analyses_project/scripts/xmfa2fasta.pl --file core_alignment.xmfa> core_alignment.fasta 
sed 's/-//g' core_alignment.fasta> core_alignment_gapless.fasta
bowtie2-build core_alignment_gapless.fasta core_alignment_gapless
"""


rule align_reads_v2:
   input:
     merged="/data/sahebkashafs2/AD_analyses_project_vf/data/00_preprocessing/processed/cat_metagenomes/{run}.fastq.gz",
     core=BIB_v2_dir+"core_alignment_gapless.fasta"
   output: 
      alpha=BIB_v2_dir+"samples_results/{run}/{run}.m_alphas",
      abund=BIB_v2_dir+"abundance/{run}/{run}_abundance.tsv",
      abundvf=BIB_v2_dir+"abundance/{run}/{run}_abundance_vf.tsv",
   params:
     name=BIB_v2_dir+"core_alignment_gapless",
     readout=BIB_v2_dir+"samples_results/{run}/{run}"
   shell:"""
module load bowtie
python scripts/BIB_analyse_reads.py {input.merged} {input.core} {params.name} {params.readout}
python scripts/BIB_parse_output.py {params.readout} > {output.abund}
awk '{{print FILENAME (NF?"\t":"") $0}}' {output.abund}>{output.abundvf}
"""


rule cat_align_v2:
   input:
      align=BIB_v2_dir+"full_alignment.xmfa",
      files=expand(BIB_v2_dir+"abundance/{run}/{run}_abundance_vf.tsv",run=metagenomes)
   output:
      abundvf=BIB_v2_dir+"abundance_all_vf.tsv", 
      deffile=BIB_v2_dir+"seq_deffile.txt",
   params:
      fold=BIB_v2_dir+"abundance/"
   shell:"""
cat {params.fold}/*/*_abundance_vf.tsv>{output.abundvf}
grep '#' {input.align}>{output.deffile}
"""


