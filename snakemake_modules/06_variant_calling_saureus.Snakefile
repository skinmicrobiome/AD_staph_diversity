# -*- coding: utf-8
import os
from os.path import join
import sys
import glob
import pandas as pd
import csv



genomes_read_dir="data/00_preprocessing/processed/genomes/AD/"
tree_snippy_folder= "data/genomes_analysis/tree/saureus_all/" #ST5/"
reference="data/genomes/renamed/aatyw.fasta" #ference/saureus/GCF_000013465.1_ASM1346v1_genomic.fna"

id="_paperfinal" #more genomes


df_run=pd.read_csv("metadata_isolates.csv")
df_run=df_run[df_run["Species"]=="s__Staphylococcus aureus"]
df_run=df_run[df_run["QC"]=="Passed QC"]
genomes=df_run["run_accession"]



rule all_snp:
    input:
        tree_snippy_folder + "results/snp_tree.treefile",
        tree_snippy_folder + "results/snp_distancematrix.tsv"



rule run_pe_snippy:
    input:
         illumina_reads_1=ancient(genomes_read_dir+"{run}_1.fastq.gz"),
         illumina_reads_2=ancient(genomes_read_dir+"{run}_2.fastq.gz"),
         ref=reference
    output:
        vcf = tree_snippy_folder+"snippy/{run}/snps.vcf2",
        csv = tree_snippy_folder+"snippy/{run}/snps.csv2",
    threads: 4
    conda: "envs/snippy.yaml"
    params:
        outdir=tree_snippy_folder+"snippy/{run}/",
        mapqual = 60, #Snippy only uses uniquely mapping reads with mapping quality 60 (this means unique in BWA MEM).
        basequal = 20,  #It also only uses bases with Q20 or higher (1 in 100 error).
        mincov = 10,  # Minimum coverage of variant site (default '10')
        minfrac = 0.9, #Minumum proportion for variant evidence (default '0.9')
        run="{run}"
    shell:
        """
        rm -rf {params.outdir}
        snippy --cpus {threads} --outdir {params.outdir} --reference {input.ref} --pe1 {input.illumina_reads_1} --pe2 {input.illumina_reads_2} --mapqual {params.mapqual} --basequal {params.basequal} --mincov {params.mincov} --minfrac {params.minfrac}
        """



# collect snippy results per sample and find core snps
rule collect_snippy:
    input:
        expand(tree_snippy_folder+"snippy/{run}/snps.vcf", run=genomes),
    output:
        txt = tree_snippy_folder + "results/core.txt",
        vcf = tree_snippy_folder + "results/core.vcf",
        aln = tree_snippy_folder + "results/core.aln",
        fullaln = tree_snippy_folder + "results/core.full.aln",
    conda: "envs/snippy.yaml"
    params:
        prefix = tree_snippy_folder + "results/core",
        ref=reference,
        inp=expand(tree_snippy_folder+"snippy/{run}/", run=genomes),
    shell:"""
        snippy-core --ref {params.ref} --prefix {params.prefix} {params.inp} 
"""


rule clean_snippy:
    input:
        fullaln = tree_snippy_folder + "results/core.full.aln",
    output:
        fullaln = tree_snippy_folder + "results/core.clean.full.aln",
    message: "Cleaning core alignment"
    conda:
       "envs/snippy.yaml"
    params:
       dir= tree_snippy_folder + "results/"
    shell:
        "snippy-clean_full_aln {input} > {output}"



rule run_gubbins:
    input:
        fullaln = tree_snippy_folder + "results/core.clean.full.aln",
    output:
        filtered_fasta = tree_snippy_folder + "results/gubbins/gubbins.filtered_polymorphic_sites.fasta",
        gff=tree_snippy_folder + "results/gubbins/gubbins.recombination_predictions.gff",
        stat = tree_snippy_folder +  "results/gubbins/gubbins.per_branch_statistics.csv"
    conda:
       "envs/gubbins.yaml"
    params:
        #gubbins_dir = tree_snippy_folder + "results/gubbins",
        #prefix = tree_snippy_folder + "results/gubbins/gubbins",
    shell:
        """
        mkdir -p {params.gubbins_dir}  
        cd {params.gubbins_dir}
        run_gubbins.py --threads {threads} -p gubbins ../core.clean.full.aln 
        """



rule call_gubbines_snpsites:
    input:
        filtered_fasta = ancient(tree_snippy_folder + "results/gubbins/gubbins.filtered_polymorphic_sites.fasta"),
    output:
        corealn = tree_snippy_folder + "results/core.clean.aln",
    message: "Calling clean core alignment with rule call_gubbines_snpsites"
    conda:
       "envs/snippy.yaml"
    log:
       "logs/gubbins_snps.log"
    shell:
        "snp-sites -c {input.filtered_fasta} > {output.corealn}"



rule compute_filtered_distanceMatrix:
    input:
        tree_snippy_folder + "results/core.clean.aln"
    output:
        tree_snippy_folder + "results/snp_distancematrix.tsv"
    message: "Computing the filtered SNP distance matrix"
    conda:
        "envs/snippy.yaml"
    log:
       "logs/compute_distanceMatrix_filtered.log"
    shell:
        "snp-dists {input} > {output} 2> {log}"



rule compute_Tree_iqtree:
    input:
        tree_snippy_folder + "results/core.clean.aln"
    output:
        tree_snippy_folder + "results/snp_tree.treefile"
    message: "Computing a tree from core alignment with iqtree"
    conda:
        "envs/iqtree.yaml"
    threads: 14
    params:
       prefix = tree_snippy_folder + "results/snp_tree"
    shell:
       "iqtree -redo -s {input} -nt 14 -pre {params.prefix}"


