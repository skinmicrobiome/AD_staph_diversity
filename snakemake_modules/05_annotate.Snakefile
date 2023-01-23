# vim::w set ft=python:


import pandas as pd
from os.path import join

genomes_fna_dir="data/genomes/renamed/"
annotation_folder="data/genomes_analysis/annotate/"

df=pd.read_csv("metadata_isolates.csv")
all_staphaureus=df["run_accession"]

rule all:
    input: expand(annotation_folder+"prokka/{isolate}/{isolate}.gff",isolate=all_staphaureus)  


rule move_to_AD:
     input: "data/genomes_analysis/assemblies/{id}/contigs.fasta"
     output: "data/genomes/AD/{id}.fasta"
     shell:"""
     scp {input} {output}
"""


rule rename_fa:
    input:
        isolate=ancient(genomes_fna_dir+"{isolate}.fasta")
    output:"data/genomes/renamed/{isolate}.fasta"
    params:
        name="{isolate}"
    group: "fast"
    shell:
        """
        scripts/rename_multifasta_prefix.py -f {input.isolate} -p {params.name} > {output}
        """


rule prokka:
    input:ancient("data/genomes/renamed/{isolate}.fasta")
    output: 
        annotation_folder+"prokka/{isolate}/{isolate}.gff",
        annotation_folder+"prokka/{isolate}/{isolate}.gbk",
    singularity:
        "docker://staphb/prokka:1.14.5"
    params:
        out_prokka=annotation_folder+"prokka/{isolate}",
        prefix="{isolate}"
    shell:
        """
        prokka {input} --kingdom Bacteria --outdir {params.out_prokka} \
        --prefix {params.prefix} --force --locustag {params.prefix} --cpus {threads}
        """

