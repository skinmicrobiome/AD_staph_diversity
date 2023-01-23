# vim: set ft=python:
import pandas as pd

amr_folder="data/genomes_analysis/amr/"
genomes_fna_dir="data/genomes/renamed/"


df_genomes=pd.read_csv("metadata_isolates.csv",sep=",")
all_genomes=df_genomes["run_accession"]


rule all:
   input:
       expand(amr_folder+"amrfinder/{id}.txt",id=all_genomes),
       expand(amr_folder+"VFDB/{id}.txt",id=all_genomes),
       expand(amr_folder+"VFDB_diamond_parsed/{id}_diamond_blastp_vfdb.tsv",id=all_genomes),


rule amr:
   input: ancient(genomes_fna_dir+"{id}.fasta")
   output:  amr_folder+"amrfinder/{id}.txt"
   group: "fast"
   shell: """
amrfinder -n {input} > {output}
"""

rule abricate_vfdb:
   input: genomes_fna_dir+"{id}.fasta"
   output:  amr_folder+"VFDB/{id}.txt"
   group: "fast"
   shell:"""
abricate --db vfdb {input}>{output}
"""

rule vfdb:
  input:
     faa="data/genomes_analysis/annotate/prokka/{id}/{id}.faa",
     db="data/databases/VFDB/VFDB_setB_pro.fas"
  output: amr_folder+"VFDB_diamond/{id}_diamond_blastp_vfdb.tsv"
  shell:"""
   diamond blastp -p {threads} -d {input.db} -q {input.faa} -o {output} --ultra-sensitive --max-target-seqs 1000 --evalue 0.000001
"""


rule vfdb_parsed:
   input: amr_folder+"VFDB_diamond/{id}_diamond_blastp_vfdb.tsv"
   output: amr_folder+"VFDB_diamond_parsed/{id}_diamond_blastp_vfdb.tsv"
   shell:"""
    Rscript scripts/parse_vfdb.R {input} {output}
"""




