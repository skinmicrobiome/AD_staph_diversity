# vim: set ft=python:
import pandas as pd
from snakemake.utils import R


extension=".fasta"
mlst_dir="data/genomes_analysis/mlst/saureus/"

genomes_fna_dir="data/genomes/renamed/"
genomes_read_dir="data/00_preprocessing/processed/genomes/AD/"


df_run=pd.read_csv("metadata_isolates.csv")
df_run=df_run[df_run["Species"]=="Staphylococcus aureus"]
all_genomes=df_run["run_accession"]


rule all:
   input:  mlst_dir+"srst2/summary/mlst_summ.txt", expand(mlst_dir+"srst2/{id}__mlst__Staphylococcus_aureus__results.txt",id=all_genomes),
           expand(mlst_dir+"getmlst/ST/{gen}.txt",gen=all_genomes), mlst_dir+"getmlst/all_ST.txt"



rule copy_over:
   input: "../AD_analyses_project_vf/data/genomes/renamed/{run}.fasta"
   output: "data/genomes/renamed/{run}.fasta"
   shell:"scp {input} {output}"


#using reads

rule gzip:
   input:
      illumina_reads_1=genomes_read_dir+"{acc}_1.fastq",
      illumina_reads_2=genomes_read_dir+"{acc}_2.fastq",
   output:
      illumina_reads_1=genomes_read_dir+"{acc}_1.fastq.gz",
      illumina_reads_2=genomes_read_dir+"{acc}_2.fastq.gz",
   shell:"""
   gzip {input.illumina_reads_1}
   gzip {input.illumina_reads_2}
""" 

#saureus

rule get_db:
   output:
     "data/databases/mlst/saureus/Staphylococcus_aureus.fasta",
     "data/databases/mlst/saureus/profiles_csv",
   params:
     outdir="data/databases/mlst/saureus",
   conda:"envs/srst2.yaml"
   shell:"""
cd {params.outdir}
getmlst.py --species "Staphylococcus aureus"
"""

rule srst2:
   input:
        fwd=genomes_read_dir+"{id}_1.fastq.gz",
        rev=genomes_read_dir+"{id}_2.fastq.gz",
        db="data/databases/mlst/saureus/Staphylococcus_aureus.fasta",
        db2="data/databases/mlst/saureus/profiles_csv"
   output:mlst_dir+"srst2/{id}__mlst__Staphylococcus_aureus__results.txt"
   params:
        out=mlst_dir+"srst2/{id}"
   conda:"envs/srst2.yaml"
   threads: 1
   shell:"""
   module load srst2/0.2.0
   srst2 --input_pe {input.fwd} {input.rev} --output {params.out} --log --mlst_db {input.db} --mlst_definitions {input.db2} --mlst_delimiter '_'
"""


rule summarize_srst2:
    input:
       isolates=ancient(expand(mlst_dir+"srst2/{id}__mlst__Staphylococcus_aureus__results.txt", id=all_genomes)),
    output:mlst_dir+"srst2/summary/mlst_summ.txt"
    shell:"""
    rm -rf {output}
    echo -e 'Sample,ST'>>{output}
    for i in {input}; do
       sed -n '2p' $i | awk '{{ print $1 "," $2 }}' >> {output}
    done
"""


###get mlst

rule mlst:
   input:
      fa=ancient(genomes_fna_dir+"{gen}.fasta"),
   output: mlst_dir+"getmlst/ST/{gen}.txt"
   group: "fast"
   shell:"""
mlst --scheme saureus {input} | cut -f1,3 > {output}
"""

rule concatenate_ST:
   input:
      files=expand(mlst_dir+"getmlst/ST/{gen}.txt",gen=all_genomes)
   output: mlst_dir+"getmlst/all_ST.txt"
   params:
      out=mlst_dir+"getmlst/ST/"
   shell:"""
  rm -rf {output}
  cd {params.out}
  find ./ -type f | xargs cat > ../all_ST.txt
"""

