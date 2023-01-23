# -*- coding: utf-8

import pandas as pd

genomes_read_dir="data/00_preprocessing/processed/genomes/AD/"
mlst_dir="data/genomes_analysis/mlst/staphepi/"

staphepi_genomes=pd.read_csv("metadata_files/genomes/sepidermidis.txt")
staphepi_genomes=staphepi_genomes["GCA"]

rule all:
   input: expand(mlst_dir+"srst2/{id}__mlst__Staphylococcus_epidermidis__results.txt", id=staphepi_genomes),mlst_dir+"summary/mlst_summ.txt"




rule get_db:
   output: 
     "data/databases/mlst/staphepi/Staphylococcus_epidermidis.fasta",
     "data/databases/mlst/staphepi/profiles_csv",
   params:
     outdir="data/databases/mlst/staphepi",
   conda:"envs/srst2.yaml" 
   shell:"""
cd {params.outdir}
getmlst.py --species "Staphylococcus epidermidis"
"""




rule srst2:
   input:
        r1=ancient(genomes_read_dir+"{sample}_1.fastq.gz"),
        r2=ancient(genomes_read_dir+"{sample}_2.fastq.gz"),
        db="data/databases/mlst/staphepi/Staphylococcus_epidermidis.fasta",
        db2="data/databases/mlst/staphepi/profiles_csv"
   output:mlst_dir+"srst2/{sample}__mlst__Staphylococcus_epidermidis__results.txt"
   params:
        db="data/databases/mlst/staphepi/Staphylococcus_epidermidis.fasta",
        db2="data/databases/mlst/staphepi/profiles_csv",
        out=mlst_dir+"srst2/{sample}",
   conda:"envs/srst2.yaml" 
   params:
   shell:"""
srst2 --input_pe {input.r1} {input.r2} --output {params.out} --log --mlst_db {params.db} --mlst_definitions {params.db2} --mlst_delimiter '_'
"""


rule summarize_srst2:
    input: 
       isolates=expand(mlst_dir+"srst2/{id}__mlst__Staphylococcus_epidermidis__results.txt", id=staphepi_genomes),
    output:mlst_dir+"summary/mlst_summ.txt"
    params:
       dir=mlst_dir+"srst2/",
       file="../summary/mlst_summ.txt"
    shell:"""
    rm -rf {output}
    cd {params.dir}
    echo -e 'Sample\tST'>>{params.file}
    for i in *__results.txt; do 
       tmp="${{i%%__mlst__Staphylococcus_epidermidis__results.txt}}";
       echo $tmp
       awk '{{print FILENAME,"\t",$2}}' $i | tail -n 1 >> {params.file}
    done 
    sed -i 's/__mlst__Staphylococcus_epidermidis__results.txt//g' {params.file}
    sed -i 's/ //g' {params.file}
"""



