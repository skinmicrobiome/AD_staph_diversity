# vim::w set ft=python:
import pandas as pd




SCC_dir="data/genomes_analysis/SCC/"
genomes_fna_dir="data/genomes/renamed/"

rule all:
   input: SCC_dir+"SCC_types.txt",


#conda activate staphopia-sccmec
rule sccmec:
   output:SCC_dir+"SCC_types.txt"
   params:
      indir=genomes_fna_dir,
   shell:"""
staphopia-sccmec --assembly {params.indir} --ext fasta > {output}
"""



