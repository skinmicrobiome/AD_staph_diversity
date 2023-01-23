# vim::w set ft=python:
import pandas as pd

genomes_fna_dir="data/genomes/renamed/"
annotation_folder="data/genomes_analysis/annotate/"
GWAS_dir="data/genomes_analysis/GWAS_epi_2022/" 


df_run=pd.read_csv("metadata_isolates.csv")
df_run=df_run[df_run["Species"]=="s__Staphylococcus epidermidis"]
df_run=df_run[df_run["QC"]=="Passed QC"]
all_genomes=df_run["run_accession"]


outfiles_all=[]
outfiles_all.append(GWAS_dir+"panaroo/gene_presence_absence.Rtab")
outfiles_all.append(GWAS_dir+"eggnog/pangenome_done")
outfiles_all.append(GWAS_dir+"panaroo/core_gene_alignment.aln.treefile")


rule all:
   input: outfiles_all


rule make_inputfile:
   input:ancient(expand(annotation_folder+"prokka/{id}/{id}.gff",id=all_genomes))
   output:GWAS_dir+"input.txt"
   params:
     indir=GWAS_dir+"prokka"
   shell:"""
   ls {input} > {output}
"""


rule panaroo:
    input: GWAS_dir+"input.txt"
    output:
      GWAS_dir+"panaroo/gene_presence_absence.Rtab",
      GWAS_dir+"panaroo/gene_presence_absence.csv",
      GWAS_dir+"panaroo/gene_presence_absence_roary.csv",
      GWAS_dir+"panaroo/core_gene_alignment.aln", 
      GWAS_dir+"panaroo/pan_genome_reference.fa"
    params:
      outdir=GWAS_dir + "panaroo/"
    threads: 10
    singularity: "docker://quay.io/biocontainers/panaroo:1.2.9--pyhdfd78af_0"
    shell:"""
rm -rf {params.outdir}
panaroo -i {input} -o {params.outdir} -t {threads} -a core --clean-mode strict --merge_paralogs --remove-invalid-genes
"""


rule run_prodigal_pan2:
   input: GWAS_dir+"panaroo/pan_genome_reference.fa"
   output: GWAS_dir+"eggnog/pangenome_reference.faa"
   shell:"""
    prodigal -p 'meta' -a {output} -i {input}
"""


rule run_eggnog_pan2:
   input: GWAS_dir+"eggnog/pangenome_reference.faa"
   output: GWAS_dir+"eggnog/pangenome_done"
   params:
      outdir=GWAS_dir+"eggnog/",
      tmpdir=GWAS_dir+"eggnog/tmp/",
      MAG="pangenome"
   conda: "envs/py27.yaml"
   shell:"""
mkdir -p {params.outdir}
mkdir -p {params.tmpdir}
/data/sahebkashafs2/SOFTWARE/eggnog-mapper/emapper.py --cpu 10 -i {input} -m diamond -o {params.MAG} --output_dir {params.outdir} --temp_dir {params.tmpdir}
touch {output}
"""


rule iqtree_core:
   input: GWAS_dir+"panaroo/core_gene_alignment.aln",
   output: GWAS_dir+"panaroo/core_gene_alignment.aln.treefile"
   threads: 30
   shell:"""
iqtree -s {input} -nt {threads}
"""
