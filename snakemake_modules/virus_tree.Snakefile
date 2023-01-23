# vim::w set ft=python:
import pandas as pd

DBDIR="data/databases/virus/virsearch_db/"
virus_dir="data/genomes_analysis/virus_tree_0.1_bigtree_finalpaper/" 


#all_runs=pd.read_csv("metadata_files/specifications/virus/virus_hits.txt")
#all_runs=all_runs["genome"]

all_runs,=glob_wildcards("data/genomes_analysis/virus/all_viruses/{run}.fa")

aln,=glob_wildcards(virus_dir+"core_genes/{run}_ren.faa")
#aln=["abfqu_5_14","adpuy_2_1","ERR1904012_1_54","acvas_2_58","ERR1753510_5_3"] #GCA_017695915.1_ASM1769591v1_genomic_2_33"]
extension="fna" #asta"
#viruses,=glob_wildcards("data/genomes_analysis/virus/all_viruses_indiv/{run}.fasta")

all_genomes=pd.read_csv("metadata_files/specifications/all_isolates/all_staph.txt")
all_genomes=all_genomes["genome"]

rule all:
   input: virus_dir+"trees/trimmed.concat_notrim.newick",virus_dir+"trees/concat_notrim.treefile"
#expand(virus_dir+"spacer_host/{id}.tsv", id=all_genomes)
#virus_dir+"trees/trimmed.concat_notrim.newick",virus_dir+"trees/concat_notrim.treefile"
#virus_dir+"eggnog/pangenome.emapper.annotations"
#virus_dir+"all_viruses/cdhit/filtered_predictions_nr_clstr.txt",
#         virus_dir+"trees/trimmed.concat_notrim.newick",virus_dir+"trees/concat_notrim.treefile"
#




#rule copy_viruses:
 #  input:"data/genomes_analysis/virus/all_viruses/{run}.fa"
   #output: virus_dir+"all_viruses/{run}.fasta"
   #group: "fast"
   #shell:"scp {input} {output}"



rule cat_viruses:
   input: expand("data/genomes_analysis/virus/all_viruses/{run}.fa",run=all_runs)
#expand(virus_dir+"all_viruses/{run}.fasta",run=all_runs)
   output: virus_dir+"all_viruses.fasta"
   shell:"""
   cat {input} > {output}
""" 



rule run_prodigal:
   input: ancient(virus_dir+"all_viruses.fasta")
   output:
     genes=virus_dir+"prodigal/genes.fa",
     faa=virus_dir+"prodigal/prot.faa",
   params:
     prodigal=virus_dir+"prodigal/"
   shell:"""
prodigal -a {output.faa} -d {output.genes} -i {input} -p meta
"""


# cluster sequences at 50% identity
rule cluster_sequences:
   input:
     faa=virus_dir+"prodigal/prot.faa"
   output:virus_dir+"mmseqs2_c50/mmseqs_cluster.tsv"
   params:
     outdir=virus_dir+"mmseqs2_c50"
   threads: 4
   shell:"""
rm -rf {params.outdir}
mkdir -p {params.outdir}
scripts/mmseqs_wf.sh -f {input} -o {params.outdir} -i 0.5 -c 0.5 -t 4
touch {output}
"""



rule mmseq:
   input: virus_dir+"mmseqs2_c50/mmseqs_cluster.tsv"
   output: virus_dir+"mmseqs2_c50/mmseqs_matrix.tsv" 
   shell:"""
scripts/mmseqs2matrix.py {input} {output}
"""



rule get_core:
   input: virus_dir+"mmseqs2_c50/mmseqs_matrix.tsv"
   output: virus_dir+"mmseqs2_c50/core_genes.txt"
   params:
     dir=virus_dir+"mmseqs2_c50/"
   shell:"""
   Rscript --vanilla scripts/gubaphage_core.R {params.dir}
"""



checkpoint extract_list:
   input:
     faa=virus_dir+"prodigal/prot.faa",
     core=virus_dir+"mmseqs2_c50/core_genes.txt",
     mmseq=virus_dir+"mmseqs2_c50/mmseqs_cluster.tsv"
   output:virus_dir+"core_genes/done" 
   params:
     core=virus_dir+"core_genes"
   shell:"""
rm -rf {params.core}
for i in $(cat {input.core}); do grep -w ${{i}} {input.mmseq} | cut -f2>${{i}}_list.txt; done
mkdir -p {params.core}
scp {input.core} {params.core}
mv *list.txt {params.core}
cd {params.core}
for i in *list.txt; do python /data/sahebkashafs2/SRP002480_singlerun/MAG_Snakemake_wf/scripts/select_seqs_by_IDs.py -i ../prodigal/prot.faa -d ${{i}} -o ${{i%%_list.txt}}.faa; done
for i in *faa; do python2 /data/sahebkashafs2/SRP002480_singlerun/MAG_Snakemake_wf/scripts/rename_multifasta_trim.py -f ${{i}} -s "_" > ${{i%%.faa}}_ren.faa; done
touch done
"""



checkpoint run_indiv:
  input:
     done=virus_dir+"core_genes/done", 
     faa=virus_dir+"core_genes/{run}_ren.faa"  
  output: 
    aln=virus_dir+"core_genes/two_{run}_ren.aln",
    trim=virus_dir+"core_genes/two_{run}_ren.trimmed.aln"
  shell:"""
  module load muscle/3.8.31
  muscle -in {input.faa} -out {output.aln}
  trimal -automated1 -in {output.aln} -out {output.trim}
"""
 


rule find_duplicates:
  input: virus_dir+"core_genes/two_{run}_ren.aln"
  output: virus_dir+"duplicates/{run}.txt"
  shell:"""
  python3 scripts/find_dup.py {input} {output}
"""




rule gather_duplicates:
  input: expand(virus_dir+"duplicates/{run}.txt",run=aln)
  output: virus_dir+"duplicates.txt"
  shell:"""
  cat {input} | cut -f2 -d '>' | sort | uniq > {output}
"""


rule genome_list:
  input: expand(virus_dir+"core_genes/two_{run}_ren.aln",run=aln)
  output: virus_dir+"genomelist"
  shell:"""
  grep -e '>' {input} | cut -f2 -d ">" | sort | uniq > {output}
"""


rule remove_duplicates:
  input:
     dup=virus_dir+"duplicates.txt",
     fa=virus_dir+"core_genes/two_{run}_ren.aln"
  output: virus_dir+"nodup_notrim/two_{run}_ren.aln"  #remove this
  shell:"""
 python2 scripts/filter_fasta_by_header.py -d {input.dup} -i {input.fa} -o {output}
"""



rule concat:
   input: expand(virus_dir+"nodup_notrim/two_{run}_ren.aln", run=aln),virus_dir+"core_genes/done"
   output: virus_dir+"concat_alignment.aln"
   params:
     dir=virus_dir+"nodup_notrim/"
   shell:"""
  cd {params.dir}
  concat -f . -e _ren.aln --Prefix "two" -o ../concat_alignment.aln -N 2
"""


rule tree:
  input: virus_dir+"concat_alignment.aln"
  output: virus_dir+"trees/trimmed.concat_notrim.newick"
  shell:"""
  FastTree -gamma -lg {input} > {output}
"""

rule iqtree:
  input: virus_dir+"concat_alignment.aln"
  output: virus_dir+"trees/concat_notrim.treefile"
  message: "Computing a tree from core alignment with iqtree"
  conda:
        "envs/iqtree.yaml"
  threads: 40
  params:
       prefix =  virus_dir+"concat_notrim"  
  shell:"""
  iqtree -redo -s {input} -nt 40 -pre {params.prefix}
"""



rule run_eggnog:
  input: virus_dir+"mmseqs2_c50/mmseqs_cluster_rep.fa"
  output: virus_dir+"eggnog/pangenome.emapper.annotations"
  params:
    outdir=virus_dir+"eggnog/",
    tmpdir= virus_dir+"eggnog/tmp/",
    MAG="pangenome"
  conda: "envs/py27.yaml"
  shell:"""
rm -rf {params.outdir} {params.tmpdir}
mkdir -p {params.outdir}
mkdir -p {params.tmpdir}
/data/sahebkashafs2/SOFTWARE/eggnog-mapper/emapper.py --cpu 30 -i {input} -m diamond -o {params.MAG} --output_dir {params.outdir} --temp_dir {params.tmpdir}
touch {output}
"""


rule vfdb:
  input: 
     faa=virus_dir+"mmseqs2_c50/mmseqs_cluster_rep.fa",
     db="data/databases/VFDB/VFDB_setB_pro.fas"
  output: virus_dir+"diamond_blastp_vfdb.tsv"
  shell:"""
   diamond blastp -p {threads} -d {input.db} -q {input.faa} -o {output} --ultra-sensitive --max-target-seqs 1000 --evalue 0.000001
"""



rule crispr:
   input: "data/genomes/renamed/{id}.fasta"
   output: virus_dir+"CRISPRCasFinder/{id}/rawCRISPRs.fna"
   params:
     indir=virus_dir+"CRISPRCasFinder/{id}/",
     outdir=virus_dir+"CRISPRCasFinder/{id}/out",
     file="{id}.fasta"
   singularity: "tools/CrisprCasFinder.simg"
   group: "fast"
   priority: 1000
   shell:"""
   rm -rf {params.indir}
   mkdir -p {params.indir}
   scp {input} {params.indir}
   cd {params.indir}
   perl /usr/local/CRISPRCasFinder/CRISPRCasFinder.pl -so /usr/local/CRISPRCasFinder/sel392v2.so -cf /usr/local/CRISPRCasFinder/CasFinder-2.0.3 -drpt /usr/local/CRISPRCasFinder/supplementary_files/repeatDirection.tsv -rpts /usr/local/CRISPRCasFinder/supplementary_files/Repeat_List.csv -out {params.outdir} -in {params.file} -cpuM {threads}
scp */raw* .
"""


rule blastn_compare_spacers:
    input:
      MAG=virus_dir+"all_viruses/cdhit/nr_predictions.fa",
      spacer=virus_dir+"CRISPRCasFinder/{id}/rawCRISPRs.fna"
    output:
      sum= virus_dir+"spacer_host/{id}.tsv"
    group: "fast"
    shell:"""
        blastn -query {input.spacer} -subject {input.MAG} -outfmt 6 -max_target_seqs 5000 -evalue 1e-7 -gapopen 10 -gapextend 2 -reward 1 -penalty -1 -word_size 5 -out {output.sum}
        """


