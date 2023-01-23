# vim::w set ft=python:
from snakemake.utils import R
#import pandas as pd
import pandas as pd

QC_output_dir="data/genomes_analysis/QC/"
genomes_fna_dir="data/genomes/renamed/"  
extension=".fasta"
ROOT="/data/sahebkashafs2/AD/"


df_genomes=pd.read_csv("metadata_isolates.csv",sep=",") 
all_genomes=df_genomes["run_accession"]


(IDS, ) = glob_wildcards("data/genomes_analysis/QC/renamed/{id,[^/]+}")


rule all:
   input: 
          expand(QC_output_dir+"checkm/{dir}/checkm_metrics_parsed.tsv",dir=IDS),expand(QC_output_dir+"quast/sepidermidis/{id}/{id}.tsv", id=all_genomes), expand(QC_output_dir+"quast/{species}/quast_summary_numbcontig.tsv",species="sepidermidis"),expand(QC_output_dir+"quast/saureus/{id}/{id}.tsv", id=all_genomes), expand(QC_output_dir+"quast/{species}/quast_summary_numbcontig.tsv",species="saureus"), expand(QC_output_dir+"quast/{species}/quast_summary_genomefrac.tsv",species=["sepidermidis","saureus"])



#checkm for contaminants
rule checkm_folder:
    input: "data/genomes_analysis/QC/renamed/{dir}/"
    output: QC_output_dir+"checkm/{dir}/checkm_metrics.tsv",
    singularity:"containers/metawrap_w_checkm_database_latest.sif"
    params:
        ext="fasta",
        outdir=QC_output_dir+"checkm/{dir}/",
    threads: 10
    shell:
        """
        checkm data setRoot ~/checkm_database/
        rm -rf {params.outdir}
        checkm lineage_wf -t {threads} -x {params.ext} --tab_table -f {output} {input} {params.outdir}
        """


#parse checkm output
rule parse_checkm:
    input: QC_output_dir+"checkm/{dir}/checkm_metrics.tsv"
    output: QC_output_dir+"checkm/{dir}/checkm_metrics_parsed.tsv"
    params:
        checkm=QC_output_dir+"checkm/{dir}/checkm_metrics_tmp.tsv",
        checkm2=QC_output_dir+"checkm/{dir}/checkm_metrics_tmp.tsv2",
        checkm3=QC_output_dir+"checkm/{dir}/checkm_metrics_tmp.tsv3",
    shell:
        """
        ###sed '1d' {input} > {params.checkm3}
        scp {input} {params.checkm3}
        cut -f1,12,13,14 {params.checkm3} | tr '\t' ','>{params.checkm}
        sed 's/,/.fasta,/' {params.checkm}>{params.checkm2}
        echo -e "genome,completeness,contamination,strain_heterogeneity" | cat - {params.checkm2} > {output}
        rm {params.checkm}
        rm {params.checkm2}
        rm {params.checkm3}
""" 


rule merge_checkm:
    input: expand(QC_output_dir+"checkm/{dir}/checkm_metrics_parsed.tsv",dir=IDS)
    output: QC_output_dir+"checkm/all_checkm_metrics_parsed.tsv2" #EDIT
    params:
       output=QC_output_dir+"checkm/all_checkm_metrics_parsed.tsv.tmp"
    shell: """
cat {input} | grep -v 'genome'>{params.output}
echo -e "genome,completeness,contamination,strain_heterogeneity" | cat - {params.output} > {output}
rm {params.output}
"""


#check if all genomes are staph epi/aureus using quast - looking for genome fraction greater than 85%
rule quast:
   input:
     reference="data/genomes/reference/{species}/tree/reference.fasta",
     fasta=ancient(genomes_fna_dir+"{id}"+extension)
   output:QC_output_dir+"quast/{species}/{id}/{id}.tsv"
   params:
     outdir=QC_output_dir+"quast/{species}/{id}",
     out=QC_output_dir+"quast/{species}/{id}/report.tsv"
   shell:"""
   module load quast/5.0.2
   quast.py -r {input.reference} -o {params.outdir} {input.fasta}
   mv {params.out} {output}
"""


rule parse_quast:
   input: expand(QC_output_dir+"quast/{{species}}/{id}/{id}.tsv", id=all_genomes)
   output: 
      frac=QC_output_dir+"quast/{species}/quast_summary_genomefrac.tsv",
      numbcontig=QC_output_dir+"quast/{species}/quast_summary_numbcontig.tsv",
      N50=QC_output_dir+"quast/{species}/quast_summary_N50.tsv",
   params:
     indir=ROOT,
     dir=QC_output_dir+"quast/{species}/"
   shell:"""
   rm -rf {output.frac} {output.numbcontig} {output.N50}
   echo -e "genome\tgenome_frac" >> {output.frac}
   echo -e "genome\tN50" >> {output.N50}
   echo -e "genome\tnumb_contig" >> {output.numbcontig}
   cd {params.dir}
   for i in */*tsv; do {params.indir}scripts/parse_quast_genomefrac.sh ${{i}} >> quast_summary_genomefrac.tsv;done
   for i in */*tsv; do {params.indir}scripts/parse_quast_N50.sh ${{i}} >> quast_summary_N50.tsv;done
   for i in */*tsv; do {params.indir}scripts/parse_quast_numbcontig.sh ${{i}} >> quast_summary_numbcontig.tsv;done
"""


rule staph_epi_frac:
   input: QC_output_dir+"quast/{species}/quast_summary.tsv"
   output: "metadata_files/genomes/{species}.txt"
   run:
     R("""
   df=read.delim("{input}",sep="\t")
   df=df[df$genome_frac>=85,]
   df$genome_frac=NULL
   write.csv(df,"{output}",quote=F,row.names=F)
""")

