# vim::w set ft=python:
import pandas as pd

metagenomesdir="data/00_preprocessing/processed/metagenomes/AD/" 


df_run=pd.read_csv("metadata_all_metagenomes.tsv",sep="\t")
metagenomes=df_run["Sample_ID"]


output_dir_saureus="data/metagenome_analyses/strain_tracking/saureus/"

#parameters
SOURMASH_DATABASE_THRESHOLD_BP =  1e5
SOURMASH_COMPUTE_KSIZES = '61'
SOURMASH_COMPUTE_SCALED = '10000'
param_str="k="+SOURMASH_COMPUTE_KSIZES+",scaled="+SOURMASH_COMPUTE_SCALED+",abund"


all_gen,=glob_wildcards(output_dir_saureus+"genomes/all_genomes/{gen}.fna")


rule gather_strains:
    input: output_dir_saureus+"summary/genbank_all.csv", expand(output_dir_saureus + "genbank/{sample}.x.genbank.gather.csv", sample=metagenomes)
#output_dir_saureus+"genomes/done_move", output_dir_saureus+"genomes/all_ST.txt",  output_dir_saureus+"summary/genbank_all.csv", expand(output_dir_saureus + "genbank/{sample}.x.genbank.gather.csv", sample=metagenomes)



rule download_acc:
    output:output_dir_saureus+"genomes/done",
    threads: 1
    resources:
      mem_mb=30000,
      runtime=300 #minutes
    params:
        out=output_dir_saureus+"genomes/",
        section="genbank"
    conda: "envs/ncbi_gen_downland.yaml"
    shell: """
     mkdir -p {params.out}
     cd {params.out} 
     ncbi-genome-download --genera "Staphylococcus aureus" --section {params.section} --formats fasta bacteria
     touch done 
"""


rule move_over:
     input: output_dir_saureus+"genomes/done",
     output: output_dir_saureus+"genomes/done_move"
     params:
       dir=output_dir_saureus+"genomes/genbank/bacteria/",
       outdir=output_dir_saureus+"genomes/all_genomes/"
     log: "logs/saureus/move.log"
     shell:"""
     mkdir {params.outdir}
     scp {params.dir}*/*_genomic.fna.gz {params.outdir} 2> {log}
     touch {output}
"""


checkpoint record_genomes:
     input: output_dir_saureus+"genomes/done_move"
     output: output_dir_saureus+"genomes/all_genomes.txt"
     params:
       outdir=output_dir_saureus+"genomes/all_genomes/"
     shell:"""
     ls {params.outdir}/* > {output}
"""


rule mlst:
   input:
      done=output_dir_saureus+"genomes/all_genomes.txt", 
      fa=output_dir_saureus+"genomes/all_genomes/{gen}.fna.gz"
   output: output_dir_saureus+"genomes/ST/{gen}.txt"
   params:
      fa=output_dir_saureus+"genomes/all_genomes/{gen}.fna"
   log: "logs/saureus/mlst/{gen}.log"
   shell:"""
mlst --scheme saureus {params.fa} | cut -f1,3 > {output}
"""


rule concatenate_ST:
   input: 
      done=output_dir_saureus+"genomes/all_genomes.txt",
      files=expand(output_dir_saureus+"genomes/ST/{gen}.txt",gen=all_gen)
   output: output_dir_saureus+"genomes/all_ST.txt"
   params:
      out=output_dir_saureus+"genomes/ST/"
   shell:"""
  rm -rf {output}
  cd {params.out}
  find ./ -type f | xargs cat > ../all_ST.txt
"""


#dplyr subset and pick ST with lowest number of characters
rule subset_ST:
   input: output_dir_saureus+"genomes/all_ST.txt"
   output: output_dir_saureus+"genomes/subset_ST.txt"
   params:
     dir=output_dir_saureus+"genomes/subset_ST" 
   shell:"""
   rm -rf {params.dir}
   mkdir -p {params.dir}
   Rscript --vanilla scripts/subset_STs_strains.R {input} {output}
   scp $(cut -f1 -d , {output}) {params.dir} 
"""


rule dereplicate_genome_db:
    input: output_dir_saureus+"genomes/subset_ST.txt"
    output: output_dir_saureus+"genomes/final_drep_2/done"
    params:
      indir=output_dir_saureus+"genomes/subset_ST/",
      outdir=output_dir_saureus+"genomes/final_drep_2/"
    singularity:
        "shub://sskashaf/MAG_wf_containers_2021:drep"
    shell:"""
   rm -rf {params.outdir} 
   #gunzip {params.indir}*.fna.gz
   dRep dereplicate -p 50 {params.outdir} -g {params.indir}/*.fna -pa 0.98 -sa 0.995 -nc 0.60 -cm larger --ignoreGenomeQuality
   touch {output}
"""  


#sketch genome database
rule sketch_genomes:
    input: output_dir_saureus+"genomes/final_drep/done"
    output: output_dir_saureus+"genomes/final_drep/dereplicated_genomes/genome.sig" 
    params:
      dir=output_dir_saureus+"genomes/final_drep/dereplicated_genomes/",
      p=param_str
    shell:"""
cd {params.dir}
sourmash sketch dna -p {params.p} * -o genome.sig
"""


# compute sourmash signature for reads
rule sketch_reads:
    input:
        r1 = metagenomesdir+"{sample}_1.fastq.gz",
        r2 = metagenomesdir+"{sample}_2.fastq.gz"
    output:
        sig = "data/metagenome_analyses/strain_tracking/sigs/{sample}.abundtrim.sig"
    conda: "envs/sourmash.yml"
    params:
        p=param_str,
        tmp = "data/metagenome_analyses/strain_tracking/sigs/{sample}.fastq.gz"
    shell: """
        cat {input.r1} {input.r2}> {params.tmp}
        sourmash sketch dna -p {params.p}\
           {params.tmp} -o {output}
        rm {params.tmp}
    """


rule sourmash_gather:
    input:
        sig = ancient("data/metagenome_analyses/strain_tracking/sigs/{sample}.abundtrim.sig"),
        db = ancient(output_dir_saureus+"genomes/final_drep/dereplicated_genomes/genome.sig")
    output:
        csv = output_dir_saureus + "genbank/{sample}.x.genbank.gather.csv",
        out = output_dir_saureus + "genbank/{sample}.x.genbank.gather.out",
    conda: "envs/sourmash.yml"
    params:
        threshold_bp = SOURMASH_DATABASE_THRESHOLD_BP,
        name = "data/metagenome_analyses/strain_tracking/sigs/{sample}.fastq.gz"
    shell: """
        sourmash gather {input.sig} {input.db} -o {output.csv} \
          --threshold-bp {params.threshold_bp}> {output.out} || true
        touch {output}
        if read -r && read -r
        then
        echo "This has more than 1 line."
        else
        echo 'intersect_bp,f_orig_query,f_match,f_unique_to_query,f_unique_weighted,average_abund,median_abund,std_abund,name,filename,md5,f_match_orig,unique_intersect_bp,gather_result_rank,remaining_bp,query_filename,query_name,query_md5,query_bp'>{output.csv}
        echo '0,0,0,0,0,0,0,0,NA,NA,NA,0,0,0,0,{params.name},NA,NA,NA'>>{output.csv}
        fi < {output.csv}
    """


rule sourmash_concatenate:
    input: expand(output_dir_saureus + "genbank/{sample}.x.genbank.gather.csv", sample=metagenomes)
    output: output_dir_saureus + "summary/genbank_all.csv"
    params:
       outdir=output_dir_saureus + "genbank/"
    shell:"""
    echo 'intersect_bp,f_orig_query,f_match,f_unique_to_query,f_unique_weighted,average_abund,median_abund,std_abund,name,filename,md5,f_match_orig,unique_intersect_bp,gather_result_rank,remaining_bp,query_filename,query_name,query_md5,query_bp'>{output}
   awk FNR!=1 {params.outdir}/*.genbank.gather.csv >> {output}
"""

