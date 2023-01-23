# -*- coding: utf-8
import pandas as pd


#metagenome dir
metagenomesdir="data/00_preprocessing/processed/metagenomes/AD/" 

#define runs
df_run=pd.read_csv("metadata_all_metagenomes.tsv",sep="\t")
#df_run=df_run[df_run["Study"]=="To be announced"]
metagenomes=df_run["Sample_ID"]

drep="0.998"

output_dir_sepidermidis="data/metagenome_analyses/strain_tracking/sepidermidis_NIH_JOh_genomes_"+drep+"/"

mlst_dir="data/genomes_analysis/mlst/staphepi/"


genomes_fna_dir="data/genomes/renamed/"

all_gen=pd.read_csv("metadata_files_final/sourmash_epi_genomes.txt")   
all_gen=all_gen["genome"]

rule gather_strains:
    input: 
        expand(output_dir_sepidermidis + "all_summary/"+drep+"_{kmer}_{scale}_{thres}__genbank_all.csv",thres="1e5",kmer="61",scale="10000"),expand(output_dir_sepidermidis + "{kmer}_{scale}_{thres}/genbank/{sample}.x.genbank.gather.csv", sample=metagenomes,thres="1e5",kmer="61",scale="10000")
       # expand(output_dir_sepidermidis + "all_summary/"+drep+"_{kmer}_{scale}_{thres}__genbank_all.csv",thres="1e6",kmer="61",scale="10000"),
      #  expand(output_dir_sepidermidis + "all_summary/"+drep+"_{kmer}_{scale}_{thres}__genbank_all.csv",thres="1e5",kmer="51",scale="10000"),
     #   expand(output_dir_sepidermidis + "all_summary/"+drep+"_{kmer}_{scale}_{thres}__genbank_all.csv",thres="1e6",kmer="51",scale="10000"),
    #    expand(output_dir_sepidermidis + "all_summary/"+drep+"_{kmer}_{scale}_{thres}__genbank_all.csv",thres="1e5",kmer="31",scale="10000"),
   #     expand(output_dir_sepidermidis + "all_summary/"+drep+"_{kmer}_{scale}_{thres}__genbank_all.csv",thres="1e6",kmer="31",scale="10000"),
     #   expand(output_dir_sepidermidis + "all_summary/"+drep+"_{kmer}_{scale}_{thres}__genbank_all.csv",thres="1e5",kmer="61",scale="1000"),
  #      expand(output_dir_sepidermidis + "all_summary/"+drep+"_{kmer}_{scale}_{thres}__genbank_all.csv",thres="1e5",kmer="51",scale="1000"),
     #   expand(output_dir_sepidermidis + "all_summary/"+drep+"_{kmer}_{scale}_{thres}__genbank_all.csv",thres="1e5",kmer="61",scale="100000"),
 #       expand(output_dir_sepidermidis + "all_summary/"+drep+"_{kmer}_{scale}_{thres}__genbank_all.csv",thres="1e5",kmer="51",scale="100000"),



#start off with these skin genomes
rule refmt_mlst:
    input:
       met="metadata_files/specifications/sourmash_epi_genomes.txt",
       isolates=ancient(expand(mlst_dir+"srst2/{id}__mlst__Staphylococcus_epidermidis__results.txt", id=all_gen)),
       mlst="data/genomes_analysis/mlst/staphepi/summary/mlst_summ.txt"
    output: output_dir_sepidermidis+"genomes/all_ST.txt"
    shell:"""
    rm -rf {output}
    while read -r i; do grep $i {input.mlst} >>{output}; done < {input.met}
    sed -i 's/\?//g' {output}
    sed -i 's/\*//g' {output}
"""



rule subset_ST:
   input: output_dir_sepidermidis+"genomes/all_ST.txt"
   output: output_dir_sepidermidis+"genomes/subset_ST.txt"
   params:
     dir=output_dir_sepidermidis+"genomes/subset_ST",
     fna=genomes_fna_dir,
     tmp=output_dir_sepidermidis+"genomes/subset_ST_ids.txt" 
   shell:"""
   rm -rf {params.dir}
   mkdir -p {params.dir}
   Rscript --vanilla scripts/subset_STs_strains.R {input} {output}
   cut -f1 -d , {output} > {params.tmp}
   while read -r line; do scp {params.fna}/${{line}}.fasta {params.dir}; done<{params.tmp}
"""



rule dereplicate_genome_db:
    input: output_dir_sepidermidis+"genomes/subset_ST.txt"
    output: output_dir_sepidermidis+"genomes/final_drep/done"
    params:
      indir=output_dir_sepidermidis+"genomes/subset_ST/",
      outdir=output_dir_sepidermidis+"genomes/final_drep/",
      drep=drep
    singularity:
        "shub://sskashaf/MAG_wf_containers_2021:drep"
    shell:"""
   rm -rf {params.outdir} 
   #gunzip {params.indir}*.fna.gz
   dRep dereplicate -p {threads} {params.outdir} -g {params.indir}*.fasta -pa 0.95 -sa {params.drep} -nc 0.6 -cm larger --ignoreGenomeQuality
   touch {output}
"""


#sketch genome database
rule sketch_genomes:
    input: output_dir_sepidermidis+"genomes/final_drep/done",
    output: output_dir_sepidermidis+"genomes_{kmer}_{scale}/genome.sig" 
    params:
      dir=output_dir_sepidermidis+"genomes/final_drep/dereplicated_genomes/",
      outdir=output_dir_sepidermidis+"genomes_{kmer}_{scale}/",
      scale="{scale}",
      kmer="{kmer}"
    shell:"""
rm -r {params.outdir}
scp -r {params.dir} {params.outdir}
cd {params.outdir}
sourmash sketch dna -p k={params.kmer},scaled={params.scale},abund * -o genome.sig
"""


# compute sourmash signature for reads
rule sketch_reads:
    input:
        r1 = metagenomesdir+"{sample}_1.fastq.gz",
        r2 = metagenomesdir+"{sample}_2.fastq.gz"
    output:
        sig = "data/metagenome_analyses/strain_tracking/sigs_{kmer}_{scale}/{sample}.abundtrim.sig"
    conda: "envs/sourmash.yml"
    params:
        scale="{scale}",
        kmer="{kmer}",
        tmp = "data/metagenome_analyses/strain_tracking/sigs_{kmer}_{scale}/{sample}.fastq.gz"
    shell: """
        cat {input.r1} {input.r2}> {params.tmp}
        sourmash sketch dna -p k={params.kmer},scaled={params.scale}\
           {params.tmp} -o {output}
        rm {params.tmp}
    """


rule sourmash_gather:
    input:
        sig = ancient("data/metagenome_analyses/strain_tracking/sigs_{kmer}_{scale}/{sample}.abundtrim.sig"),
        db = output_dir_sepidermidis+"genomes_{kmer}_{scale}/genome.sig"
    output:
        csv = output_dir_sepidermidis + "{kmer}_{scale}_{thres}/genbank/{sample}.x.genbank.gather.csv",
        out = output_dir_sepidermidis + "{kmer}_{scale}_{thres}/genbank/{sample}.x.genbank.gather.out",
    conda: "envs/sourmash.yml"
    params:
        threshold_bp = "{thres}",
        name = "data/metagenome_analyses/strain_tracking/sigs_{kmer}_{scale}/{sample}.fastq.gz"
    group: "fast"
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
    input: expand(output_dir_sepidermidis + "{{kmer}}_{{scale}}_{{thres}}/genbank/{sample}.x.genbank.gather.csv", sample=metagenomes)
    output: output_dir_sepidermidis + "{kmer}_{scale}_{thres}/summary/genbank_all.csv"
    params:
      out=output_dir_sepidermidis + "{kmer}_{scale}_{thres}/genbank/",
      sum=output_dir_sepidermidis + "{kmer}_{scale}_{thres}/summary/"
    shell:"""
    mkdir -p {params.sum}
    rm -rf {output}
    echo 'intersect_bp,f_orig_query,f_match,f_unique_to_query,f_unique_weighted,average_abund,median_abund,std_abund,name,filename,md5,f_match_orig,unique_intersect_bp,gather_result_rank,remaining_bp,query_filename,query_name,query_md5,query_bp'>{output}
    cd {params.out}
    awk FNR!=1 *genbank.gather.csv >> ../summary/genbank_all.csv
"""



rule group_output:
    input: output_dir_sepidermidis + "{kmer}_{scale}_{thres}/summary/genbank_all.csv"
    output: output_dir_sepidermidis + "all_summary/"+drep+"_{kmer}_{scale}_{thres}__genbank_all.csv"
    shell:"""
    scp {input} {output}
"""

