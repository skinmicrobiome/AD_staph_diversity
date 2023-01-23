#vim:w set ft=python:
import pandas as pd

DBDIR="data/databases/virus/virsearch_db/"
genomes_fna_dir="data/genomes/renamed/"


all_genomes=pd.read_csv("metadata_isolates.csv",sep=",")
all_genomes=all_genomes["run_accession"]

virus_dir="data/genomes_analysis/virus/"


#all_genomes,=glob_wildcards(virus_dir+"{id}/cdhit/nr_predictions.fa")
all_genomes,=glob_wildcards(virus_dir+"{id}/checkv/viral_sequences.fa")

extension="fasta"
rule all:
    input: 
     expand(virus_dir+"all_viruses/{id}.fa",id=all_genomes), 
     expand(virus_dir+"{id}/checkv/checkv_output/quality_summary.tsv",id=all_genomes)
     expand(virus_dir+"{id}/cdhit/nr_predictions.fa"
     expand(virus_dir+"{id}/demovir/DemoVir_assignments.txt",id=all_genomes)
     virus_dir+"all_viruses/cdhit/eggnog/done"
 



rule rename_virus:
    input: ancient(genomes_fna_dir+"{id}."+extension)
    output: virus_dir+"renamed/{id}.fasta"
    params:
        name="{id}"
    group: "fast"
    shell:
        """
        scripts/rename_multifasta_prefix_virus.py -f {input} -p {params.name} > {output}
        """


rule vibrant_detect:
    input:
        ancient(virus_dir+"renamed/{id}.fasta")
    output:
        virus_dir+"{id}/vibrant/VIBRANT_{id}/VIBRANT_phages_{id}/{id}.phages_combined.fna"
    params:
        outdir = virus_dir+"{id}/vibrant",
        database = DBDIR+"databases"
    conda:
        "envs/vibrant.yml"
    threads: 1
    shell:"""
        rm -rf {params.outdir}
        VIBRANT_run.py -t 4 -d {params.database} -i {input} -folder {params.outdir} -l 10000 -no_plot #|| mkdir $(dirname {output}) && touch {output}
"""


rule vibrant_rename:
    input:
        ancient(virus_dir+"{id}/vibrant/VIBRANT_{id}/VIBRANT_phages_{id}/{id}.phages_combined.fna")
    output:
        virus_dir+"{id}/vibrant/viral_sequences.fa"
    conda:
        "envs/vibrant.yml"
    group: "fast"
    threads: 1
    shell:
        "scripts/rename_multifasta_prefix.py -f {input} -p VIBRANT > {output}"



rule vs2_detect:
    input:
        ancient(virus_dir+"renamed/{id}.fasta")
    output:
        virus_dir+"{id}/virsorter2/final-viral-combined.fa"
    params:
        outdir = virus_dir+"{id}/virsorter2",
        database = "data/databases/virus/vs2/" #data/databases/virus/db/" 
    singularity: "docker://staphb/virsorter2:2.1" #
    threads: 1
    shell:
        "virsorter run -j 8 -i {input} --db-dir {params.database} -w {params.outdir} --min-length 10000 all"


rule vs2_rename:
    input:
        virus_dir+"{id}/virsorter2/final-viral-combined.fa"
    output:
        virus_dir+"{id}/virsorter2/viral_sequences.fa"
    group: "fast"
    threads: 1
    shell:
        "scripts/rename_multifasta_prefix.py -f {input} -p VS2 > {output}"



rule dvf_detect:
    input:
        virus_dir+"renamed/{id}.fasta"
    output:
        virus_dir+"{id}/deepvirfinder/{id}.fasta_gt10000bp_dvfpred.txt"
    params:
        outdir = virus_dir+"{id}/deepvirfinder",
        database = DBDIR+"models"
    conda:
        "envs/deepvirfinder.yml"
    threads: 1
    shell:
        "tools/DeepVirFinder/dvf.py -i {input} -m {params.database} -o {params.outdir} -l 10000 -c 4"


rule dvf_filter:
    input:
        dvf = virus_dir+"{id}/deepvirfinder/{id}.fasta_gt10000bp_dvfpred.txt",
        fa = virus_dir+"renamed/{id}.fasta"
    output:
        virus_dir+"{id}/deepvirfinder/viral_sequences.fa"
    params:
        contigs = virus_dir+"{id}/deepvirfinder/viral_contigs.txt",
        fa_ori = virus_dir+"{id}/deepvirfinder/viral_contigs.fa",
    conda:
        "envs/deepvirfinder.yml"
    threads: 1
    shell:
        """
        awk '{{if($3 > 0.9 && $4 < 0.01)print$1}}' {input.dvf} > {params.contigs}
        scripts/select_seqs_by_IDs.py -i {input.fa} -d {params.contigs} -o {params.fa_ori}
        scripts/rename_multifasta_prefix.py -f {params.fa_ori} -p DVF > {output}
        rm {params.fa_ori} {params.contigs}
        """


rule cat_input: 
    input:
        vs2 = virus_dir+"{id}/virsorter2/viral_sequences.fa",
        dvf = virus_dir+"{id}/deepvirfinder/viral_sequences.fa",
        vibrant = virus_dir+"{id}/vibrant/viral_sequences.fa"
    output:
        virus_dir+"{id}/checkv/viral_sequences.fa2"
    threads: 1
    group: "fast"
    priority: 10000
    shell:
        "cat {input.vs2} {input.dvf} {input.vibrant} > {output}"



rule checkv_analysis:
    input:
        virus_dir+"{id}/checkv/viral_sequences.fa"
    output:
        virus_dir+"{id}/checkv/checkv_output/quality_summary.tsv"
    params:
        outdir = virus_dir+"{id}/checkv/checkv_output/",
        database = DBDIR+"checkv/checkv-db-v1.0"
    conda:
        "envs/checkv.yml"
    threads: 1
    shell:
        "checkv end_to_end -t 8 -d {params.database} {input} {params.outdir}"


rule checkv_filter:
    input:
        virus_dir+"{id}/checkv/checkv_output/quality_summary.tsv"
    output:
        virus_dir+"{id}/checkv/checkv_output/filtered_predictions.fa"
    params:
        indir = virus_dir+"{id}/checkv/checkv_output/"
    threads: 1
    conda:
        "envs/checkv.yml"
    shell:
        """
        scripts/filter_checkv.py {params.indir}proviruses.fna {input} {params.indir}proviruses_filt.fna pro
        scripts/filter_checkv.py {params.indir}viruses.fna {input} {params.indir}viruses_filt.fna vir
        cat {params.indir}proviruses_filt.fna {params.indir}viruses_filt.fna > {output}
        """


rule checkv_filter2:
    input: virus_dir+"{id}/cdhit/nr_predictions.fa"
    output: virus_dir+"all_viruses/{id}.fa"
    params:
        name="{id}"
    shell:"""
    scripts/rename_multifasta_prefix.py -f {input} -p {params.name} > {output}
"""


# cluster sequences
rule dereplicate:
    input:
        virus_dir+"{id}/checkv/checkv_output/filtered_predictions.fa"
    output:
        pred = virus_dir+"{id}/cdhit/nr_predictions.fa"
    params:
        clst_dir = directory(virus_dir+"{id}/cdhit")
    threads: 1
    shell:
        """
        tools/cdhit/cd-hit-est -c 0.99 -i {input} -o {params.clst_dir}/filtered_predictions_nr -T 0 -M 0 -d 0
        mv {params.clst_dir}/filtered_predictions_nr {output.pred}
        """


rule cdhit:
    input: ancient(virus_dir+"all_viruses.fasta")
    output:
        pred = virus_dir+"all_viruses/cdhit/nr_predictions.fa",
        clstr = virus_dir+"all_viruses/cdhit/filtered_predictions_nr.clstr"
    params:
        clst_dir = directory(virus_dir+"all_viruses/cdhit/")
    threads: 1
    shell:
        """
        tools/cdhit/cd-hit-est -i {input} -o {params.clst_dir}/filtered_predictions_nr -c 0.90 -aS 0.75 -G 0 -T 0 -M 50000
        mv {params.clst_dir}/filtered_predictions_nr {output.pred}
        """


rule parse_cdhit_cluster:
    input: virus_dir+"all_viruses/cdhit/filtered_predictions_nr.clstr"
    output: virus_dir+"all_viruses/cdhit/filtered_predictions_nr_clstr.txt"
    threads: 1
    shell:"""
    perl tools/cdhit/clstr2txt.pl {input}>{output}
"""


rule run_prodigal_1:
   input: virus_dir+"all_viruses/cdhit/nr_predictions.fa",
   output: virus_dir+"all_viruses/cdhit/prodigal/nr_predictions.faa"
   threads: 1
   shell:"""
    prodigal -p 'meta' -a {output} -i {input}
"""

rule run_eggnog:
   input: virus_dir+"all_viruses/cdhit/prodigal/nr_predictions.faa"
   output: virus_dir+"all_viruses/cdhit/eggnog/done"
   params:
      outdir=virus_dir+"all_viruses/cdhit/eggnog/",
      tmpdir=virus_dir+"all_viruses/cdhit/eggnog/tmp/",
      MAG="nr_predictions"
   conda: "envs/py27.yaml"
   shell:"""
mkdir -p {params.outdir}
mkdir -p {params.tmpdir}
/data/sahebkashafs2/SOFTWARE/eggnog-mapper/emapper.py --cpu {threads} -i {input} -m diamond -o {params.MAG} --output_dir {params.outdir} --temp_dir {params.tmpdir}
touch {output}
"""



rule tax_class:
    input:
        virus_dir+"{id}/cdhit/nr_predictions.fa"
    output:
        tax = virus_dir+"{id}/demovir/DemoVir_assignments.txt",
        contig_ids = virus_dir+"{id}/demovir/trembl_ublast.viral.u.contigID.txt"
    params:
        demovir_dir = directory(virus_dir+"{id}/demovir"),
        database = DBDIR+"/demovir"
    conda:
        "envs/checkv.yml"
    threads: 1
    shell:
        """
        prodigal -a {params.demovir_dir}/proteins.faa -i {input} -p meta &> /dev/null
        tools/usearch -ublast {params.demovir_dir}/proteins.faa -db {params.database}/uniprot_trembl.viral.udb -evalue 1e-5 -trunclabels -blast6out {params.demovir_dir}/trembl_ublast.viral.txt -threads 4 &> /dev/null
        sort -u -k1,1 {params.demovir_dir}/trembl_ublast.viral.txt > {params.demovir_dir}/trembl_ublast.viral.u.txt
        cut -f 1,2 {params.demovir_dir}/trembl_ublast.viral.u.txt | sed 's/_[0-9]\+\t/\t/' | cut -f 1 | paste {params.demovir_dir}/trembl_ublast.viral.u.txt - > {output.contig_ids}
        scripts/demovir.R {params.demovir_dir} {params.database}
        """


