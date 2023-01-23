#vim: set ft=python:


metadata_folder="snakemake_modules/gwas_metadata_saureus/"
output_folder="data/genomes_analysis/GWAS/saureus/GWAS_aur_final_paper_all_vf/"
tree="data/genomes_analysis/tree/saureus_all/results_paperfinal/gubbins/gubbins.final_tree.tre" 
cpus=40
GWAS_dir="data/genomes_analysis/GWAS_2022/"
dimensions="4"

phenos_ratio=["scorad","scorad_cate"] 





rule all:
   input: 
            expand(output_folder+"gwas_unitig_{pheno_column}/unitig_significance_annotated.txt",pheno_column=phenos_ratio),
            expand(output_folder+"gwas_unitig_{pheno_column}/gene_hits.txt", pheno_column=phenos_ratio),
            expand(output_folder+"gwas_unitig_{pheno_column}/qq_plot_{pheno_column}.png",pheno_column=phenos_ratio),
            output_folder+"eggnog/pangenome_done",
            expand(output_folder+"gwas_unitig_{pheno_column}/unitigs_position.txt",pheno_column=phenos_ratio)




#similarity matrix
rule similarity_matrix:
    input: tree
    output: output_folder+"gubbins/similarity_matrix.txt",
    shell:
        "python tools/pyseer/scripts/phylogeny_distance.py --lmm {input} > {output}"



#count unitigs
rule unitigs:
    input: metadata_folder+"strains.txt"
    output: output_folder+"unitigs/unitigs.txt"
    params:
        outdir=output_folder+"unitigs/"
    threads: cpus
    shell:
        """
        rm -rf {params.outdir}
        unitig-counter -strains {input} -output {params.outdir} -nb-cores {threads}
        """


#copy reference
rule copy_ref:
  input:"data/genomes/renamed/aatyw.fasta" #data/genomes/reference/saureus/GCF_000013465.1_ASM1346v1_genomic.fna"
  output:
      fna=output_folder+"reference/ref.fna",
  shell:"""
  scp {input} {output.fna}
  """


#reference
rule rename_fa:
    input: output_folder+"reference/ref.fna"
    output: output_folder+"reference/ref_renamed.fna"
    params:
        name="ref"
    group: "fast"
    shell:
        """
        scripts/rename_multifasta_prefix.py -f {input} -p {params.name} > {output}
        """



#run prokka on reference
rule prokka:
    input: output_folder+"reference/ref_renamed.fna"
    output:
        output_folder+"reference/prokka/ref.gff",
        output_folder+"reference/prokka/ref.faa",
    #singularity:
     #   "docker://staphb/prokka:1.14.5"
    params:
        out_prokka=output_folder+"reference/prokka/",
        prefix="ref"
    threads: cpus
    shell:
        """
        module load prokka
        prokka {input} --kingdom Bacteria --outdir {params.out_prokka} \
        --prefix {params.prefix} --force --locustag {params.prefix} --cpus {threads}
        """



# make ref file
rule make_ref:
   input:
      strains=metadata_folder+"strains.txt",
      fna=output_folder+"reference/ref_renamed.fna",
      gff=output_folder+"reference/prokka/ref.gff",
   output: output_folder+"reference/references.txt"
   params:
      name="ref",
      fnadir="data/genomes/renamed",
      annotate="data/genomes_analysis/annotate/prokka",
      strains=output_folder+"reference/references.txt_tmp"
   shell:"""
   rm -rf {output}
   echo -e "{input.fna}\t{input.gff}\t{params.name}">>{output}
"""


rule make_ref_multiple:
   input:
      strains=metadata_folder+"strains.txt",
      fna=output_folder+"reference/ref_renamed.fna",
      gff=output_folder+"reference/prokka/ref.gff",
   output: output_folder+"reference/references_multiple.txt"
   params:
      name="ref",
      fnadir="data/genomes/renamed",
      annotate="data/genomes_analysis/annotate/prokka",
      strains=output_folder+"reference/references_multiple.txt_tmp"
   shell:"""
   rm -rf {output}
   echo -e "{input.fna}\t{input.gff}\t{params.name}">>{output}
   cut -f1 {input.strains} | grep -v 'run_accession' > {params.strains}
   while read -r line; do VAR1={params.fnadir}/${{line}}.fasta; VAR2={params.annotate}/${{line}}/${{line}}.gff; echo -e "$VAR1\t$VAR2\tdraft" >>{output};done<{params.strains}
   """


#resize
rule reduce_size:
  input:
    pheno=metadata_folder+"pheno.txt",
    m=output_folder+"gubbins/similarity_matrix.txt",
    i=metadata_folder+"pheno.txt",
    lineage=metadata_folder+"clusters.txt",
  output:
    pheno=metadata_folder+"pheno_{pheno_column}.txt",
    sim=output_folder+"gubbins/similarity_matrix_{pheno_column}.txt",
    lineage=metadata_folder+"clusters_{pheno_column}.txt"
  params:
    col="{pheno_column}"
  shell:"""
    module load R
    Rscript --vanilla scripts/resize_sim_matrix.R {input.i} {params.col} {input.m} {output.sim} {input.lineage} {output.lineage} {input.pheno} {output.pheno}
"""



#lmm
rule lmm_gwas_lineage:
    input:
        pheno=metadata_folder+"pheno_{pheno_column}.txt",
        unitigs=output_folder+"unitigs/unitigs.txt",
        sim=output_folder+"gubbins/similarity_matrix_{pheno_column}.txt",
        metadata=metadata_folder+"met.txt",
    output:
        patterns=output_folder+"gwas_unitig_{pheno_column}/unitig_patterns.txt",
        associations=output_folder+"gwas_unitig_{pheno_column}/unitig_significance.txt",
    params:
        col="{pheno_column}",
        dimensions=dimensions
    log:output_folder+"gwas_unitig_{pheno_column}/log"
    conda:
        "envs/pyseer.yml"
    threads: cpus
    shell:
        """
        pyseer --print-samples --phenotypes {input.pheno} --phenotype-column {params.col} --covariates {input.metadata} --use-covariates 2 --uncompressed --kmers {input.unitigs} --max-dimensions {params.dimensions} --output-patterns {output.patterns} --lmm --similarity {input.sim} > {output.associations} 2> {log}
        """


rule pyseer_count_patterns:
    input:output_folder+"gwas_unitig_{pheno_column}/unitig_patterns.txt",
    output: output_folder+"gwas_unitig_{pheno_column}/significance_limits.txt",
    params:
    resources:
        cpus=1,
        mem_mb=1000,
        time=10,
    conda:
        "envs/pyseer.yml"
    shell:
        """
        python tools/pyseer/scripts/count_patterns.py {input} > {output}
        """


rule filter_significant:
    input:
        script="scripts/filter_significant_unitigs.R",
        limit= output_folder+"gwas_unitig_{pheno_column}/significance_limits.txt",
        unitig_significance=output_folder+"gwas_unitig_{pheno_column}/unitig_significance.txt",
    output: output_folder+"gwas_unitig_{pheno_column}/unitig_significance_filtered.txt",
    shell:
        """
        module load R
        Rscript scripts/filter_significant_unitigs.R {input.limit} {input.unitig_significance} {output}
        """


rule phandango_input:
    input:
        unitigs=output_folder+"gwas_unitig_{pheno_column}/unitig_significance_filtered.txt",
        reference_genome=output_folder+"reference/ref_renamed.fna",
    output:output_folder+"gwas_unitig_{pheno_column}/unitigs_position.txt"
    shell:
        """
        phandango_mapper {input.unitigs} {input.reference_genome} {output}
        """


rule annotate:
    input:
        unitig_filtered=output_folder+"gwas_unitig_{pheno_column}/unitig_significance_filtered.txt",
        unitig=output_folder+"gwas_unitig_{pheno_column}/unitig_significance.txt",
        reference=output_folder+"reference/references.txt"
    output:
        annotated=output_folder+"gwas_unitig_{pheno_column}/unitig_significance_annotated.txt",
        annotated_filtered=output_folder+"gwas_unitig_{pheno_column}/unitig_significance_annotated_filtered.txt",
    params:
        out_dir=output_folder+"gwas_unitig_{pheno_column}/",
    shell:
        """
        annotate_hits_pyseer {input.unitig_filtered} {input.reference} {output.annotated_filtered}
        annotate_hits_pyseer {input.unitig} {input.reference} {output.annotated}
"""
   #    if [ -f remaining_kmers.txt ]; then
    #    mv remaining_kmers.fa {params.out_dir}/extended_remaining_kmers.fa
     #   mv remaining_kmers.txt {params.out_dir}/extended_remaining_kmers.txt
      #  fi
       # annotate_hits_pyseer {input.unitig} {input.reference} {output.annotated}
       # if [ -f remaining_kmers.txt ]; then
       # mv remaining_kmers.fa {params.out_dir}/extended_remaining_kmers.fa
       # mv remaining_kmers.txt {params.out_dir}/extended_remaining_kmers.txt
       # fi
       # """


rule annotate_multiple:
    input:
        unitig_filtered=output_folder+"gwas_unitig_{pheno_column}/unitig_significance_filtered.txt",
        unitig=output_folder+"gwas_unitig_{pheno_column}/unitig_significance.txt",
        reference=output_folder+"reference/references" #_multiple.txt"
    output:
        annotated=output_folder+"gwas_unitig_{pheno_column}/unitig_significance_annotated_multiple.txt",
        annotated_filtered=output_folder+"gwas_unitig_{pheno_column}/unitig_significance_annotated_filtered_multiple.txt",
    params:
        out_dir=output_folder+"gwas_unitig_{pheno_column}/",
    shell:
        """
        annotate_hits_pyseer {input.unitig_filtered} {input.reference} {output.annotated_filtered}
        annotate_hits_pyseer {input.unitig} {input.reference} {output.annotated}
"""


rule summarise:
    input: output_folder+"gwas_unitig_{pheno_column}/unitig_significance_annotated_filtered.txt", #_multiple.txt",
    output: output_folder+"gwas_unitig_{pheno_column}/gene_hits.txt"
    shell:"""
    python tools/pyseer/scripts/summarise_annotations.py --nearby {input} > {output}
"""


rule qqplot:
    input: output_folder+"gwas_unitig_{pheno_column}/unitig_significance.txt"
    output: output_folder+"gwas_unitig_{pheno_column}/qq_plot_{pheno_column}.png"
    params:
       indir=output_folder+"gwas_unitig_{pheno_column}/",
       infile="unitig_significance.txt",
       outfile="qq_plot_{pheno_column}.png"
    shell:"""
    cd {params.indir}
    python /data/sahebkashafs2/AD_clean/tools/pyseer/scripts/qq_plot.py {params.infile} > qq_plot.png
    mv qq_plot.png {params.outfile}
"""


rule run_eggnog_pan2:
   input: output_folder+"reference/prokka/ref.faa",
   output:  output_folder+"eggnog/pangenome_done"
   params:
      outdir=output_folder+"eggnog/",
      tmpdir=output_folder+"eggnog/tmp/",
      MAG="pangenome"
   conda: "envs/py27.yaml"
   threads: cpus
   shell:"""
mkdir -p {params.outdir}
mkdir -p {params.tmpdir}
/data/sahebkashafs2/SOFTWARE/eggnog-mapper/emapper.py --cpu 10 -i {input} -m diamond -o {params.MAG} --output_dir {params.outdir} --temp_dir {params.tmpdir}
touch {output}
"""



