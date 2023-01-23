# vim: set ft=python:

import pandas as pd

outdir="data/genomes_analysis/plasmids/"
genomes_dir="data/genomes/renamed/"
plasmids_dir="data/genomes_analysis/plasmids/mobsuite/all_plasmids/"

extension=".fasta"
all_plasmids,=glob_wildcards(plasmids_dir+"{plasmids}.fasta")


#genomes for plasmid analysis
df_genomes=pd.read_csv("metadata_isolates.csv",sep=",")
all_genomes=df_genomes["run_accession"]


#genomes to compare via dnadiff
all_genomes_compare=pd.read_csv("metadata_files/2080_plasmids_all.txt")
all_genomes_compare=all_genomes_compare["sample_id"]


#genomes where plasmids were detected
mob,=glob_wildcards(outdir+"mobsuite/{gen}/mobtyper_results.txt")
eggnog,=glob_wildcards(outdir+"platon/{gen}/eggnog/pangenome.emapper.annotations")


rule all:
   input: expand(outdir+"mobsuite/{gen}/mobtyper_results_renamed.txt",gen=all_genomes), 
           expand(outdir+"plasmid_compare/sortmashhit/logs/{bins}_done.txt",bins=all_genomes_compare),
           expand("data/genomes_analysis/plasmids/mobsuite/plasmid_amr/{plas}.txt",plas=all_plasmids)
           #outdir+"plasmid_compare/dnadiff_summary.tsv"


###############Reconstruct the plasmids#######################################

rule recon_plasmid:
    input: ancient(genomes_dir+"{gen}"+extension)
    output: outdir+"mobsuite/{gen}/mobtyper_results.txt", outdir+"mobsuite/{gen}/contig_report.txt"
    singularity:
        "docker://kbessonov/mob_suite:3.0.3"
    params:
      outdir=outdir+"mobsuite/{gen}/",
      gen="{gen}"
    threads: 5
    shell:"""
    mob_recon -i {input} -o {params.outdir} -n {threads} --force
    cd {params.outdir}
    for i in *fasta; do mv "${{i}}" {params.gen}_"${{i}}"; done
"""


rule recon_parse:
   input: outdir+"mobsuite/{gen}/mobtyper_results.txt"
   output: outdir+"mobsuite/{gen}/mobtyper_results_renamed.txt"
   shell:"""
  awk '{{print FILENAME (NF?"\\t":"") $0}}' {input}>{output}
"""


rule aggregate_recon:
   input: expand(outdir+"mobsuite/{gen}/mobtyper_results_renamed.txt",gen=mob)
   output: outdir+"all_mobsuite_summ.txt"
   shell:"""
   awk FNR!=1 {input} >> {output}
"""





###############Compare to COMPASS database###################################

rule make_mash_db:
    input: "data/databases/COMPASS/data/COMPASS.fasta"
    output: outdir+"ref.msh"
    singularity:
        "docker://quay.io/biocontainers/mash:2.2.1--h3d38be6_0"
    group:"fast"
    params:
      indir="data/databases/COMPASS/data/indiv/"
    shell:"""
    mash sketch -o {output} {params.indir}/*fasta
"""



rule mash_dist:
    input:
        bins=plasmids_dir+"{i}",
        db=outdir+"ref.msh",
    output: outdir+"plasmid_compare/mashdist/{i}.tab"
    threads:1
    singularity:"docker://quay.io/biocontainers/mash:2.2.1--h3d38be6_0"
    group:"fast"
    shell:"""
    mash dist -p {threads} {input.db} {input.bins} > {output}
"""


rule best_mash:
    input:
        mashdist=outdir+"plasmid_compare/mashdist/{i}.tab"
    output: outdir+"plasmid_compare/best_mash/{i}.tab"
    threads:1
    singularity:
        "shub://sskashaf/Containers:isolatescompare"
    group:"fast"
    shell:
        """
        sort -gk3 {input.mashdist}|sed -n 1p >{output}
        """


rule sort_besthits:
    input: outdir+"plasmid_compare/mashdist/{bins}.tab"
    output: outdir+"plasmid_compare/sortmashhit/{bins}.tab"
    threads:1
    params:
      tmp= outdir+"plasmid_compare/sortmashhit/tmp_{bins}.tab"
    shell: """
       sort -gk3 {input}>{params.tmp}
       head -n 15 {params.tmp}> {output}
       rm {params.tmp}
"""


rule dna_diff:
      input: outdir+"plasmid_compare/sortmashhit/{bins}.tab"
      output: outdir+"plasmid_compare/sortmashhit/logs/{bins}_done.txt"
      params:
         indir="/data/sahebkashafs2/AD_clean/",
         outdir= outdir+"plasmid_compare/dnadiff_all/",
      #singularity: "docker://quay.io/biocontainers/mummer:3.23--pl526_7"
      priority: 1000
      shell: """
            mkdir -p {params.outdir}
            cd {params.outdir}
            while read col1 col2 rem
            do
               echo 'dnadiff ${{col2}} ${{col1}} -p ${{col1%%.fa}}_${{col2%%.fasta}}_'
               one=$(basename ${{col2}} .fa)
               two=$(basename ${{col1}} .fasta)
               dnadiff {params.indir}${{col1}} {params.indir}${{col2}} -p ${{one}}_${{two}}_
            done < {params.indir}{input}
            cd {params.indir}
            touch {output}
            """

rule parse_dnadiff:
    input:
        dnadiff=outdir+"plasmid_compare/dnadiff_all/{i}.report"
    group:"fast"
    output:dnadiff=outdir+"plasmid_compare/dnadiff_all_parsed/{i}_parsed.tsv"
    run:
        outfile = str(output)
        f = open(input.dnadiff)
        data = f.read()
        first_line = data.split("\n", 1)[0]
        a = first_line.split(" ")
        ref = a[0]
        quer = a[1]
        with open(outfile, "w") as outf:
            path_dna = input.dnadiff
            base = os.path.basename(path_dna)
            base = base.split(".report")[0]
            with open(path_dna) as f:
                for line in f:
                    if "TotalBases" in line:
                        cols = line.split()
                        lenref = int(cols[1])
                        lenquer = int(cols[2])
                    if "AlignedBases" in line:
                        cols = line.split()
                        aliref = cols[1].split("(")[-1].split("%")[0]
                        alique = cols[2].split("(")[-1].split("%")[0]
                    if "AvgIdentity" in line:
                        cols = line.split()
                        ident = float(cols[1])
            line = "%s\t%s\t%i\t%.2f\t%i\t%.2f\t%.2f" % (ref, quer, lenref, float(aliref), lenquer, float(alique), float(ident))
            outf.writelines(line + "\n")


rule aggregate_dnadiff:
    input: expand(outdir+"plasmid_compare/dnadiff_all_parsed/{i}_parsed.tsv",i=all_genomes_compare)
    output: outdir+"plasmid_compare/dnadiff_summary.tsv"
    shell:
        """
        cat {input}>{output}
        """


###############AMR###################################

rule amr:
   input:plasmids_dir+"{id}.fasta"
   output:"data/genomes_analysis/plasmids/mobsuite/plasmid_amr/{id}.txt"
   shell: """
amrfinder -n {input} > {output}
"""
