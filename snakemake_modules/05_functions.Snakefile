# vim::w set ft=python:
import pandas as pd

all_genomes=pd.read_csv("metadata_isolates.csv")
all_ids=all_genomes["run_accession"]


outdir="out_functions/"


rule all:
   input: expand(outdir+"pathways/{id}.summary.kegg_contigs.tsv",id=all_ids)



rule prokka:
    input:ancient("data/genomes/renamed/{isolate}.fasta")
    output:
        outdir+"data/genomes_analysis/annotate/prokka/{isolate}/{isolate}.faa",
        outdir+"data/genomes_analysis/annotate/prokka/{isolate}/{isolate}.gbk",
    singularity:
        "docker://staphb/prokka:1.14.5"
    params:
        out_prokka="data/genomes_analysis/annotate/prokka/{isolate}",
        prefix="{isolate}"
    shell:
        """
        prokka {input} --kingdom Bacteria --outdir {params.out_prokka} \
        --prefix {params.prefix} --force --locustag {params.prefix} --cpus {threads}
        """


rule run_hmm:
   input: 
      faa="data/genomes_analysis/annotate/prokka/{id}/{id}.faa",
   output: outdir+"hmmscan/{id}.txt"
   params:
      db="tools/PG_function-KEGG/ref-dbs/db_kofam.hmm",
   threads: 15
   shell:"""
   hmmscan --noali --cut_ga --domtblout {output} {params.db} {input.faa} --cpu {threads}
"""


rule parse:
  input: 
     faa="data/genomes_analysis/annotate/prokka/{id}/{id}.faa",
     hmm=outdir+"hmmscan/{id}.txt"
  output:
     tab=outdir+"hmmscan_tab/{id}.txt",
     parsed=outdir+"hmmscan_parsed/{id}.txt"
  group: "fast"
  params:
     out="{id}.txt_parsed",
  shell:"""
  python3 scripts/hmmscan_tab.py -i {input.hmm} -o {output.tab}
  python3 scripts/parsing_hmmscan.py -i {output.tab} -f {input.faa}
  mv {params.out} {output.parsed}
"""


rule give_pathways:
  input: outdir+"hmmscan_parsed/{id}.txt"
  output: outdir+"pathways/{id}.summary.kegg_contigs.tsv"
  params:
    out=outdir+"pathways/{id}" 
  group: "fast"
  shell:"""
  python3 tools/PG_function-KEGG/Tools/give_pathways.py -i {input} -g tools/PG_function-KEGG/help_files/graphs.pkl -c tools/PG_function-KEGG/help_files/all_pathways_class.txt -n tools/PG_function-KEGG/help_files/all_pathways_names.txt -o {params.out}
"""







