# vim: ft=python
import sys
import os
import csv
import os.path
from os import path

shell.prefix("set -eo pipefail; ")
configfile:"config/config.yaml"
localrules: all

run = os.getcwd() + "/"
raw = config["raw"]
fastq = config["fastq"]
out = config["out"]

snRNA = config["snRNA"]
flowcell=config["flowcells"]
table = config["table"]
medalt = config["medalt"]

genome = config["genome"]
fasta = config["fasta"]
gtf = config["gtf"]
genes = config["genes"]
scRNA_ref = config["scRNA_ref"]

# cytoscape style files
scCNV_group = run + "cytoscape_styles/scCNV_group.xml"

with open(table,'r') as f:
     reader = csv.reader(f)
     sample = list(zip(*reader))[1]
     sample = sample[1:]

if config["patient"]=="ready":
   patient = open('patient.txt','r').read().splitlines()
else:
    patient=""

wildcard_constraints:
     sample= '|'.join([re.escape(x) for x in sample]),
     patien= '|'.join([re.escape(x) for x in patient])

if config["patient"]=="ready":
   include: "modules/Snakefile_cellranger"
   include: "modules/Snakefile_gen"
   include: "modules/Snakefile_patient"
   rule all:
       input:
             expand(out + "link/counts/{sample}_filtered_feature_bc_matrix.csv",sample=sample),
             expand(out + "reanalysis/link/loup/{sample}_cloupe.cloupe",sample=sample),
             expand(out + "reanalysis/link/cogaps/{sample}_cogaps.html",sample=sample),
             expand(out + "reanalysis/{sample}/outs/conicsmat/conicsmat.html",sample=sample),
             expand(out + "patient/{patient}/MEDALT_group/medalt.group.force.directed.pdf",patient=patient)

else:
     include: "modules/Snakefile_cellranger"
     include: "modules/Snakefile_gen"
     rule all:
         input:
               expand(out + "link/counts/{sample}_filtered_feature_bc_matrix.csv",sample=sample),
               expand(out + "reanalysis/link/loup/{sample}_cloupe.cloupe",sample=sample),
               expand(out + "reanalysis/link/cogaps/{sample}_cogaps.html",sample=sample),
               expand(out + "reanalysis/{sample}/outs/conicsmat/conicsmat.html",sample=sample)
