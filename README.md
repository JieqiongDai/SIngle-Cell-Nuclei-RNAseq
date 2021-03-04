# Single Cell/Nuclei RNA-sequencing Analysis
## Description
This snakemake pipeline is for single cell/nuclei RNA-sequencing analysis using 10 X single cell RNA-seq generated data. The pipeline may be run on an HPC or in a local environment.

Major steps in the workflow include:
1) Primary QC and general analysis using cellranger
2) Secondary QC ([DIEM](https://github.com/marcalva/diem)) and reanalysis
3) Futher QC filtering and generaly analysis using [Seurat](https://satijalab.org/seurat/)
4) Reference-based cell type annotation using [singleR](https://bioconductor.org/packages/release/bioc/html/SingleR.html)
5) Coordinated gene association analysis using [CoGAPS](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-020-03796-9)
6) Copy number variation (CNV) analysis using [CONICSmat](https://github.com/diazlab/CONICS/wiki/Tutorial---CONICSmat;---Dataset:-SmartSeq2-scRNA-seq-of-Oligodendroglioma)
7) CNV evolutionary analysis at cell cluster leveland result visualization using [infercnv (https://www.bioconductor.org/packages/release/bioc/html/infercnv.html), [MEDALT](https://github.com/KChen-lab/MEDALT) and [Cytoscape](https://cytoscape.org/) in tumor-normal paired samples 

Expected results include:
* Cellular transcriptome profiles: 
  * Accessible to expression patterns of interested genes across cell clusters
  * Analysis of cell cluster related genes and pathways
  * Coordinated gene association detection
* Intratumor heterogeneity analysis in single-cell resolution:
  * Major cell clusters estimation, 
  * Rare subclone detection
  * Cell type identification
* Cellular level relative copy number profiles: 
  * Accessible to CNAs in the resolution of chromosome arm level
* Reconstructing tumor copy number evolution lineages (in tumor-normal paired samples): 
  * Lineage tracing of major tumor cell clusters 
  * Lineage tracing of CNA (deep diving)
  * Analysis of CNA related genes and pathways (deep diving)

## Software Requirements
* [Snakemake](https://snakemake.readthedocs.io/en/stable/)
* [cellranger](https://support.10xgenomics.com/single-cell-gene-expression/software/overview/welcome)
* [R](https://www.r-project.org)
* [MEDALT](https://github.com/KChen-lab/MEDALT)
* [Xvfb](https://www.x.org/releases/X11R7.6/doc/man/man1/Xvfb.1.xhtml)
* [Cytoscape](https://cytoscape.org/)

## Run modes
The pipeline has two run modes available; The general run mode is basic and the CNV evolutionary run mode is dependent on it; The detail of how-to-run is described in User's guider:
* General analysis: general single cell/nuclei RNA-seq analysis 
* CNV evolutionary analysis: CNV evolutionary analysis at cell cluster level in tumor-normal paired samples

## User's guide
### I. Input requirements
Basic:
* Edited config/config.yaml
* 10 X single cell CNV raw data
* 10 X simple sample sheet csv file
* [cellramger-dna reference](https://support.10xgenomics.com/single-cell-dna/software/pipelines/latest/advanced/references)

Optional:
* Required for the run mode of CNV evolutionary analysis:
  * {working_directory}/patient.txt
  * gene position file required by the tool infercnv

### II. 10 X simple sample sheet csv file format
Three columns with headers: Lane,Sample,Index

Example:
```bash
Lane,Sample,Index
1,A,SI-GA-A4
1,B,SI-GA-B4
2,A,SI-GA-A4
2,B,SI-GA-B4
```

### III. patient.txt file format
One headerless column: patient ID

Example:
```bash
SI-A
```
* Note, if original sample IDs of the patient related samples are not labeled as: {patientID}-N and {patientID}-T, create symbolic inks of the seurat data of original samples in the same directory with modified names as: {output_directory}/reanalysis/{patientID}-N/outs/seurat/seurat.Rdata and {output_directory}/reanalysis/{patientID}-T/outs/seurat/seurat.Rdata  

### IV. Editing the config.yaml
Basic parameters:
* snRNA: If single nuclei RNA-seq, yes or no
* flowcells: Flowcell ID
* raw: Path to the raw 10 X data stored directory
* table: Path to 10 X sample sheet
* fastq: Path to desired directory to store fastq files
* out: Path to desired directory to store output files
* medalt: Path to MEDALT package installed directory
* genome: hg38, hg19, mm10, etc
* fasta: Path to reference fasta file
* gtf: Path to reference gtf file
* scRNA_ref: Path to cellranger reference stored directory

Optional parameters:
* genes: Path to gene position file required for tool infercnv
* patient: Input 'ready' to initiate the run mode of CNV evolutionary analysis in tumor-normal paired samples when the basic run mode is complete, and the require patient.txt file and seurat data in the right naming format are ready

### V. To run
* Clone the repository to your working directory
```bash
git clone https://github.com/JieqiongDai/SingleCellCNV-Evolution.git
```
* Install required software; To run on NIH biowulf (an HPC using slurm job scheduler), you only need to download the MDEALT package and module load other required software.
* Edit and save config/config.yaml 
* To run on an HPC using slurm job scheduler: 
  Edit config/cluster_config.yaml according to your HPC information
  Run sbatch.sh to initiate running of the pipeline 
* To run in a local environment:
  ```bash
  snakemake -p --cores 14 --keep-going --rerun-incomplete --jobs 300 --latency-wait 120 all
  ```
* Look in log directory for logs for each rule
* To view the snakemkae rule graph:
```bash
snakemake --rulegraph | dot -T png > scDNA_whole.png
```
![dag](https://github.com/JieqiongDai/SingleCellCNV-Evolution/blob/master/scDNA_whole.png)


### V. Example output
```bash
. user/defined/output_dir/{}
├── link # main output files of cellranger-dna of all samples
│   ├── alarm # alarm files
│   │   ├── {sample_A}_alarms_summary.txt
│   │   ├── {sample_B}_summary.txt
│   │   └── {sample_C}_alarms_summary.txt 
│   ├── bam # indexed bam files
│   │   ├── {sample_A}_possorted_bam.bam
│   │   ├── {sample_A}_possorted_bam.bam.bai
│   │   ├── {sample_B}_possorted_bam.bam
│   │   ├── {sample_B}_possorted_bam.bam.bai
│   │   ├── {sample_C}_possorted_bam.bam
│   │   └── {sample_C}_possorted_bam.bam.bai
│   ├── loup # loup files
│   │   ├── {sample_A}_dloupe.dloupe 
│   │   ├── {sample_B}_dloupe.dloupe 
│   │   └── {sample_C}_dloupe.dloupe 
│   └── summary # summary files
│       ├── {sample_A}_web_summary.html 
│       ├── {sample_B}_web_summary.html 
│       └── {sample_C}_web_summary.html 
├── MEDALT # CNV evolutionary analysis results from the basic run mode
│   ├── {sample_A}
│   │   ├── gene.LSA.txt # list of genes associated with CNA 
│   │   ├── LSA.tree.pdf # lineage tracing of CNA
│   │   ├── medalt.cys # cytoscpe accessible file
│   │   ├── medalt.pdf # lineage tracing of single cell 
│   │   ├── segmental.LSA.txt # list of CNA
│   │   └── other output files
│   ├── {sample_B}
│   │   ├── medalt.cys
│   │   ├── medalt.pdf
│   │   └── other output files
│   ├── {sample_C}
│   │   ├── medalt.cys
│   │   ├── medalt.pdf
│   │   └── other output files
├── MEDALT_group # CNV evolutionary analysis results from the cluster based run mode
│   ├── {sample_A}
│   │   ├── gene.LSA.txt # list of genes associated with CNA 
│   │   ├── LSA.tree.pdf # lineage tracing of CNA
│   │   ├── medalt.group.cys # cytoscpe accessible file
│   │   ├── medalt.group.force.directed.cys # cytoscpe accessible file with force directed layout
│   │   ├── medalt.group.force.directed.pdf # lineage tracing of cluster with force directed layout
│   │   ├── medalt.group.pdf # lineage tracing of cluster
│   │   ├── segmental.LSA.txt # list of CNA
│   │   └── other output files
│   ├── {sample_B}
│   │   ├── medalt.group.cys
│   │   ├── medalt.group.force.directed.cys
│   │   ├── medalt.group.force.directed.pdf
│   │   ├── medalt.group.pdf
│   │   └── other output files
│   ├── {sample_C}
│   │   ├── medalt.group.cys
│   │   ├── medalt.group.force.directed.cys
│   │   ├── medalt.group.force.directed.pdf
│   │   ├── medalt.group.pdf
│   │   └── other output files
├── MEDALT_patient # CNV evolutionary analysis results from the cluster based and patient merged run mode
│   ├── {patient_A}
│   │   ├── gene.LSA.txt # list of genes associated with CNA 
│   │   ├── group_filter.txt
│   │   ├── LSA.tree.pdf # lineage tracing of CNA
│   │   ├── medalt.patient.cys # cytoscpe accessible file
│   │   ├── medalt.patient.force.directed.cys # cytoscpe accessible file with force directed layout
│   │   ├── medalt.patient.force.directed.pdf # lineage tracing of patient merged clusters with force directed layout
│   │   ├── medalt.patient.pdf # lineage tracing of patient merged clusters
│   │   ├── segmental.LSA.txt # list of CNA
│   │   └── other output files
├── {sample_A}
│   ├── cellranger-dan output files
├── {sample_B}
│   ├── cellranger-dan output files
├── {sample_C}
│   ├── cellranger-dan output files
└── reanalysis # cellranger-dna renalysis after noise filtering
    ├── link # main output files of cellranger-dna reanalysis of all samples
    │   ├── alarm # alarm files
    │   │   ├── {sample_A}_alarms_summary.txt 
    │   │   ├── {sample_B}_alarms_summary.txt 
    │   │   └── {sample_C}_alarms_summary.txt 
    │   ├── loup # loup files
    │   │   ├── {sample_A}_dloupe.dloupe 
    │   │   ├── {sample_B}_dloupe.dloupe 
    │   │   └── {sample_C}_dloupe.dloupe 
    ├── {sample_A}
    │   ├── cellranger-dan renalysis output files
    ├── {sample_B}
    │   ├── cellranger-dan renalysis output files
    ├── {sample_C}
    │   ├── cellranger-dan renalysis output files  
    └── tsne # tsne plots
        ├── {sample_A}_plotly_tsne_group.html # Rmarkdown file of tsne ploting with cluster information
        ├── {sample_A}_plotly_tsne.html # Rmarkdown file of tsne ploting
        ├── {sample_A}_tsne_group.pdf # tsne plot with cluster information
        ├── {sample_A}_tsne.pdf # tsne plot
        ├── {sample_B}_plotly_tsne_group.html 
        ├── {sample_B}_plotly_tsne.html 
        ├── {sample_B}_tsne_group.pdf 
        ├── {sample_B}_tsne.pdf
        ├── {sample_C}_plotly_tsne_group.html 
        ├── {sample_C}_plotly_tsne.html 
        ├── {sample_C}_tsne_group.pdf 
        └── {sample_C}_tsne.pdf 

```
