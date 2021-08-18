# RNA_seq_analysis_pipe  
This is for the QC matrix construction, data analysis and visualization for RNA-seq data.  
This is forked RNA-seq pipeline from NIEHS TaRGETII Environmental Epigenetics project.
Current version: `v4`   

Contributor: Cheng Lyu, Shaopeng Liu, Bo Zhang, and Bongsoo Park

For any question, please contact genomicspark@gmail.com  

## Overview
The TaRGET II RNA-seq pipeline is designed for standardized data processing of mouse RNA-seq samples. It incorporates automatic quality control, generates user friendly files for computational analysis and outputs genome browser tracks for data visualization. To ensure consistent and reproducible data processing, the entire workflow, associated software and libraries are built into a singularity image, which can be run on computational clusters with job submission as well as on stand-alone machines. Pipeline installation requires minimal user input. All the software and genome references used for TaRGET II RNA-seq data processing are included in the pipeline image. The pipeline supports both single- and paired-end data, it accepts FASTQ files, performs alignments, gene features summary and data visualization.

## Software Used in the Pipeline
cutadapt (v1.16) was used to find and remove adapter sequence, and other types of unwanted sequence from high-throughput sequencing reads.

fastqc (v0.11.7) was used to provide quality control check s on raw sequence data.

star (v2.5.4b) was used to map sequence reads against the mouse reference genome.

RSeQC was used to comprehensively evaluate high throughput sequence data.

featureCounts (v1.5.1) was used to count reads to genomic features such as genes.

RSEM (v1.3.0) was used to estimate gene and isoform expression levels from RNA-Seq data.

samtools (v1.9) and bedtools (v2.29) were used to manipulate the sequence alignments and genome features.

## Genomic References Used in the Pipeline
The following mouse genome and gene references were built into the singularity image and used for TaRGET RNA-seq data processing:

STAR indexes of mouse genome (mm10).
Gtf file of vM22 gene annotation from GENCODE.
RSEM references using GENCODE vM22 annotations.
Chromosome size file of mm10.

## Usage: 
Singularity 2-step solution (easiest way)  

Step1. download singularity container (you only need download the containcer for once, then you can use them directly):  
####  
```bash
# download image from local server:  
wget http://regmedsrv1.wustl.edu/Public_SPACE/shaopengliu/Singularity_image/rna-seq/rna-seq_mm10_v4.simg  
```

or for human

```bash
# download image from local server:
wget http://regmedsrv1.wustl.edu/Public_SPACE/shaopengliu/Singularity_image/rna-seq/hg38_rna-seq.simg
```

Step2. process data by the singularity image: 
#### Please run at same directory with your data OR the soft link of your data    
```bash
singularity run -H ./:/scratch rna-seq_mm10_v4.simg -r <SE/PE> -g <mm10>  -o <read_file1>  -p <read_file2>    
```

That's it!

#parameters:  
`-h`: help information  
`-r`: SE for single-end, PE for paired-end  
`-g`: genome reference, one simg is designed for ONLY one species due to the file size. For now the supported genoms are: <mm10/mm9/hg19/hg38/danRer10> (only mm10 in the example).  
`-o`: reads file 1 or the SE reads file, must be ended by .fastq or .fastq.gz or .sra (for both SE and PE)  
`-p`: reads file 2 if input PE data, must be ended by .fastq or .fastq.gz  
`-a`: ADAPT1 for cutadapt  
`-b`: ADAPT2 for cutadapt, if not specified there will be ONLY quality trimming rather than adapter trimming    
`-t`: (optional) specify number of threads to use, default 24  

e.g:
a) mm10 SE data A.fastq  
```bash
singularity run -H ./:/scratch <path_to_simg>  -r SE -g mm10 -o A.fastq  
```
b) hg38 PE data B_1.fastq B_2.fastq  
```bash
singularity run -H ./:/scratch <path_to_simg>  -r PE -g hg38 -o B_1.fastq  -p B_2.fastq  
```
c) danRer10 PE data in sra file C.sra  
```bash
singularity run -H ./:/scratch <path_to_simg> -r PE -g danRer10 -o C.sra  
```


