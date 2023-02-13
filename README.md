<h1 align="center">TranScripts</h1>
<div align="center">

<i>A suite of tools for reads cleanup, de novo assembly, and quality control of transcriptomes</i>
</div>

### About
The _clean_and_assemble_ script executes a pipeline for cleaning Illumina RNAseq data reads (single-end or paired-end) and assembly with rnaSPAdes, TransABYSS and TransLiG.

The _assebly_qc_ script executes a pipeline to summarize statistics for a *de novo* assembled transcriptome

### Requirements:
The following programs must be installed on your system:
- Trimommatic (v0.39)
- RCorrector
- TranscriptomeAssemblyTools
- Bowtie2
- rnaSPAdes
- TransABYSS
- TransLiG
- R (if an automated search for k-mers is done)

In addition, an rRNA database (e.g. obtained from [SILVA](https://www.arb-silva.de/)) is required.
For assebly_qc analysis is required salmon

### How to use:
Clone this directory to your destination folder. Users must edit the config.file and set the paths to k-mer programs and arguments. Once the file has been edited with the appropriate arguments, just run the script as follow:

#### 1) clean_and_assemble script:

```
bash clean_and_assemble [flags] [<options>]
```
The allowed flags are the following:

```
Usage:

bash clean_and_assemble [flags] [<options>]

  -t threads [N] (12 default)
  
  -m max memory
  
  -d path to directory containing subdirectories (folders) for each set of libraries 
                    (e.g. maindir/sp1/*.fastq; maindir/sp2/*.fastq => -d maindir)
  
  -f alternatively users can set folder with libraries whitout subdirectories 
                    (e.g., sp1/*.fastq => -f sp1)
  
  -l library format paired-end ['pe'] or single-end ['se']
  
  -s add argument to strand-specific libraries ['ss']
  
  -r run mode (optional): ['a'] to 'assembly' pipeline only, or 
                          ['c'] to 'clean' pipeline only. 
                          Default (without -r) runs both 'clean' and 'assembly' pipelines
  
  -o path to output directory
  
  -h help

```
The users can set a directory containing subdirectories (folders) for each set of libraries or, alternatively, users can set folder with libraries whitout subdirectories. Example:

-d `dir` option to:
```
dir
 |_subdir1_
 |         |-> *rep1.fastq/.fastq.gz files (PE or SE)
 |         |-> *rep2.fastq/.fastq.gz files (PE or SE)
 |         |-> *rep3.fastq/.fastq.gz files (PE or SE)
 |
 |_subdir2_
 |         |-> *rep1.fastq/.fastq.gz files (PE or SE)
 |         |-> *rep2.fastq/.fastq.gz files (PE or SE)
 |         |-> *rep3.fastq/.fastq.gz files (PE or SE)
 |
 |_subdir3_
           |-> *rep1.fastq/.fastq.gz files (PE or SE)
           |-> *rep2.fastq/.fastq.gz files (PE or SE)
           |-> *rep3.fastq/.fastq.gz files (PE or SE)

```
-f `folder1` option to:

```
folder1
   |_____
         |-> *rep1.fastq/.fastq.gz files (PE or SE)
         |-> *rep2.fastq/.fastq.gz files (PE or SE)
         |-> *rep3.fastq/.fastq.gz files (PE or SE)

```
Paired-end reads must **mandatorily** contain `_R1` and `_R2` indicators at the end of the name (e.g., `sample_R1.fq.gz`, `sample_R2.fq.gz`). Supported extensions are: `.fq`, `.fastq`, `.fq.gz` or `.fastq.gz`.

The default option executes both pipelines (clean & assembly), however the user can execute the pipelines individually with the `-r` argument (e.g., `-r a` to 'assembly' pipeline only, or `-r c` to 'clean' pipeline only).



#### 2) assembly_qc script:

For _assembly_qc_ analysis, all `.fastq` files must be in the same folder and run:

```
Usage:  
${green}assembly_qc.sh${reset} ${yellow}[flags]${reset}

Flags:
  -f fasta files comma separated (without spaces)
  -r library format ['pe']['se'] 
  -l libaries folder with fastq files
  -p numeric vector from 0 to 1 indicating proportion of reads to map (default 1:100%)
  -t threads [N] (1 default)  
  -o Output directory used in quant
  -h  help
```

To run the analysis on different transcriptomes, users only have to declare the paths with the `-f` argument separated by commas without spaces. The _assembly_qc_ script returns the following summarized statistics:

```
Statistics:
              Total 'genes'   Total number of genes
         Total transcripts    Total number of transcripts
           Genes >500bp(*)    Genes with more than 500bp in length
          Genes >1000bp(*)    Genes with more than 1000bp in length
        Transcripts >500bp    Transcripts with more than 500bp in length
       Transcripts >1000bp    Transcripts with more than 1000bp in length
                Percent CG    Average of percentage CG in each transcript
                Percent AT    Average of percentage AT in each transcript
              N50 genes(*)    Number of genes with lengths greater than the N50 gene
             N50L genes(*)    N50 gene length
              N90 genes(*)    Number of genes with lengths greater than the N90 gene
             N90L genes(*)    N90 gene length
           N50 transcripts    Number of genes with lengths greater than the N50 transcript
          N50L transcripts    N50 transcript length
           N90 transcripts    Number of genes with lengths greater than the N90 transcript
          N90L transcripts    N90 transcript length
    Median gene lengths(*)    Median of gene lengths
 Median transcript lengths    Median of transcript lengths
               Ex90N50(**)    N50 gene length for subset data containing 90% of expression
      Gene percentage Ex90    Gene percentage containing 90% expression
              Mapping rate    Mean of Salmon libraries quantification mapping rate 
```
(*) the gene lengths is taken as the longest transcript for each gene

(**) the gene length is taken as the expression-weighted mean of isoform lengths


