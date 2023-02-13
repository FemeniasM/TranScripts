#!/bin/bash
# Using getopt

abort() {
  echo -e >&2 '
===============
/// ABORTED ///
===============
'
  echo -e "${red}[$(printf '%(%F %T)T\n')][ERROR]:An error occurred, check the input files. Exiting...${reset}" >&2
  exit 1
}

trap 'abort' 0
set -e

#########################################################################################
# 	    This script executes quality control fro assemblers for RNA-seq data		
#   For questions or suggestions please contact to Martin Femenias (mmfemenias@gmail.com)
#########################################################################################


CONFDIR=$(dirname `readlink -f $0`) && export CONFDIR=${CONFDIR}

red=`tput setaf 1`
green=`tput setaf 2`
yellow=`tput setaf 3`
reset=`tput sgr0`


print_usage () {
  echo "$0
Usage:  
tranScripts ${green}[assembly_qc]${reset} ${yellow}[flags]${reset}

Flags:
   ${yellow}-f${reset} fasta files comma separated (without spaces)
   ${yellow}-r${reset} library format ['pe']['se'] 
   ${yellow}-l${reset} libaries folder with fastq files
   ${yellow}-p${reset} numeric vector from 0 to 1 indicating proportion of readings to map (default 1:100%)
   ${yellow}-t${reset} threads [N] (1 default)  
   ${yellow}-o${reset} Output directory used in quant
   ${yellow}-help${reset}

${yellow}NOTE: Users must edit the config.file before running the program
For this analysis, all .fastq files must be in the same folder 
${reset}
"
  exit 1
}


print_help () {
    echo "$0
Usage:  
tranScripts ${green}[assembly_qc]${reset} ${yellow}[flags]${reset}

Flags:
   ${yellow}-f${reset} fasta file
   ${yellow}-r${reset} Folder with reads
   ${yellow}-o${reset} Output directory used in quant
   ${yellow}-help${reset}

This function returns summary statistics for a de novo assembled transcriptome
Estatistics:
           Total 'genes'  Total number of genes
       Total transcripts  Total number of transcripts
           Genes* >500bp  Genes with more than 500bp in length
          Genes* >1000bp  Genes with more than 1000bp in length
      Transcripts >500bp  Transcripts with more than 500bp in length
     Transcripts >1000bp  Transcripts with more than 1000bp in length
              Percent CG  Average of percentage CG in each transcript
              Percent AT  Average of percentage AT in each transcript
              N50 genes*  Number of genes with lengths greater than the N50 gene
             N50L genes*  N50 gene length
              N90 genes*  Number of genes with lengths greater than the N90 gene
             N90L genes*  N90 gene length
         N50 transcripts  Number of genes with lengths greater than the N50 transcript
        N50L transcripts  N50 transcript length
         N90 transcripts  Number of genes with lengths greater than the N90 transcript
        N90L transcripts  N90 transcript length
    Median gene lengths*  Median of gene lengths
Median transcript lengths Median of transcript lengths
                 Ex90N50  N50 gene length for subset data containing 90% of expression
    Gene percentage Ex90  Gene percentage containing 90% expression
            Mapping rate  Mean of Salmon libraries quantification mapping rate 

" 
}


abort() 
{
    echo -e >&2 "
${red}
_____________________
  summarise aborted 
=====================
${reset}"
    echo -e "${red}An error occurred. Check the arguments, and paths of the config.file file
    Exiting...${reset}" >&2
    exit 1
}

trap 'abort' 0
set -e

echo "${green} Summarizing counts in TESSA with arguments:${reset}"
while getopts ":f:r:l:p:t:o:h:" opt; do
    case $opt in

        f)
            fasta="$OPTARG"
            echo "-f <fasta file> = $fasta_files"
            ;;
        r)
            lib_format="$OPTARG"
            echo "-r <libraries format 'pe' | 'se'> = $lib_format"
            ;;
        l)
            lib_folder=`realpath $OPTARG`
            echo "-l <fastq libraries path> = $lib_folder"
            ;;

        p)
            val=$OPTARG ;;
            echo "-p <proportion sample of reads> = $val"
            ;;
        t)
            threads="$OPTARG"
            echo "-t <threads> = $threads"
            ;;

        o)
            outfolder=`realpath $OPTARG`
            echo "-o <output directory> = $outfolder"
            ;;
        h)  
            print_help && exit 1
            ;;
        \?)
            echo "Invalid option: -$OPTARG"
            print_usage && exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument."
            print_usage && exit 1
            ;;
            
    esac
done

source $CONFDIR/config.file


if [ ! -s $CONFDIR/scripts/run_stats.sh ]; then

  echo -e "${red}[$(printf '%(%F %T)T\n')][ERROR] could not find file run_stats.sh $fn ${reset}"
  exit 1

fi

if [ ! -s $CONFDIR/scripts/run_stats.r ]; then

  echo -e "${red}[$(printf '%(%F %T)T\n')][ERROR] could not find file run_stats.r $fn ${reset}"
  exit 1

fi


if [ ! -s $CONFDIR/scripts/cleanseqs.sh ]; then

  echo -e "${red}[$(printf '%(%F %T)T\n')][ERROR] could not find file cleanseqs.sh $fn ${reset}"
  exit 1

fi

for sn in $(echo -e "${fasta_files}" | sed 's/,/ /g'); do

tme_nam=$(basename $sn | sed 's/.fasta\|.fa//g')

echo -e "[$(printf '%(%F %T)T\n')][INFO] run analysis to $tme_nam transcriptome"

$CONFDIR/scripts/run_stats.sh -f $sn -r $lib_format -l $lib_folder -p $val -t $threads -o $outfolder/${tme_nam}

done


trap : 0
echo >&2 "${green}
  Assemby summary statistics is in $outfolder            
"

