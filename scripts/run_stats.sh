#!/bin/bash
red=`tput setaf 1`
green=`tput setaf 2`
yellow=`tput setaf 3`
blue=`tput setaf 4`
magenta=`tput setaf 5`
reset=`tput sgr0`
val=1
print_usage () {
    echo "$0
Usage:  
tranScripts ${green}[assembly_qc]${reset} ${yellow}[flags]${reset}

Flags:
   ${yellow}-f${reset} fasta file
   ${yellow}-r${reset} library format ['pe']['se'] 
   ${yellow}-l${reset} libaries folder with fastq files
   ${yellow}-p${reset} numeric vector from 0 to 1 indicating proportion of readings to map (default 1:100%)
   ${yellow}-t${reset} threads [N] (1 default)  
   ${yellow}-o${reset} Output directory used in quant
   ${yellow}-help${reset}
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
Statistics:
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
            fasta=`realpath $OPTARG`
            echo "-f <fasta file> = $summarise"
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

mkdir -p $outfolder/temp
cd $outfolder/temp

source $CONFDIR/config.file
seqkit=$CONFDIR/scripts/seqkit

echo "${blue}[$(printf '%(%F %T)T\n')] ${magenta}[INFO] ${yellow} importing files and summarizing lengths"


$seqkit fx2tab --length --name -i -i -B AT -B CG -B N  $fasta | 
awk -v OFS="\t" '{j=$0; split($1,a,"t"); $6=a[1]; $7="t"a[2] } 1' > $outfolder/temp/lengths_by_transcripts.txt


echo "${blue}[$(printf '%(%F %T)T\n')] ${magenta}[INFO] ${yellow} running salmon"

mkdir -p $outfolder/temp/salmon_index

  echo -e "${blue}[$(printf '%(%F %T)T\n')] ${magenta}[INFO] ${yellow}Make salmon index"

  $salmon index -p $threads -t $fasta -k 31 -i $outfolder/temp/salmon_index

  echo -e "${blue}[$(printf '%(%F %T)T\n')] ${magenta}[INFO] ${yellow}Index finished successfully"


echo -e "${blue}[$(printf '%(%F %T)T\n')] ${magenta}[INFO] ${yellow} Quantifying with Salmon${reset}"


mkdir -p $outfolder/temp/salmon_quant

if [[ $lib_format == 'se' ]];then
    for fn in $lib_folder/*; do
    if [[ (${fn#*.} == "fastq.gz") || (${fn#*.} == "fq.gz") || (${fn#*.} == "fastq") || (${fn#*.} == "fq") ]]
        then
        
        sample_name=$(basename ${fn} | sed 's/.fastq.gz\|.fq.gz\|.fastq\|.fq//g')
        

        echo "${blue}[$(printf '%(%F %T)T\n')] ${magenta}[INFO] ${yellow} processing sample $sample_name"
        
            if [[ $val -ne 1 ]];then
            
            echo "${blue}[$(printf '%(%F %T)T\n')] ${magenta}[INFO] ${yellow} subsampling $sample_name with proportion $val"
            $seqkit sample -s 1234 -p 0.1 $fn | gzip > ${sample_name}_${val}.fastq.gz
            
            echo "${blue}[$(printf '%(%F %T)T\n')] ${magenta}[INFO] ${yellow} running salmon quant"
            $salmon quant -i $outfolder/temp/salmon_index -l A --gcBias --useVBOpt -r ${sample_name}_${val}.fastq.gz -p $threads --validateMappings 0.65 -o $outfolder/temp/salmon_quant/$sample_name
            rm ${sample_name}_${val}.fastq.gz
            
            else

            echo "${blue}[$(printf '%(%F %T)T\n')] ${magenta}[INFO] ${yellow} running salmon quant"
            $salmon quant -i $outfolder/temp/salmon_index -l A --gcBias --useVBOpt -r ${fn} -p $threads --validateMappings 0.65 -o $outfolder/temp/salmon_quant/$sample_name
    
            fi

        else 
        file_extension_error
        exit 1
    fi
    done

else 
 if [[ $lib_format == 'pe' ]]; then
    for fn in $lib_folder/*R1*; do
        if [[ (${fn#*.} == "fastq.gz") || (${fn#*.} == "fq.gz") || (${fn#*.} == "fastq") || (${fn#*.} == "fq") ]]; then
         sample_name=$(basename ${fn} | sed 's/_R1.fastq.gz\|_R1.fq.gz\|_R1.fastq\|_R1.fq//g')
          if ls $lib_folder/$sample_name* | grep -q -e "_R1" -e "_R2"; then
            
             ext=$(basename ${fn} | awk -F'_R1' '{print $2}')
             R1=${sample_name}_R1${ext}
             R2=${sample_name}_R2${ext}
            
             echo "${blue}[$(printf '%(%F %T)T\n')] ${magenta}[INFO] ${yellow} processing sample $sample_name"

            if [[ $val -ne 1 ]];then
            
            echo "${blue}[$(printf '%(%F %T)T\n')] ${magenta}[INFO] ${yellow} subsampling $sample_name with proportion $val"
            $seqkit sample -s 1234 -p 0.1 $lib_folder/$R1 | gzip > ${sample_name}_${val}_R1${ext}
            $seqkit sample -s 1234 -p 0.1 $lib_folder/$R2 | gzip > ${sample_name}_${val}_R2${ext}
            
            echo "${blue}[$(printf '%(%F %T)T\n')] ${magenta}[INFO] ${yellow} running salmon quant"
            $salmon quant -i $outfolder/temp/salmon_index -l A --gcBias --useVBOpt -1 ${sample_name}_${val}_R1${ext} -2 ${sample_name}_${val}_R2${ext} -p $threads --validateMappings 0.65 -o $outfolder/temp/salmon_quant/$sample_name
            rm ${sample_name}_${val}_R*
            
            else
            
            echo "${blue}[$(printf '%(%F %T)T\n')] ${magenta}[INFO] ${yellow} running salmon quant"
            $salmon quant -i $outfolder/temp/salmon_index -l A --gcBias --useVBOpt -1 $lib_folder/$R1 -2 $lib_folder/$R2 -p $threads --validateMappings 0.65 -o $outfolder/temp/salmon_quant/$sample_name
            
            fi
        
          else

          file_name_error
          exit 1

          fi
        
        else

        file_extension_error
        exit 1

        fi

    done

 else

    echo -e "\n${red} Wrong library format, it should be 'pe' or 'se'${reset}\n"    
    exit 1
 fi

fi

mapping_rate=cat $(find -name 'salmon_quant.log' -type f) | grep "Mapping rate" | awk -F"=" '{print $2}' | awk -F"%" '{print $1}' | awk '{ sum+=$1 }END { print sum/NR }'

echo "${blue}[$(printf '%(%F %T)T\n')] ${magenta}[INFO] ${yellow} calculating statistics for transcriptome $fasta"
Rscript $CONFDIR/scripts/run_stats.r $outfolder/temp $mapping_rate


echo "${blue}[$(printf '%(%F %T)T\n')] ${magenta}[INFO] ${yellow} removing temporary files"
rm -r $outfolder/temp/


trap : 0
echo >&2 "${green}
======================================================
  Assemby summary statistics finished successfully            
======================================================
"
