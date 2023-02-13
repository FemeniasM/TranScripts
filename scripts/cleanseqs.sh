#!/bin/bash

red=`tput setaf 1`
green=`tput setaf 2`
yellow=`tput setaf 3`
reset=`tput sgr0`


if [ "X" = "X$threads" ]; then threads=12; fi
if [ "X" = "X$lib_format" ]; then lib_format="pe"; fi
if [ "X" = "X$outdir" ]; then outdir=.; fi

adapters=$trimmomatic_dir/adapters/TruSeq3-PE.fa

abort()
{
    echo >&2 "
${yellow}        -- !!!!!!! -- ${reset}
${red} -- CLEANING SEQUENCES ABORTED  -- ${reset}
${yellow}        -- !!!!!!! -- ${reset}
"
    echo "An error occurred while cleaning the sequences..." >&2
    exit 1
}

trap 'abort' 0

set -e

file_extension_error ()
{
echo -e "${red}[$(printf '%(%F %T)T\n')][ERROR] FileExtensionError: Invalid extension${reset}"
echo -e "${green}Supported extensions are: <.fq> or <.fastq> or <.fq.gz> or <.fastq.gz>${reset}"
return     
}

file_name_error ()
{
echo -e "${red}[$(printf '%(%F %T)T\n')][ERROR] Filename Error: Paired-end file names should contain _R1 _R2${reset}"
echo -e "${green}Example: sample_R1.fq.gz, sample_R2.fq.gz${reset}"
return 
}


if [ -z "$lib_folder" -o -z "$rRNAdb" -o -z "$trimmomatic_dir" -o -z "$rcorrector_dir" -o -z "$TrAssTools" -o -z "$lib_format"  ]; then
    echo -e "${red}[$(printf '%(%F %T)T\n')][ERROR] missing required argument(s)${reset}"
    abort && exit 1
fi


echo "[$(printf '%(%F %T)T\n')][INFO] checking fastq files extensions.."
if [[ $lib_format == 'se' ]];then
    for fn in $lib_folder/*; do
    if [[ (${fn#*.} == "fastq.gz") || (${fn#*.} == "fq.gz") || (${fn#*.} == "fastq") || (${fn#*.} == "fq") ]]
        then
        echo "[$(printf '%(%F %T)T\n')][INFO] libraries extension format ok..."
        else 
        file_extension_error
        echo -e "${red}[$(printf '%(%F %T)T\n')][ERROR] check extension format: $fn ${reset}"
        exit 1
    fi
    done
else 
 if [[ $lib_format == 'pe' ]]; then
    for fn in $lib_folder/*R1*; do
        if [[ (${fn#*.} == "fastq.gz") || (${fn#*.} == "fq.gz") || (${fn#*.} == "fastq") || (${fn#*.} == "fq") ]]; then
         echo "[$(printf '%(%F %T)T\n')][INFO] $fn: extension format ok..."
            sample_name=$(basename ${fn} | sed 's/.fastq.gz\|.fq.gz\|.fastq\|.fq//g')
            if ls $lib_folder/$sample_name* | grep -q -e "_R1" -e "_R2"; then
            echo "[$(printf '%(%F %T)T\n')][INFO] $fn: name format ok..."
            else
               file_name_error
               echo -e "${red}[$(printf '%(%F %T)T\n')][ERROR] check name format: $fn ${reset}"
               exit 1
            fi
         else
         file_extension_error
         echo -e "${red}[$(printf '%(%F %T)T\n')][ERROR] check extension format: $fn ${reset}"
         exit 1
         fi
     done
 else
    echo -e "${red}[$(printf '%(%F %T)T\n')][ERROR] Wrong library format, it should be 'pe' or 'se'${reset}\n"    
    exit 1
 fi
fi


mkdir -p $outdir/temp


echo -e "[$(printf '%(%F %T)T\n')][INFO] collecting PE reads"
echo -e "${yellow}
Running Trimmomatic
======================${reset}"

for fn in $lib_folder/*R1*; do
  if [[ (${fn#*.} == "fastq.gz") || (${fn#*.} == "fq.gz") || (${fn#*.} == "fastq") || (${fn#*.} == "fq") ]]; then
    sample_name=$(basename ${fn} | sed 's/_R1.fastq.gz\|_R1.fq.gz\|_R1.fastq\|_R1.fq//g')
  if ls $lib_folder/$sample_name* | grep -q -e "_R1" -e "_R2"; then
    ext=$(basename ${fn} | awk -F'_R1' '{print $2}')

R1=${sample_name}_R1${ext}
R2=${sample_name}_R2${ext}

echo -e "[$(printf '%(%F %T)T\n')][INFO] Running Trimmomatic on ${sample_name}"

java -jar $trimmomatic_dir/trimmomatic-0.39.jar PE -threads $threads -phred33 -summary $outdir/temp${sample_name}_trimm.sum $lib_folder/$R1 $lib_folder/$R2 $outdir/temp/paired_${R1} $outdir/temp/unpaired_${R1} $outdir/temp/paired_${R2} $outdir/temp/unpaired_${R2} ILLUMINACLIP:$adapters:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

    else
    file_name_error
    exit 1
  fi
    else
    file_extension_error
    exit 1
  fi
done

rm $outdir/temp/unpaired_*

echo -e "${yellow}
Running Rcorrector
=====================${reset}"

for fn in $outdir/temp/paired_*_R1*; do
    sample_name=$(basename ${fn} | sed 's/_R1.fastq.gz\|_R1.fq.gz\|_R1.fastq\|_R1.fq//g')
    ext=$(basename ${fn} | awk -F'_R1' '{print $2}')

R1=${sample_name}_R1${ext}
R2=${sample_name}_R2${ext}

echo -e "[$(printf '%(%F %T)T\n')][INFO] Running Rcorrector on ${sample_name}"
$rcorrector_dir/run_rcorrector.pl -1 $outdir/temp/$R1 -2 $outdir/temp/$R2 -t $threads -od $outdir/temp/

done


echo -e "${yellow}[$(printf '%(%F %T)T\n')][INFO] Removing unfixable sequences${reset}"

for fn in $outdir/temp/*_R1.cor*; do
    sample_name=$(basename ${fn} | sed 's/_R1.cor.fq.gz//g')

R1=${sample_name}_R1.cor.fq.gz
R2=${sample_name}_R2.cor.fq.gz

echo -e "[$(printf '%(%F %T)T\n')][INFO] cleaning the library: ${sample_name}(_R1/_R2).cor.fq.gz"

python2 $TrAssTools/FilterUncorrectabledPEfastq.py -1 $outdir/temp/$R1 -2 $outdir/temp/$R2

mv $(pwd)/corrected_* $outdir/temp

done

rm $outdir/temp/paired_*

echo -e "[$(printf '%(%F %T)T\n')][INFO] compressing .fasta files"
gzip $outdir/temp/corrected_*

echo -e "${yellow}[$(printf '%(%F %T)T\n')][INFO] removing reads that map to rRNA${reset}"

if [ -s $rRNAdb ]; then

echo -e "[$(printf '%(%F %T)T\n')][INFO] Running Bowtie2"

  for fn in $outdir/temp/corrected_paired_*_R1.cor.fq.gz; do
  sample_name=$( basename ${fn} | sed 's/_R1.cor.fq.gz//g' | sed 's/corrected_un_paired_//g')
  R1=corrected_un_paired_${sample_name}_R1.cor.fq.gz
  R2=corrected_un_paired_${sample_name}_R2.cor.fq.gz

  echo -e "[$(printf '%(%F %T)T\n')][INFO] Running Bowtie2 on ${sample_name}"
  mkdir -p $outdir/clean_reads
  bowtie2 -p $threads --very-sensitive-local --fr -x $outdir/index/$(basename $rRNAdb | awk -F'.' '{print $1}' ) -1 $outdir/temp/$R1 -2 $outdir/temp/$R2 --un-conc-gz $outdir/clean_reads/${sample_name}_R%.fastq.gz -S $outdir/temp/${sample_name}.sam
  rm $outdir/temp/${sample_name}.sam
  done

else
  echo -e "[$(printf '%(%F %T)T\n')][INFO] Running Bowtie2 with index: $rRNAdb"

  for fn in $outdir/temp/corrected_un_paired_*_R1.cor.fq.gz; do
  sample_name=$( basename ${fn} | sed 's/_R1.cor.fq.gz//g' | sed 's/corrected_un_paired_//g')
  R1=corrected_un_paired_${sample_name}_R1.cor.fq.gz
  R2=corrected_un_paired_${sample_name}_R2.cor.fq.gz

  echo -e "[$(printf '%(%F %T)T\n')][INFO] Running Bowtie2 on ${sample_name}"
  mkdir -p $outdir/clean_reads
  bowtie2 -p $threads --very-sensitive-local --fr -x $rRNAdb -1 $outdir/temp/$R1 -2 $outdir/temp/$R2 --un-conc-gz $outdir/clean_reads/${sample_name}_R%.fastq.gz -S $outdir/temp/${sample_name}.sam
  rm $outdir/temp/${sample_name}.sam
  done

fi

rm -r $outdir/temp
#--quiet --very-sensitive-local --phred33 #bowtie parameters

trap : 0
echo >&2 "${green}
=====================================================
>     The reads have been successfully cleaned:     <
<   The files in the 'clean_reads' folder are ready >
>             to run the assemblers                 <
=====================================================
"

