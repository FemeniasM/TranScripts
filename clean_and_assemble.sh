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
# 	    This script cleans reads and executes assemblers for RNA-seq data		
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
bash clean_and_assemble ${green}[flags]${reset} [<options>]
 Flags:
  ${green}-t${reset} threads [N] (12 default)
  
  ${green}-m${reset} max memory
  
  ${green}-d${reset} path to directory containing subdirectories (folders) for each set of libraries 
                    (e.g. maindir/sp1/*.fastq; maindir/sp2/*.fastq => -d maindir)
  
  ${green}-f${reset} alternatively users can set folder with libraries whitout subdirectories 
                    (e.g., sp1/*.fastq => -f sp1)
  
  ${green}-l${reset} library format ['pe']['se']
  
  ${green}-s${reset} add argument to strand-specific libraries ['ss']
  
  ${green}-r${reset} run mode (optional): ['a'] to 'assembly' pipeline only, or 
                          ['c'] to 'clean' pipeline only. 
                     Default (without -r) runs both 'clean' and 'assembly' pipelines
  
  ${green}-o${reset} path to output directory
  
  ${green}-h${reset} help

${yellow}NOTE: Users must edit the config.file before running the program${reset}
"
  exit 1
}



while getopts ":t:d:s:r:m:l:o:f:h:" opt; do
  case $opt in
    t)
    threads="$OPTARG"
    echo "-t <threads> = $threads"
    ;;
    d)
    dir_libraries=`realpath $OPTARG`
    echo "-d <directory with libraries> = $dir_libraries"
    ;;
    f)
    lib_folder=`realpath $OPTARG`
    echo "-f <folder with libraries> = $lib_folder"
    ;;
    m)
    memmax="$OPTARG"
    echo "-m <max memory> = $memmax"
    ;;
    l)
    lib_format="$OPTARG"
    echo "-l <libraries format 'pe' | 'se'> = $lib_format"
    ;;
    s)
    SS="TRUE"
    echo "-s <strand-specific>"
    ;;
    r)
    run="$OPTARG"
    echo "-r <run mode> = $run"
    ;;
    o)
    outdir=`realpath $OPTARG`
    echo "-o <Output files Path> = $outdir"
    ;;
    h)
    print_usage && exit 1
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

#dir_libraries=~/2022/parthenos/assembly 
#rRNAdb=~/2022/parthenos/index/rRNAdb 


if [ -z "$trimmomatic_dir" -o -z "$rcorrector_dir" -o -z "$TrAssTools" -o -z "$rnaspades_dir" -o -z "$transabyss_dir" -o -z "$translig_dir"  ]; then
echo -e "${red}[$(printf '%(%F %T)T\n')][ERROR] missing required program path(s) in config file ${reset}"
print_usage && exit 1
fi

if [ -z "$kmerset" ]; then 

  if [ -z "$readsize" -o -z "$KMIN" -o -z "$STEPS"  ]; then
  echo -e "${red}[$(printf '%(%F %T)T\n')][ERROR] if kmerset is not defined, then you must define: readsize (reads lenght), KMIN (k-mer min) and STEPS (num of k-mer steps) ${reset}"
  print_usage && exit 1
  fi
  
  kmeropt="readsize=$readsize KMIN=$KMIN STEPS=$STEPS"

else
  
  kmeropt="kmerset=$kmerset"
  
fi  

if [ ! -s $CONFDIR/scripts/assembly.sh ]; then

  echo -e "${red}[$(printf '%(%F %T)T\n')][ERROR] could not find file assembly.sh $fn ${reset}"
  exit 1

fi

if [ ! -s $CONFDIR/scripts/cleanseqs.sh ]; then

  echo -e "${red}[$(printf '%(%F %T)T\n')][ERROR] could not find file cleanseqs.sh $fn ${reset}"
  exit 1

fi


if [ -s $rRNAdb ]; then

echo -e "[$(printf '%(%F %T)T\n')][INFO] checking if the rRNA index is in .fasta format"
perl -ne \
    '$id = />.+/;
    die "Empty $.\n" if $id && $p || $id && eof;
    $p = $id;
    die "[ERROR] Index error: Must be in fasta format... Invalid char $1 ($.)\n" if !$id && /([^A-Z\n])/
    ' $rRNAdb
    
echo -e "[$(printf '%(%F %T)T\n')][INFO] indexing rRNA database for Bowtie2"

mkdir $outdir/rRNAdb_index
bowtie2-build --threads $threads -f $rRNAdb $outdir/rRNAdb_index

fi

if [ -d $dir_libraries ]; then

  nam_folder=$(find $dir_libraries/* -type d -printf "%f\n")
  
  if [ -z "$run"  ]; then
  
    for n in $nam_folder; do
    echo -e "[$(printf '%(%F %T)T\n')][INFO] run step 1 (sequences cleaning) to $n sample"
    
    env lib_folder=$dir_libraries/${n} rRNAdb=$rRNAdb trimmomatic_dir=$trimmomatic_dir adapters=$adapters rcorrector_dir=$rcorrector_dir \
    TrAssTools=$TrAssTools threads=$threads lib_format=$lib_format outdir=$outdir/step1_${n} $CONFDIR/scripts/cleanseqs.sh
    
    echo -e "[$(printf '%(%F %T)T\n')][INFO] run step 2 (assembly) to $n sample"
    
    env clean_reads=$outdir/step1_${n}/clean_reads rnaspades_dir=$rnaspades_dir transabyss_dir=$transabyss_dir translig_dir=$translig_dir \
    threads=$threads lib_format=$lib_format outdir=$outdir/step2_${n} $kmeropt $CONFDIR/scripts/assembly.sh
    
    done
  fi
  
  
  if [[ $run == 'c'  ]]; then
  
    for n in $nam_folder; do
    echo -e "[$(printf '%(%F %T)T\n')][INFO] run sequences cleaning to $n sample"
    
    env lib_folder=$dir_libraries/${n} rRNAdb=$rRNAdb trimmomatic_dir=$trimmomatic_dir adapters=$adapters rcorrector_dir=$rcorrector_dir \
    TrAssTools=$TrAssTools threads=$threads lib_format=$lib_format outdir=$outdir/${n} $CONFDIR/scripts/cleanseqs.sh
    
    done
  fi
  
  if [[ $run == 'a'  ]]; then
  
    for n in $nam_folder; do
    echo -e "[$(printf '%(%F %T)T\n')][INFO] run assembly to $n sample. 
    ${yellow}[$(printf '%(%F %T)T\n')][WARNING] input sequences without cleaning process${reset}"
    
    env clean_reads=$dir_libraries/${n} rnaspades_dir=$rnaspades_dir transabyss_dir=$transabyss_dir translig_dir=$translig_dir \
    threads=$threads lib_format=$lib_format outdir=$outdir/${n} $kmeropt $CONFDIR/scripts/assembly.sh
  
    done
  fi

else

  if [ -d $lib_folder ]; then
  
    if [ -z "$run"  ]; then
  
    echo -e "[$(printf '%(%F %T)T\n')][INFO] run step 1 (sequences cleaning) to $n sample"
    
    env lib_folder=$lib_folder rRNAdb=$rRNAdb trimmomatic_dir=$trimmomatic_dir adapters=$adapters rcorrector_dir=$rcorrector_dir \
    TrAssTools=$TrAssTools threads=$threads lib_format=$lib_format outdir=$outdir $CONFDIR/scripts/cleanseqs.sh
    
    echo -e "[$(printf '%(%F %T)T\n')][INFO] run step 2 (assembly) to $n sample"
    
    env clean_reads=$outdir/clean_reads rnaspades_dir=$rnaspades_dir transabyss_dir=$transabyss_dir translig_dir=$translig_dir \
    threads=$threads lib_format=$lib_format outdir=$outdir $kmeropt $CONFDIR/scripts/assembly.sh
    
    fi
  
    if [[ $run == 'c'  ]]; then
  
    echo -e "[$(printf '%(%F %T)T\n')][INFO] run sequences cleaning to $n sample"
    
    env lib_folder=$lib_folder rRNAdb=$rRNAdb trimmomatic_dir=$trimmomatic_dir adapters=$adapters rcorrector_dir=$rcorrector_dir \
    TrAssTools=$TrAssTools threads=$threads lib_format=$lib_format outdir=$outdir $CONFDIR/scripts/cleanseqs.sh
    
    fi
  
    if [[ $run == 'a'  ]]; then
  
    echo -e "[$(printf '%(%F %T)T\n')][INFO] run assembly to $n sample. 
    ${yellow}[$(printf '%(%F %T)T\n')][WARNING] input sequences without cleaning process${reset}"
    
    env clean_reads=$lib_folder rnaspades_dir=$rnaspades_dir transabyss_dir=$transabyss_dir translig_dir=$translig_dir \
    threads=$threads lib_format=$lib_format outdir=$outdir $kmeropt $CONFDIR/scripts/assembly.sh
  
    else
  
    echo -e "${red}[$(printf '%(%F %T)T\n')][ERROR] An input directory (-d) or folder (-f) must be specified: $fn ${reset}"
    exit 1
  
    fi

  fi

fi

trap : 0
echo >&2 "${green}
=================================================
>       PIPELINE COMPLETED SUCCESSFULLY         <
<              ------------------               >
=================================================
"
