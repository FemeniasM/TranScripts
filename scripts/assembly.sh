#!/bin/bash

red=`tput setaf 1`
green=`tput setaf 2`
yellow=`tput setaf 3`
reset=`tput sgr0`

#################
## This script runs rnaSPAdes, TransABYSS and TransLiG to assembly clean reads
##################

if [ "X" = "X$threads" ]; then threads=18; fi
if [ "X" = "X$memmax" ]; then memmax=90000; fi
if [ "X" = "X$lib_format" ]; then lib_format="pe"; fi
if [ "X" = "X$outdir" ]; then outdir=.; fi

abort()
{
    echo >&2 "
${yellow}          -- !!!!!!! -- ${reset}
${red}    -- ASSEMBLY ABORTED  -- ${reset}
${yellow}          -- !!!!!!! -- ${reset}
"
    echo "An error occurred while assembling the sequencess..." >&2
    exit 1
}

trap 'abort' 0
set -e


if [ -s $kmerset ]; then
echo -e "[$(printf '%(%F %T)T\n')][INFO] defining k-mer set"

  if [ -z "$readsize" -o -z "$KMIN" -o -z "$STEPS"  ]; then
    echo -e "${red}[$(printf '%(%F %T)T\n')][ERROR] missing required argument(s) to define k-mer set${reset}"
    print_usage && exit 1
  fi

kmerset=$(Rscript -e '
readsize=as.numeric(commandArgs(TRUE)[1]) #125
KMIN=as.numeric(commandArgs(TRUE)[2])#27
STEPS=as.numeric(commandArgs(TRUE)[3])#5
KMAX=(round(readsize/10) -1 )* 10 + 5

halfr1=round(readsize/2) - 6
halfr2= halfr1 + 12

KSET <- integer() 

for(k in seq(KMIN, KMAX, round((KMAX-KMIN)/STEPS, 0))) {  
  if(k>halfr1 & k<halfr2){
    kh=k+4
    KSET=c(KSET,kh)
  }else{
    KSET=c(KSET,k)
  }
}
cat(KSET,sep = ",")' $readsize $KMIN $STEPS)

fi

KMAX=$(echo "$kmerset" | awk -F',' '{for (i=1;i<=NF;i++) if ($i>max) max=$i} END{print max}')

kmer_translig=$(echo "$kmerset" | awk -F',' 'BEGIN { ORS="," }; { for (i=1;i<=NF;i++) if ($i<32) print $i}' |  sed '$s/,$//')

mkdir -p $outdir/temp


if [[ $lib_format == 'se' ]];then
echo -e "[$(printf '%(%F %T)T\n')][INFO] collecting reads SE"

single_reads=$(ls $clean_reads/* )
echo -e "[$(printf '%(%F %T)T\n')][INFO] escribiendo archivo YAML para rnaSPAdes"

echo "
    [
      {
        orientation: \"fr\",
        type: \"single\",
        single reads: [
          $single_reads 
        ]
      }
    ]
" > $outdir/temp/YAML.yaml


echo -e "[$(printf '%(%F %T)T\n')][INFO] defining arguments for trans-ABYSS"

s_reads=$( find $clean_reads/* | sort | awk 'BEGIN { ORS=" " }; {print $0}')
opt_reads_transabyss="--se $s_reads"

echo -e "[$(printf '%(%F %T)T\n')][INFO] defining arguments for TransLig"

$translig_dir/plugins/fastool/fastool --to-fasta $single_reads > $outdir/temp/single.fa

opt_reads_translig="-s fa -p single -u $outdir/temp/single.fa"

 else
     if [[ $lib_format == 'pe' ]]; then
  echo -e "[$(printf '%(%F %T)T\n')][INFO] collecting reads PE"

  right_reads=$(find "$clean_reads" -name "*R1*" | sort | awk '{print "\"" $0 "\","}' |  sed '$s/,$//')
  left_reads=$(find "$clean_reads" -name "*R2*" | sort | awk '{print "\"" $0 "\","}' |  sed '$s/,$//')

  echo -e "[$(printf '%(%F %T)T\n')][INFO] writing YAML file for rnaSPAdes"

echo "[
      {
        orientation: \"fr\",
        type: \"paired-end\",
        right reads: [
$right_reads
        ],
        left reads: [
$left_reads
        ]
      }
    ]
  " > $outdir/temp/YAML.yaml

  echo -e "[$(printf '%(%F %T)T\n')][INFO] defining arguments for trans-ABYSS"

  reads_transabyss=$(find "$clean_reads" -name "*R1*" | awk -F'_R1' 'BEGIN { ORS=" " }; {print $1"_R1"$2" "$1"_R2"$2}')

  opt_reads_transabyss="--pe $reads_transabyss"

  echo -e "[$(printf '%(%F %T)T\n')][INFO] defining arguments for TransLig"

  r_reads_translig=$(find "$clean_reads" -name "*R1*" | sort | awk 'BEGIN { ORS=" " }; {print $0}' )
  l_reads_translig=$(find "$clean_reads" -name "*R2*" | sort | awk 'BEGIN { ORS=" " }; {print $0}' )

$translig_dir/plugins/fastool/fastool --to-fasta $l_reads_translig > $outdir/temp/reads.left.fa
$translig_dir/plugins/fastool/fastool --to-fasta $r_reads_translig > $outdir/temp/reads.right.fa

opt_reads_translig="-s fa -p pair -l $outdir/temp/reads.left.fa -r $outdir/temp/reads.right.fa"
 
     fi
 fi


echo -e "${yellow}
Assembly with rnaSPAdes
=======================${reset}"

echo -e "[$(printf '%(%F %T)T\n')][INFO] Starting rnaSPAdes assembly with k-mers: $kmerset"

rnaspades_opt="-k $kmerset --dataset $outdir/temp/YAML.yaml -t $threads -m $memmax -o $outdir/temp/rnaspades"

if [ ! -z "$SS" ]; then rnaspades_opt="--ss-fr $rnaspades_opt"; fi

$rnaspades_dir/rnaspades.py $rnaspades_opt

echo -e "[$(printf '%(%F %T)T\n')][INFO] Finishing rnaSPAdes assembly"

echo -e "[$(printf '%(%F %T)T\n')][INFO] Moving files to output directory"

mkdir -p $outdir/rnaspades
mv $outdir/temp/rnaspades/*.fasta $outdir/rnaspades
rm -r $outdir/temp/rnaspades

echo -e "${yellow}
Assembly with trans-ABYSS
=========================${reset}"

echo -e "[$(printf '%(%F %T)T\n')][INFO] Running step 1 of trans-ABYSS: assembly by k-mer"


for k in $(echo $kmerset | sed "s/,/ /g"); do

transabyss_opt="--threads $threads -k $k --outdir $outdir/temp/transabyss_assembly --name transabyss_${k} $opt_reads_transabyss"

if [ ! -z "$SS" ]; then transabyss_opt="--SS $transabyss_opt"; fi

echo -e "[$(printf '%(%F %T)T\n')][INFO] Running trans-ABYSS with parameters: $transabyss_opt"
$transabyss_dir/transabyss $transabyss_opt

done

echo -e "[$(printf '%(%F %T)T\n')][INFO] Running step 2 of trans-ABYSS: merge assemblies"

fasta_transabyss=$( find "$outdir/temp/transabyss_assembly" -name "*-final.fa" | sort | awk 'BEGIN { ORS=" " }; {print $0}' )
echo "$fasta_transabyss"

trab_outdir="$outdir/temp/transabyss_assembly/transabyss_"
prefix_transabyss=$( echo $fasta_transabyss | awk -v p1="$outdir/temp/transabyss_assembly/transabyss_" '{ gsub(p1,"k",$0); print $0 }' | sed 's/-final.fa/./g')

transabyss_merge_opt="$fasta_transabyss --out $outdir/transabyss_assembly/transabyss_assembly_merged.fa --threads $threads --mink $KMIN --maxk $KMAX --prefix $prefix_transabyss --force"

if [ ! -z "$SS" ]; then transabyss_merge_opt="--SS $transabyss_merge_opt"; fi

mkdir -p $outdir/transabyss_assembly

$transabyss_dir/transabyss-merge $transabyss_merge_opt

echo -e "[$(printf '%(%F %T)T\n')][INFO] Finishing the trans-ABYSS pipeline"

echo -e "[$(printf '%(%F %T)T\n')][INFO] Moving files to output directory"


rm -r $outdir/temp/transabyss_assembly


echo -e "${yellow}
Assembly with TransLiG
======================${reset}"

echo -e "[$(printf '%(%F %T)T\n')][INFO] Starting TransLig assembly"

if [ ! -z "$SS" ]; then opt_reads_translig="-m FR $opt_reads_translig"; fi

$translig_dir/TransLiG $opt_reads_translig -k 31 -o $outdir/temp/translig_assembly

#-k $kmer_translig

echo -e "[$(printf '%(%F %T)T\n')][INFO] Finishing TransLig assembly"

echo -e "[$(printf '%(%F %T)T\n')][INFO] Moving files to output directory"

mkdir -p $outdir/translig
mv $outdir/temp/translig_assembly/TransLiG.fa $outdir/translig

echo -e "[$(printf '%(%F %T)T\n')][INFO] Removing temporary files"

rm -r $outdir/temp/translig_assembly
#rm $outdir/temp/*.left.fa

rm -r $outdir/temp

echo "
Assemblies executed with:
TransABYSS v$(python3 $transabyss_dir/transabyss --version)
rnaSPADES $($rnaspades_dir/rnaspades.py --version | awk '{print $4}')
TransLiG $($translig_dir/TransLiG --version | awk 'NR>1 {print $8}')
k-mer set:
$kmerset
" > $outdir/VERSION_ASSEMBLY_PROGRAMS.txt

trap : 0
echo >&2 "${green}
===============================================
>     Assemblies completed successfully       <
===============================================
${reset}"
cat $outdir/VERSION_ASSEMBLY_PROGRAMS.txt

