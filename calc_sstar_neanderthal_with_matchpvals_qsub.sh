#!/bin/bash
#SBATCH --get-user-env
#SBATCH --mem=5G
#SBATCH --qos=1day
#SBATCH --time=22:00:00
#SBATCH --output=/scratch/tmp/abwolf/Sstar/chr1_variable_ref/calc_sstar_neanderthal_with_matchpvals_qsub.sh.%A_%a.o
source ~/.bashrc

date
echo $SLURM_JOB_NAME
echo $SLURM_ARRAY_TASK_ID
# Import paramsfile value from sbatch command line
paramsfile=$paramsfile
#paramsfile=~/Sstar_files/bin/params_abw.sh
echo $paramsfile
source $paramsfile
echo ''

## SETUP ##
read -a ARGV <<< $(head -n $SLURM_ARRAY_TASK_ID $regionsfile | tail -n 1) # grab line from regionsfile
chr=${ARGV[0]}
s=${ARGV[1]}
e=${ARGV[2]}
regiontag="$tag"_chr_"$chr"_start_"$s"

# Import refInds value from sbatch command line
refInds=$refInds

## BUILD PYTHON COMMAND##
cmd=$(echo " ~/software/anaconda2/bin/python $windowcalcscript \n
    --vcf-has-illumina-chrnums \n
    -vcfz $vcf_modern \n
    -indf $ind_info \n
    -target-pops $targetpops \n
    -ref-inds $refInds \n
    --archaic-vcf $vcf_arch \n
    -p 10 \n
    -s-star \n
    -ancbsg $bsg_ancestral \n
    -winlen 50000 \n
    -winstep 10000 \n
    -ptable $ptable \n
    -table-query "mh" "len" "mapped" \n
    -winchrom $chr \n
    -range $s $e \n
    -x $x1 $x2 $x3 \n
    -r $r1 \n
    -ir $ir1 $ir2 $ir3 \n
    | gzip -c - \n
    > $outputdir/RegionFiles/$regiontag.$(echo $refInds | cut -f 1 -d ' ').windowcalc_out.gz")
#    > ./RegionFiles/$regiontag.$(echo $refInds | cut -f 1 -d ' ').windowcalc_out.gz")

echo -e $cmd

eval $( echo -e $cmd )

date
echo FIN
