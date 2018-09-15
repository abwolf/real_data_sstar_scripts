#!/bin/bash
#SBATCH --get-user-env
source ~/.bashrc

paramsfile=~/Sstar_files/bin/params_YRI_1_abw.sh
source $paramsfile

mkdir -p $outputdir
mkdir -p $outputdir/RegionFiles
mkdir -p $outputdir/SstarSigFiles

# refInds=$( for i in $(seq 1 $num_YRI); do cat $ind_info | awk 'BEGIN {OFS="\t"} NR!=1 {if($2=="YRI") print$0}' | shuf -n 1 - | cut -f 1 ; done )
# echo $refInds >> $outputdir/Ref.sample_list

refInds=$( cat ./Ref.sample_list | awk 'BEGIN {OFS="\t"} NR=='$1' {print $0}')
echo $refInds

sbatch --array=1-250%25 --export=refInds="$refInds",paramsfile="$paramsfile" ~/Sstar_files/bin/calc_sstar_neanderthal_with_matchpvals_qsub.sh

echo FIN


# # 2. Calculate S* p-values
# # and determine significance at several reasonable thresholds
#
# ls -C RegionFiles/ | cut -f 6 -d '_' | cut -f 1 -d '.' \
#    | sort -n - | uniq - | grep -v start \
#    | awk 'BEGIN {OFS="\t"} {print "1",$1,($1+999999)}' \
#    > Target.regions_list
#
# reflistfile=Ref.sample_list
# regionslistfile=Target.regions_list
#
# for aref in $( cat $reflistfile | cut -f 1 -d ' ' );  do
#     if [ ! -d SstarSigFiles/$aref ]; then
#         for astart in $(cat $regionslistfile | cut -f 2); do
#             echo $aref      $astart
#             R --slave --args \
#             $tag \
#             $neand_callable_regions_file \
#             $recomb_rates_file \
#             $outputdir \
#             glm \
#             $sstar_null_model \
#             $astart \
#             $chr \
#             $aref \
#             $n_req_window_snps \
#         < $scriptdir/compute_pvalues_from_indfiles.r
#         done
#     elif [ -d SstarSigFiles/$aref ]; then
#         echo "$aref directory already exists"
#         continue
#     fi
# done
