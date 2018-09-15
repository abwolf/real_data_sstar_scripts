
chr=1
num_YRI=1

# paths
awdir=/Genomics/grid/users/abwolf/
resources=$awdir/Sstar_files/real_data_files.sv/
scriptdir=$awdir/Sstar_files/bin/
windowcalcscript=$awdir/SimulatedDemographic/Sstar/s_star_git/bin/windowed_calculations.py
projdir=/Genomics/akeylab/abwolf/SimulatedDemographic/Sstar/chr1_variable_ref/


tag=YRI_$num_YRI
outputdir=$projdir/Output_$tag
regionsfile=$resources/regions_1Mb_from_hapmap_range.txt
recomb_rates_file=$resources/recombination_rates_per_50Kb_window_allchr.txt


bsg_ancestral=$resources/chimp_chrAll.bsg
#vcf_arch=/Genomics/grid/users/limingli/data/Altai_mq25/chr"$chr"_mq25_mapab100.vcf.gz
vcf_arch=$resources/chr"$chr".altai_neand_mpi_minimal_filtered_lowqual.vcf.gz
vcf_modern=/Genomics/grid/users/limingli/data/1kg-phase3/ALL.chr"$chr".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
ind_info=$resources/integrated_call_samples_v3.20130502.ALL.panel

targetpops="EUR"
refpop="YRI"

samplelistfile=$resources/EUR.sample_list



#refInds=$( cat $ind_info | awk 'BEGIN {OFS="\t"} NR!=1 {if($2=="YRI") print$0}' | shuf -n $num_YRI - | cut -f 1 )

sstar_null_model=$resources/gravel_asn_scale_60k.model.Rdata
neand_callable_regions_file=$resources/windows.50kb.10kb.bed

ptable=$resources/archaic_match_table_files.neand_table.fields_8-10.12-.gz.db
n_req_window_snps=3


## filter files to be provided to -x flag
x1=$resources/cpg2_windows_hg19.bed.zerobased.bbg
x2=$resources/genomicSuperDups.txt.bed.bbg
x3=$resources/hs37m_mask35_50.flt.bed.fixed.bbg

## filter files to be provided with -r flag
r1=$resources/neand_called_bases_x_indels.bbg

## filter files to be provided with -ir flag
## could be more files that should be provided
ir1=$resources/chimp_chrAll.mapped.bbg
ir2=$resources/20141020.pilot_mask.whole_genome.bed.bbg
ir3=$resources/AltaiNea.map35_50.MQ30.Cov.indels.TRF.bed.bbg

# ######################################
# ######################################
# ## ORIGIONAL SCRIPT FROM SELINA
# # files that exist
# samplelistfile=$projdir/ASW_sample_ids.txt
# regionsfile=$projdir/regions_for_matchpval_1kg.txt
#
#
# ## for files and directories that will be created
# tag="ref_LWK"
# outputdir=$projdir/Output_$tag
#
#
# # S-calc options
# # some filenames include chromosome, so when using them need to make
# # sure that $chr is defined
# bsg_ancestral=$resources/Chimp_bsg/chimp_chrAll.bsg
# vcf_arch=$resources/Archaic_genotypes/Neanderthal/chr$chr.altai_neand_mpi_minimal_filtered_lowqual.vcf.gz
# vcf_modern=/net/akey/vol2/wqfu/nobackup/1KGP/ALL.chr$chr.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz
# ind_info=$projdir/integrated_call_samples_v3.20130502.ALL.panel
# targetpops="ASW"
# refpop="LWK"
#
#
# # S* p-value calc options
# sstar_null_model=$resources/Lookup_tables/gravel_asn_scale_60k.model.Rdata
# neand_callable_regions_file=$resources/Lookup_tables/neand.callable.bed
# recomb_rates_file=$svdir/Saudi/Recombrates/sstar6_per_interval_recombination_rates_allchr.txt
#
# # archaic match p-value calc options
# ptable=$resources/Lookup_tables/archaic_match_table_files.neand_table.fields_8-10.12-.gz.db
# n_req_window_snps=6
#
#
# ## filter files to be provided to -x flag
# #/tigress/AKEY/akey_vol1/home/bvernot//tishkoff/filter_files/cpg2_files/cpg2_windows_hg19.bed.zerobased.bbg
# #/tigress/AKEY/akey_vol1/home/bvernot//archaic_exome/data/segdups/2013.01.04/genomicSuperDups.txt.bed.bbg
# #/tigress/AKEY/akey_vol1/home/bvernot//archaic_exome/data/snp_mappability_reich/2013.01.04/hs37m_mask35_50.flt.bed.fixed.bbg
#
#
# ## filter files to be provided with -r flag
# #/tigress/AKEY/akey_vol1/home/bvernot//archaic_exome/data/neanderthal_altai_vcfs/2014.09.29/neand_called_bases_x_indels.bbg
#
#
# ## filter files to be provided with -ir flag
# ## could be more files that should be provided
# #/tigress/AKEY/akey_vol1/home/bvernot//archaic_exome/data/chimp_from_den_epo/2013.05.02/chimp_chrAll.mapped.bbg
