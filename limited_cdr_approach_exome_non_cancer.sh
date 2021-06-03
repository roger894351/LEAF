cdr=$1

# argument1 test can be anything remove this later
python /rsrch3/home/epi/jchen15/apps/roger_script/VAAST-CASM/real_data/weight/ov/test_limited_cdr.py test $cdr 2> test_limited_cdr.err
#bedtools intersect -a /rsrch3/home/epi/jchen15/database/gnomAD/script/AF_exome_bed -b temp_001.bed > test_output2.af
#use a non_cancer external allele frequency

#bedtools intersect -a /rsrch3/home/epi/jchen15/database/gnomAD/script/all_non_cancer_exome.af.bed -b temp_001.bed > test_output2.af
#exome control only bed bedtools has bugs that casn do intercept correctly even SNP
bedtools intersect -a /rsrch3/home/epi/jchen15/database/gnomAD/script/all_exome_non_cancer.af.bed -b temp_001.bed > test_output2.af
#use exome_controls
#/rsrch3/home/epi/jchen15/apps/roger_script/VAAST-CASM/real_data/weight/ov/limited_cdr_approach_exome_controls.sh


#gnomAD_control_only
#gnomad_controls_exome_region.af.bed

#all_non_cancer_exome.af.bed
#python /rsrch3/home/epi/jchen15/database/ExAC_data/read_MAF_counts_chr.py /rsrch3/home/epi/jchen15/database/gnomAD/script/all_exome_non_cancer.af.bed
python /rsrch3/home/epi/jchen15/database/ExAC_data/read_MAF_counts_chr.py test_output2.af

#rm temp_001.bed
#rm test_output2.af
