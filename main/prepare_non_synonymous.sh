external_af=$1 #the external bed file
external_control_annotation=$2


grep -w nonsynonymous $external_control_annotation|sort -k 1,2|uniq > nonsynonymous_bed_nononsense.txt
cut -f 1-3 nonsynonymous_bed_nononsense.txt >nonsynonymous_bed_nononsense.bed
bedtools intersect -a $external_af -b nonsynonymous_bed_nononsense.bed >all_exome_non_cancer.af_nonsynonymous_nonon
sense.bed
sort -k 1,2  all_exome_non_cancer.af_nonsynonymous_nononsense.bed|uniq >af_nonsynonymous_nononsense_uniq.bed
python site_frequency_spectrum.py af_nonsynonymous_nononsense_uniq.bed
