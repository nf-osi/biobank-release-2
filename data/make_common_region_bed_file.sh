
## first run convert_bed_file_to_latest_build.sh

## Find common regions between JH_batch1 BED file and WU_batchX BED file, with atleast 50% overlap, and retain the positions from the JH_batch1 BED file
bedtools intersect -a Baits_BED_Files_AgilentV6_REVISED_S07604514_ALLBED_merged_020816_withChr_GRCh38_sorted.bed -b xgen-exome-research-panel-v2-probes-hg3862a5791532796e2eaa53ff00001c1b3c.bed -f 0.5 -wa > JHbatch1_WUbatches_BED_GRCh38_normalBEDposition_50percentoverlap_noreciprocal.bed
