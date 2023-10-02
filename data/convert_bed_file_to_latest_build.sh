#!/bin/bash

## install client if not done already
#pip install synapseclient

## log into Synapse (add your credentials here before running)
synapse login -u "" -p ""

## download relevant BED files from synapse
## JH_batch1 BED file
synapse get syn18078826

## WU_batchX BED file
synapse get syn26376829

#Install crossmap that can convert between builds
pip3 install CrossMap

#Download chain file from UCSC (https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/)

#Convert existing old bed file to txt file
cp Baits_BED_Files_AgilentV6_REVISED_S07604514_ALLBED_merged_020816_noCHR.bed Baits_BED_Files_AgilentV6_REVISED_S07604514_ALLBED_merged_020816_noCHR.txt

#Do the conversion to GRCh38 and also add "Chr" to the contigs
CrossMap.py bed --chromid l  hg19ToHg38.over.chain.gz Baits_BED_Files_AgilentV6_REVISED_S07604514_ALLBED_merged_020816_noCHR.txt Baits_BED_Files_AgilentV6_REVISED_S07604514_ALLBED_merged_020816_withChr_GRCh38.bed

#check the result
head -n 10 Baits_BED_Files_AgilentV6_REVISED_S07604514_ALLBED_merged_020816_withChr_GRCh38.bed

#Looks like it is not sorted yet. So sort the file
cat Baits_BED_Files_AgilentV6_REVISED_S07604514_ALLBED_merged_020816_withChr_GRCh38.bed | sort -k1,1V -k2,2n > Baits_BED_Files_AgilentV6_REVISED_S07604514_ALLBED_merged_020816_withChr_GRCh38_sorted.bed