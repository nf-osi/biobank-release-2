# 2024-01-23

# originally used a c5.xlarge instance with 1000 GB of memory; looks like now the equivalent available instance is c7a.xlarge
# CostCenter	Other / 000001
# CostCenterOther	NTAP NF Amendment 6 / 301102

# switch to ec2-user (I found this necessary in order to be able to work in the instance, not sure if that is the case for others)
sudo su - ec2-user

# install compiler
sudo yum groupinstall "Development Tools"
sudo yum update

# install Synapse client
pip3 install --upgrade pip
pip install synapseclient

# install mamba
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(uname)-$(uname -m).sh"
bash Mambaforge-$(uname)-$(uname -m).sh

# install samtools, bcftools, and somalier
mamba create -n ase -c bioconda samtools
mamba activate ase
mamba install -c bioconda bcftools
mamba install somalier


# log in to Synapse, then get files
# get reference HG38 fasta files
synapse get syn51706561 --downloadLocation '/home/ec2-user/'
synapse get syn51706450 --downloadLocation '/home/ec2-user/'


# get data that has already been extracted from cram/bam files with Somalier (from dataset I made on Synapse)
mkdir extracted
synapse get -q "SELECT * FROM syn53425758" --downloadLocation '/home/ec2-user/extracted/'


# get additional bam/cram files you want to add to the analysis (replace $SYNID with relevant Synapse ID)
mkdir aligned_data
synapse get $SYNID --downloadLocation '/home/ec2-user/aligned_data/'


# use somalier to extract polymorphic sites
for f in aligned_data/*.bam
do
	./somalier extract -d extracted --sites somalier_sites/sites.hg38.vcf.gz -f /home/ec2-user/reference/Homo_sapiens_assembly38.fasta.gz $f
done

for f in aligned_data/*.cram
do
	./somalier extract -d extracted --sites somalier_sites/sites.hg38.vcf.gz -f /home/ec2-user/reference/Homo_sapiens_assembly38.fasta.gz $f
done

# get relatedness from somalier
./somalier relate extracted/*.somalier




########################################################################################
# from here the rest is documented in somalier_manuscript_text_v1.Rmd (https://github.com/nf-osi/biobank-release-2/blob/dev/manuscript/somalier_manuscript_text_v1.Rmd)
# but I will also add here for reference, with some extra steps you may need to take to remove unwanted samples/fix sample IDs from other projects

# fix ids from other project
cat somalier_relatedness.txt | awk '{OFS="\t"; if($1=="2-009_pt_rna_jhk"){$1="JH-2-009-HH11A-C47BD"} else if($1=="2-055_pt_rna"){$1="JH-2-055-FC9GH-CG8EA"} else if($1=="2-079_pt"){$1="JH-2-079-CB92G-FG226"}; if($2=="2-009_pt_rna_jhk"){$2="JH-2-009-HH11A-C47BD"} else if($2=="2-055_pt_rna"){$2="JH-2-055-FC9GH-CG8EA"} else if($2=="2-079_pt"){$2="JH-2-079-CB92G-FG226"}; print }' > somalier_relatedness.id.txt
cat somalier.samples.tsv | awk '{OFS="\t"; if($1=="2-009_pt_rna_jhk"){$1="JH-2-009-HH11A-C47BD"} else if($1=="2-055_pt_rna"){$1="JH-2-055-FC9GH-CG8EA"} else if($1=="2-079_pt"){$1="JH-2-079-CB92G-FG226"}; if($2=="2-009_pt_rna_jhk"){$2="JH-2-009-HH11A-C47BD"} else if($2=="2-055_pt_rna"){$2="JH-2-055-FC9GH-CG8EA"} else if($2=="2-079_pt"){$2="JH-2-079-CB92G-FG226"}; print }' > somalier.samples.id.tsv
cat somalier.pairs.tsv | awk '{OFS="\t"; if($1=="2-009_pt_rna_jhk"){$1="JH-2-009-HH11A-C47BD"} else if($1=="2-055_pt_rna"){$1="JH-2-055-FC9GH-CG8EA"} else if($1=="2-079_pt"){$1="JH-2-079-CB92G-FG226"}; if($2=="2-009_pt_rna_jhk"){$2="JH-2-009-HH11A-C47BD"} else if($2=="2-055_pt_rna"){$2="JH-2-055-FC9GH-CG8EA"} else if($2=="2-079_pt"){$2="JH-2-079-CB92G-FG226"}; print }' > somalier.pairs.id.tsv

# fix other ids (some of the extractions can add on an extra sample identifier - not sure what caused this)
cat somalier.pairs.id.tsv | awk '{OFS="\t"; if($1~/\_JH/){split($1,a,"_"); $1=a[2]}; if($2~/\_JH/){split($2,b,"_"); $2=b[2]}; print }' > somalier.pairs.fixedid.tsv
cat somalier.samples.id.tsv | awk '{OFS="\t"; if($1~/\_JH/){split($1,a,"_"); $1=a[2]}; if($2~/\_JH/){split($2,b,"_"); $2=b[2]}; print }' > somalier.samples.fixedid.tsv
cat somalier_relatedness.id.txt | awk '{OFS="\t"; if($1~/\_JH/){split($1,a,"_"); $1=a[2]}; if($2~/\_JH/){split($2,b,"_"); $2=b[2]}; print }' > somalier_relatedness.fixedid.txt

# extract needed statistics from somalier output and determined if samples come from same individual
cut -f 1-3,5-6 somalier.pairs.fixedid.tsv \
	| awk '{OFS="\t"; if($1=="sample_a"){print $0,"paired"} else{split($1,a,"-"); split($2,b,"-"); if(a[3]==b[3]){print $0,"yes"} else{print $0,"no"}}}' \
	> somalier_relatedness.txt

grep -v 'JH-2-084' somalier_relatedness.txt | grep -v 'JH-2-106' | grep -v 'JH-2-009-2578C-A' | grep -v 'JH−2−009−2578C−9FG16' > tmp.txt
mv tmp.txt somalier_relatedness.txt

# extract relevant information from sample map
awk -F "\t" '{OFS="\t"; print $12"-"$13,$5,$15}' Biobank-samples-map-clinical-102523.txt | tr " " "_" > biobank_info.txt

# combine somalier statistics and sample map information, keep only RNA to WES comparisons, and add column with type of tissues compared
zjoin -a somalier_relatedness.txt -b biobank_info.txt -1 1 -2 1 \
	| cut -f 1-6,8- \
	| zjoin -a stdin -b biobank_info.txt -1 2 -2 1 \
	| cut -f 1-8,10- \
	| awk '($8~/RNA/ && $10~/Exome/) || ($8~/Exome/ && $10~/RNA/)' \
	| awk '{OFS="\t"; if($7=="cell_line" || $9=="cell_line"){print $0,"cell_line"} else if($7=="xenograft_passage" || $9=="xenograft_passage"){print $0,"xenograft"} else if($7=="blood" || $9=="blood"){print $0,"tumor/normal"} else{print $0,"tumor"} }' \
	> somalier_relatedness.rna_to_WES.comparison_type.tsv

# Convert relatedness statistics from tabular to matrix format
cut -f 1-3,8,10 somalier_relatedness.rna_to_WES.comparison_type.tsv > rna_to_WES.type.txt
awk '{OFS="\t"; if($4=="Whole_Exome_Sequencing"){print $2,$1,$3} else{print $1,$2,$3}}' rna_to_WES.type.txt > rna_to_WES.ordered.txt
cat rna_to_WES.ordered.txt | grep -v JH-2-009-2578C-A > rna_to_WES.ordered.removed_samples.txt
sort -k1,1 -k2,2 rna_to_WES.ordered.removed_samples.txt > rna_to_WES.ordered.removed_samples.sorted.txt

# note that ths python script can be found here: https://github.com/nf-osi/biobank-release-2/blob/dev/manuscript/make_matrix_rna_vs_WES.py
python matrix/make_matrix_rna_vs_WES.py -i rna_to_WES.ordered.removed_samples.sorted.txt | awk '!/^[[:space:]]*$/' | sort -k1,1 | awk -F "\t" '{OFS="\t"; if(NR==1){$1="sample"; print} else{print}}' > relatedness_matrix.rna_to_WES.txt
# re-order matrix to have all individuals together, and all WES and RNA-seq sampels together within individual's sample cluster. Using biobank info file generated in Bash code block above.
zjoin -a biobank_info.txt -b relatedness_matrix.rna_to_WES.txt -1 1 -2 1 \
	| awk '{OFS="\t"; split($1,a,"-"); print a[3],$0}' \
	| sort -k1,1 -k4,4 \
	| cut -f 5- \
	| cat <(head -n 1 relatedness_matrix.rna_to_WES.txt) - \
	> relatedness_matrix.rna_to_WES.ordered.txt


cut -f 1 relatedness_matrix.rna_to_WES.ordered.txt | grep -v sample | zjoin -a stdin -b biobank_info.txt -1 1 -2 1 | cut -f 2- > rna_to_WES.row_info.txt
head -n 1 relatedness_matrix.rna_to_WES.ordered.txt | tr '\t' '\n' | grep -v sample | zjoin -a stdin -b biobank_info.txt -1 1 -2 1 | cut -f 2- > rna_to_WES.col_info.txt
awk '{OFS="\t"; if($2=="blood" || $2=="normal_tissue"){$2="normal"}; print }' rna_to_WES.row_info.txt > tmp.txt
mv tmp.txt rna_to_WES.row_info.txt



# the above should give you all of the files you need to re-generate the plots with the code from here: https://github.com/nf-osi/biobank-release-2/blob/dev/manuscript/somalier_manuscript_text_v1.Rmd

