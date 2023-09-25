#!/bin/bash
# Author: Mia Sungeun Lee (sungeun.lee@ec-lyon.fr)
# Description: AMOA-SEQ bash script for analysing AMOA amplicon sequencing data
# Date: 20230306
# dada_AMO.R script from DADA2 pipeline (Callahan et al. 2016)

# Set default values
exp_name=""
data_directory=""
FWD=""
REV=""
min_length=""
trunc_length=""
just_concat=""
select_length=""
num_nucleotide=""
num_nucleotide_gap=""
organism=""

function usage {
  echo "
  ############################################################################################
  #     __      ___    ___     ____        __               _______    _______    _______    #
  #    /  \    |   \  /   |  /  __  \     /  \             / _______\ |  _____|  / ______ \  #
  #   / /\ \   | |\ \/ /| | / /    \ \   / /\ \     ____  | |_______  | |_____  | |______| | #
  #  / /__\ \  | | \__/ | || |     | |  / /__\ \   |____|  \______  \ |  _____|  \_______  | #
  # / ______ \ | |      | | \ \ __ / / / ______ \          _______| | | |_____          |  | #
  #/_/      \_\|_|      |_|  \ ____ / /_/      \_\         \_______/  |_______|         |__| #
  ############################################################################################

  ";
echo "Exploring the AMOA diversity";
echo " ";
echo "Usage: $0 [-e OUTPUT_FILE_NAME] [-i DATA_DIRECTORY] [-f FORWARD_READ] [-r REVERSE_READ] [-m MIN_LENGTH] [-l TRUNC_LENGTH] [-c JUST_CONCAT] [-t SELECT_LENGTH] [-n NUM_NUCLEOTIDE] [-o ORGANISM] [-h]"
}

# Parse command-line options
while getopts "e:i:f:r:m:l:c:t:n:o:h" flag; do
case "${flag}" in
e) exp_name="${OPTARG}" ;;
i) data_directory="${OPTARG}" ;;
f) FWD="${OPTARG}" ;;
r) REV="${OPTARG}" ;;
m) min_length="${OPTARG}" ;;
l) trunc_length="${OPTARG}" ;;
c) just_concat="${OPTARG}" ;;
t) select_length="${OPTARG}" ;;
n) num_nucleotide="${OPTARG}" ;;
o) organism="${OPTARG}" ;;
h) usage; exit 0 ;;
*) usage; exit 1 ;;
esac
done

while getopts "e:i:f:r:m:l:c:t:n:o:" flag; do
    case "${flag}" in
        e) exp_name="${OPTARG}" ;;
        i) data_directory="${OPTARG}" ;;
        f) FWD="${OPTARG}" ;;
        r) REV="${OPTARG}" ;;
        m) min_length="${OPTARG}" ;;
        l) trunc_length="${OPTARG}" ;;
        c) just_concat="${OPTARG}" ;;
        t) select_length="${OPTARG}" ;;
        n) num_nucleotide="${OPTARG}" ;;
        o) organism="${OPTARG}" ;;
    esac
done


for opt in "${required_opts[@]}"; do
if [[ -z "${!opt}" ]]; then
echo "ERROR: ${opt} is required."
usage
exit 1
fi
done


echo "
############################################################################################
#     __      ___    ___     ____        __               _______    _______    _______    #
#    /  \    |   \  /   |  /  __  \     /  \             / _______\ |  _____|  / ______ \  #
#   / /\ \   | |\ \/ /| | / /    \ \   / /\ \     ____  | |_______  | |_____  | |______| | #
#  / /__\ \  | | \__/ | || |     | |  / /__\ \   |____|  \______  \ |  _____|  \_______  | #
# / ______ \ | |      | | \ \ __ / / / ______ \          _______| | | |_____          |  | #
#/_/      \_\|_|      |_|  \ ____ / /_/      \_\         \_______/  |_______|         |__| #
############################################################################################

";
echo "============================================================================================";
echo "#Ouput_directory_name: $exp_name";
echo "#Data_directory: $data_directory";
echo "#Foward_primer: $FWD";
echo "#Reverse_primer: $REV";
echo "#Min_length: $min_length";
echo "#Trunc_length: $trunc_length";
echo "#Just_concat: $just_concat";
echo "#Trim_length: $select_length";
echo "#Num_nucleotide: $num_nucleotide";
echo "#Type_organism: $organism";
workng_directory=$(pwd)
echo "#workng_directory: $workng_directory";
echo "============================================================================================";

#####################################################
echo "============================================================================================";
echo "### STEP 0. Making AMOA database in $exp_name working directory ###"
mkdir $exp_name
cp dada_AMO.R ./$exp_name/
gunzip AMO_database.faa.gz
diamond makedb --in AMO_database.faa  --db AMO.dmnd
diamond makedb --in ref.$organism.amoA.faa --db ref.$organism.amoA.dmnd
mv AMO.dmnd ./$exp_name/
mv ref.$organism.amoA.dmnd ./$exp_name/
cp ref.$organism.amoA.faa ./$exp_name/
cp *.py ./$exp_name/
cd $exp_name
echo "### STEP 0. Done ###"
echo "============================================================================================";
#####################################################
echo "============================================================================================";
echo "### STEP 1. DADA2 analysis in $exp_name working directory ###"
# Call R script with command-line options as arguments
Rscript dada_AMO.R "${exp_name}" "${data_directory}" "${FWD}" "${REV}" "${min_length}" "${trunc_length}" "${just_concat}" "${select_length}" "${num_nucleotide}" "${num_nucleotide_gap}" "${organism}"
mv out.ASVs.fa out.DADA2.$organism.ASVs.fa
mv out.ASVs.counts.tsv out.DADA2.$organism.ASVs.counts.tsv
mv track-summary.tsv out.DADA2.$organism.ASVs.track-summary.tsv
echo "### STEP 1. DADA2 analysis done ###"
echo "============================================================================================";
#####################################################

#####################################################
### Although DADA2 tool pipeline removed the ambigous sequences using its error model, it can produce several ambigous sequences depending on the datasets.
### Thus, the ambigous sequences were removed by
###   - checking expected amplicon sizes
###   - blasting several AMO databases
###   - checking conserved regions of manually curated AMO databases and several AMO datasets at protein sequence level using final.$organism.py script
###
#####################################################
#####################################################
echo "============================================================================================";
echo "### STEP 2. Select the ASV sequences according to expected amplicon size ###"
seqkit seq -m $select_length -M $select_length out.DADA2.$organism.ASVs.fa > out.DADA2.correct-size.$organism.ASVs.fa
grep ">" out.DADA2.correct-size.$organism.ASVs.fa | sed 's/>//g' > ASV-ID
awk 'FNR==NR{a[$1];next} $1 in a{print}' ASV-ID out.DADA2.$organism.ASVs.counts.tsv > out.DADA2.correct-size.$organism.ASVs.counts.tsv
awk 'NR==1 {print}' out.DADA2.$organism.ASVs.counts.tsv > ID
cat ID out.DADA2.correct-size.$organism.ASVs.counts.tsv > tmp && mv tmp out.DADA2.correct-size.$organism.ASVs.counts.tsv
echo "### STEP 2. Correct size of the DADA2 generated ASV sequences was selected ###"
echo "============================================================================================";

echo "============================================================================================";
echo "### STEP 3. Correct AmoA sequence curation using the AMO databases and conserved protein sequence ###"
diamond blastx --db AMO.dmnd  --query out.DADA2.$organism.ASVs.fa --out diamond.output.DADA2.$organism.ASVs.tsv --evalue 0.00001  --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle qlen slen qcovhsp
awk '!x[$1]++' diamond.output.DADA2.$organism.ASVs.tsv > besthit.diamond.output.DADA2.$organism.ASVs.tsv
awk '{print $1}' diamond.output.DADA2.$organism.ASVs.tsv | sort -u > Annotated-ASV-ID
seqkit grep -n -f Annotated-ASV-ID out.DADA2.correct-size.$organism.ASVs.fa > annotated.DADA2.$organism.ASVs.fa


# In order to correct translation of AMO sequences, remove one 1st nucleotide and one N were removed for AOA and two 1st nucleotides were removed for AOB
sed 's/NNNNNNNNNN/NNNNNNNNN/g' annotated.DADA2.AOA.ASVs.fa > tmp && mv tmp annotated.DADA2.AOA.ASVs.fa
seqkit seq -w 0 annotated.DADA2.$organism.ASVs.fa > tmp && mv tmp annotated.DADA2.$organism.ASVs.fa
awk '/^>/ {print $0} /^[^>]/ {print substr($0,'$num_nucleotide')}' annotated.DADA2.$organism.ASVs.fa > tmp && mv tmp annotated.DADA2.$organism.ASVs.fa
seqkit translate -f 1 annotated.DADA2.$organism.ASVs.fa > annotated.DADA2.$organism.ASVs.faa
### removing the ambigous sequences using final.$organism.py python scripts:
python final.$organism.py -i annotated.DADA2.$organism.ASVs.faa -o $organism.err-seq-ID
sed 's/X//g' annotated.DADA2.$organism.ASVs.faa > tmp && mv tmp annotated.DADA2.$organism.ASVs.faa
seqkit grep -v -n -f $organism.err-seq-ID annotated.DADA2.$organism.ASVs.fa > AMOA-SEQ-curated.$organism.ASVs.fa
seqkit grep -v -n -f $organism.err-seq-ID annotated.DADA2.$organism.ASVs.faa > AMOA-SEQ-curated.$organism.ASVs.faa
seqkit grep -n -f $organism.err-seq-ID annotated.DADA2.$organism.ASVs.fa > AMOA-SEQ-ambigous.$organism.ASVs.fa
seqkit grep -n -f $organism.err-seq-ID annotated.DADA2.$organism.ASVs.faa > AMOA-SEQ-ambigous.$organism.ASVs.faa

grep ">" AMOA-SEQ-curated.$organism.ASVs.fa | sed 's/>//g' > AMOA-SEQ-curated.$organism.ASVs-ID
awk 'FNR==NR{a[$1];next} $1 in a{print}' AMOA-SEQ-curated.$organism.ASVs-ID out.DADA2.$organism.ASVs.counts.tsv > AMOA-SEQ-curated.$organism.ASVs.counts.tsv
awk 'NR==1 {print}' out.DADA2.$organism.ASVs.counts.tsv > ID
cat ID AMOA-SEQ-curated.$organism.ASVs.counts.tsv > tmp && mv tmp AMOA-SEQ-curated.$organism.ASVs.counts.tsv
echo "### STEP 3. AMOA-SEQ ASV curation done ###"
echo "============================================================================================";

#####################################################
echo "============================================================================================";
echo "### STEP 4. Comparing the AMOA-SEQ curated ASV sequences to curated AMOA databases ###"
# Annotating curated ASVs using manually curated databases
diamond blastx --db ref.$organism.amoA.dmnd  --query AMOA-SEQ-curated.$organism.ASVs.fa --out diamond.output.curateddb.AMOA-SEQ-curated.$organism.ASVs.tsv --evalue 0.00001  --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle qlen slen qcovhsp
awk '!x[$1]++' diamond.output.curateddb.AMOA-SEQ-curated.$organism.ASVs.tsv > besthit.diamond.output.curateddb.AMOA-SEQ-curated.$organism.ASVs.tsv
echo "### STEP 4. AMOA-SEQ ASV annotation using curated AMOA databases done ###"
echo "============================================================================================";
#####################################################

#####################################################
echo "============================================================================================";
echo "### STEP 5. ASV clustering into OTUs and generating OTU count table ###"
cd-hit-est -i AMOA-SEQ-curated.$organism.ASVs.fa -o AMOA-SEQ.$organism.OTUs.fa -c 0.97 -n 5
python ASV-to-OTU.py -i AMOA-SEQ.$organism.OTUs.fa.clstr -o OTU_ASV_ID.txt
python OTU-table.py -i AMOA-SEQ-curated.$organism.ASVs.counts.tsv -t OTU_ASV_ID.txt -o AMOA-SEQ.$organism.OTUs.counts.tsv
awk '{print $2, "\t", $1}' OTU_ASV_ID.txt | sort -u > ASV_OTU_ID.txt
sed 's/ //g' ASV_OTU_ID.txt > tmp && mv tmp ASV_OTU_ID.txt
sed 's/ //g' ASV_OTU_ID.txt > tmp && mv tmp ASV_OTU_ID.txt
awk 'BEGIN { FS="\t" }
     NR==FNR { map[$1]=$2; next }
     /^>/ { print ">" map[substr($0,2)]; next }
     { print }' ASV_OTU_ID.txt AMOA-SEQ.$organism.OTUs.fa > tmp && mv tmp AMOA-SEQ.$organism.OTUs.fa
diamond blastx --db ref.$organism.amoA.dmnd --query AMOA-SEQ.$organism.OTUs.fa --out diamond.output.curateddb.AMOA-SEQ.$organism.OTUs.tsv --evalue 0.00001  --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle qlen slen qcovhsp
awk '!x[$1]++' diamond.output.curateddb.AMOA-SEQ.$organism.OTUs.tsv > besthit.diamond.output.curateddb.AMOA-SEQ.$organism.OTUs.tsv
awk '{print $1, "\t", $13}' besthit.diamond.output.curateddb.AMOA-SEQ.$organism.OTUs.tsv > ID-Taxa
grep ">" AMOA-SEQ.$organism.OTUs.fa | sed 's/>//' > OTU-ID
awk 'FNR==NR{a[$1];next} $1 in a{print; delete a[$1]} END{for (i in a) print i, "NA"}' OTU-ID ID-Taxa > AMOA-SEQ.$organism.OTUs.taxa.tsv
echo "### STEP 5. AMOA-SEQ AMO-OTU count table generated and annotation done ###"
echo "============================================================================================";
#####################################################

#####################################################
echo "============================================================================================";
echo "### STEP 6. translating the ASV sequences to PSV sequences ###"
cd-hit -i AMOA-SEQ-curated.$organism.ASVs.faa -o AMOA-SEQ.$organism.PSVs.faa -c 1 -n 5
python ASV-to-PSV.py -i AMOA-SEQ.$organism.PSVs.faa.clstr -o PSV_ASV_ID.txt
python PSV-table.py -i AMOA-SEQ-curated.$organism.ASVs.counts.tsv -t PSV_ASV_ID.txt -o AMOA-SEQ.PSVs.counts.tsv
awk '{print $2, "\t", $1}' PSV_ASV_ID.txt | sort -u > ASV_PSV_ID.txt
sed 's/ //g' ASV_PSV_ID.txt > tmp && mv tmp ASV_PSV_ID.txt
sed 's/ //g' ASV_PSV_ID.txt > tmp && mv tmp ASV_PSV_ID.txt
awk 'BEGIN { FS="\t" }
     NR==FNR { map[$1]=$2; next }
     /^>/ { print ">" map[substr($0,2)]; next }
     { print }' ASV_PSV_ID.txt AMOA-SEQ.$organism.PSVs.faa > tmp && mv tmp AMOA-SEQ.$organism.PSVs.faa
echo "### STEP 6. translating the ASV sequences to PSV sequences and dereplication of PSVs, done ###"
echo "============================================================================================";
#####################################################

#####################################################
echo "============================================================================================";
echo "### STEP 7. Annotating the PSV sequences against curated AMOA database using BLASTp ###"
blastp -query AMOA-SEQ.$organism.PSVs.faa -subject ref.$organism.amoA.faa -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore salltitles' -out blastp.output.AMOA-SEQ.$organism.PSVs.tsv -num_threads 16 -evalue 0.00001
awk '!x[$1]++' blastp.output.AMOA-SEQ.$organism.PSVs.tsv  > besthit.blastp.output.AMOA-SEQ.$organism.PSVs.tsv
echo "### STEP 7. AMOA-SEQ AMO-PSV Annotation done ###"
echo "============================================================================================";
#####################################################
echo "============================================================================================";
echo "### STEP 8. Aligning of the PSV sequences and curated AMOA sequences for generating phylogenetic tree ###"
cat AMOA-SEQ.$organism.PSVs.faa ref.$organism.amoA.faa > tree.$organism.faa
muscle -super5 tree.$organism.faa -output tree.$organism.afa
trimal -in tree.$organism.afa -out tree.$organism.trim.afa -nogaps
FastTree tree.$organism.trim.afa > tree.$organism.nwk
echo "### STEP 8. Phylogenetic tree generated ### "
echo "============================================================================================";
#####################################################
rm ID ASV-ID Annotated-ASV-ID ID-Taxa OTU-ID ASV_OTU_ID.txt
mkdir $organism.ASV-analysis
mkdir $organism.PSV-analysis
mkdir $organism.Phylogenetic-analysis
mkdir $organism.OTU-analysis
rm *dmnd *R
rm *.py
mv *ASVs* ./$organism.ASV-analysis
mv *PSV* ./$organism.PSV-analysis
mv *tree* ./$organism.Phylogenetic-analysis
mv *OTUs* ./$organism.OTU-analysis
#####################################################
