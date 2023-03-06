# AMOA-SEQ
````
############################################################################################
#     __      ___    ___     ____        __               _______    _______    _______    #
#    /  \    |   \  /   |  /  __  \     /  \             / _______\ |  _____|  / ______ \  #
#   / /\ \   | |\ \/ /| | / /    \ \   / /\ \     ____  | |_______  | |_____  | |______| | #
#  / /__\ \  | | \__/ | || |     | |  / /__\ \   |____|  \______  \ |  _____|  \_______  | #
# / ______ \ | |      | | \ \ __ / / / ______ \          _______| | | |_____          |  | #
#/_/      \_\|_|      |_|  \ ____ / /_/      \_\         \_______/  |_______|         |__| #
############################################################################################
````
AMOA-SEQ is a bash pipeline using an expert-guided apparoch to analyze diverse AMOA sequences.
AMOA-SEQ utilizes non-redundant AMOA sequences from NCBI as well as curated archaeal, bacterial and comammox AMOA sequences from Elves et al. 2018, Lee et al. 2023 (this study), Palomo et al 2022 for taxonomic and phylogenetic analyses. 

• work with AMOA amplicon sequencing (Miseq reagent kit V2 (500-cycles), 2 x 250 bp) with archaeal, bacterial and comammox amoA gene amplified using the primer CrenamoA23F/CrenamoA616R (Tourna et al. 2008), AmoA1F/AmoA2R (Rotthauwe et al. 1997), and comamoA-F/comamoA-R (Zhao et al. 2018)

• For archaeal amoA amplicons that are non-overlapping reads (> 600 bp), the ‘gap’ pipeline was used (Aigle et al. 2019).

• apply DADA2 tool to generate amplicon sequence variants (ASVs), enabling a higher-resolution alternative to conventional operational taxonomic units (OTUs), which records the exact number of times a specific amplicon sequence variant is observed in each sample. Additionally, DADA2 incorporates quality information into its error model, making the algorithm more resilient to lower quality sequences (Callahan et al. 2016).

• provide phylogenetic tree of AMOA sequences using a curated AMOA database from Elves et al. 2018, Palomo et al 2022, Lee et al. 2023.

See more details in the publication (Lee et al. 2023 hopefully publish soon!).

## Installation with Conda
````
conda create -y -n AMOA-SEQ -c bioconda seqkit fasttree muscle blast diamond trimal diamond=0.9.19 cd-hit cutadapt
````

## DADA2 tool installation 
````
https://benjjneb.github.io/dada2/dada-installation.html
````

## Clone AMOA-SEQ directory
````
git clone https://github.com/miasungeunlee/AMOA-SEQ.git
cd AMOA-SEQ
chmod u+x AMOA-SEQ.sh # make the script executable
````

## AMOA sequence databases
•	AMO_database.faa: all AMOA sequences downloaded from JGI IMG (https://img.jgi.doe.gov/) and NCBI (https://www.ncbi.nlm.nih.gov/). List of accessions and detail of the gene set were shown in AMO_database.tsv 

•	ref.AOA.amoA.faa : curated archaeal amoA sequences with defined lineage from Elves et al. 2018

•	ref.COM.amoA.faa : curated archaeal amoA sequences from Palomo et al. 2022

•	ref.AOB.amoA.faa  : bacterial amoA sequences were downloaded from JGI IMG site and curated (Lee et al. 2023) 

## Quick run
````
source activate AMOA-SEQ
sh AMOA-SEQ.sh [-h help] [-e output_directory_name] [-i fastq_directory] [-f forward_primer] [-r reverse_primer] [-m minimum_read_length] [-l truncation_read_length] [-c just_concatenating_option] [-t expected merged sequence_length] [-t expected_merged_sequence_length] [-n number_nucleotide] [-o AO_type]

### option variable explanation ###
-h: help
-e: output directory name (e.g. AOA-output)
-i: fastq.gz file path (e.g. /home/ampere/slee/COMICON-Projet-2022/TEST-AOA) 
# fastq file must end with either _R1_001.fastq.gz or _R2_001.fastq.gz pattern (directly from MiSeq sequencing output fastq #
-f: forward primer (e.g. ATGGTCTGGCTWAGACG for AOA forward primer)
-r: reverse primer (e.g. GCCATCCATCTGTATGTCCA for AOA reverse primer)
-m: minimum read length (e.g. 200 bp for AOA)
-l: truncation read length (e.g. 200 bp for AOA)
-c: TRUE for AOA AMOA, just concatenating forward and reverse reads, for other AMOA, merging forward and reverse reads are possible, use FALSE option
-t: expected merged sequence length (410 bp for AOA)
-n: number of nucleotides to be removed prior to correct translation 
-o: organism; it can be either AOA or AOB or COM depending on your dataset


````

### Example of the run the AMOA-SEQ.sh for AOA, AOB and Comammox AMOA amplicon sequencing:
````
sh AMOA-SEQ.sh -e AOA-output -i /home/ampere/slee/COMICON-Projet-2022/AOA -f ATGGTCTGGCTWAGACG -r GCCATCCATCTGTATGTCCA -m 200 -l 200 -c TRUE -t 410 -n 2 -o AOA
sh AMOA-SEQ.sh -e AOB-output -i /home/ampere/slee/COMICON-Projet-2022/AOB -f GGGGTTTCTACTGGTGGT -r CCCCTCKGSAAAGCCTTCTTC -m 231 -l 250 -c FALSE -t 452 -n 3 -o AOB
sh AMOA-SEQ.sh -e COM-output -i /home/ampere/slee/COMICON-Projet-2022/COM -f AGGNGAYTGGGAYTTCTGG -r CGGACAWABRTGAABCCCAT -m 204 -l 250 -c FALSE -t 396 -n 3 -o COM
````

### Output directory and files
````
{organism}.ASV-analysis

•	track.{organism}.summary.tsv: DADA2 output summary containing quality control, denoising, number of merged sequences, number of chimeras in each sample. 
•	out.{organism}.ASVs.fa: generated amplicon sequence variants (however, some ASVs may be not real AMOA sequences due to sequencing error, not recommended to directly use these ASV sequences)
•	out.{organism}.ASVs.counts.tsv: ASV count table from different samples
•	correct.{organism}.ASVs.fa: selected ASVs according to expected amplicon size (option -t) 
•	correct.{organism}.ASVs.counts.tsv: selected ASV count table from different samples
•	diamond.output.{organism}.ASVs.tsv: annotation of ASVs using total AMOA database
•	besthit.diamond.output.{organism}.ASVs.tsv: besthit of annoated ASVs using total AMOA database
•	diamond.output.curateddb.{organism}.ASVs.tsv: annotation of ASVs using curated AMOA database  
•	besthit.diamond.output.curateddb.{organism}.ASVs.tsv: besthit of annoated ASVs using curated AMOA database  
•	annotated.{organism}.ASVs.fa: ASVs matched to AMOA database (those are confident and genuine AMOA sequences). 
•	annotated.{organism}.ASVs.counts.tsv annotated ASV count table from different samples (recommended for alpha-diversity analysis (e.g. diversity indexes)) 

{organism}.PSV-analysis

•	annotated.{organism}.ASVs.faa: translate ASV sequence to amino acid sequences
•	{organism}.PSV.faa.clstr: clustering of translated ASV sequences with 100% identity
•	{organism}.PSV.faa: Unique protein sequence variants
•	blastp.output.{organism}.PSVs.tsv : (recommended for beta-diversity analysis (e.g. phylogenetic tree)): annotation of PSVs using curated AMOA database
•	besthit.blastp.output.{organism}.PSVs.tsv : besthit of annotated PSVs using curated AMOA database 

{organism}.phylogeny-analysis

•	tree.{organism}.faa: PSV sequences + curated AMOA sequence for phylogeny analysis
•	tree.{organism}.afa: AMOA-aligments for each PSV and curated AMOA sequences using MUSCLE with -super5 option
•	tree.{organism}.trim.afa: Ambiguous regions and gaps removed using trimal with -nogaps option
•	tree.{organism}.nwk: AMOA FastTREE file in Newick format  
