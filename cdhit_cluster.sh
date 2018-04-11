#!/bin/bash

set -e

#=============================================================
# HEADER
#=============================================================

#INSTITUTION:ISCIII
#CENTRE:BU-ISCIII
#AUTHOR: Pedro J. Sola
VERSION=1.0 
#CREATED: 6 April 2018
#REVISION:
#DESCRIPTION:Script that uses cd-hit/psi-cd-hit to clusterize a FASTA file
#
#-d length of description in .clstr file, default 20. if set to 0, 
#	it takes the fasta defline and stops at first space
#-s length difference cutoff, default 0.0
#	if set to 0.9, the shorter sequences need to be
#	at least 90% length of the representative of the cluster
#-B 1 or 0, default 0, by default, sequences are stored in RAM
#	if set to 1, sequence are stored on hard drive
#	it is recommended to use -B 1 for huge databases
#-g 1 or 0, default 0
#	By cd-hit’s default algorithm, a sequence is clustered to the first
#	cluster that meet the threshold (fast mode). If set to 1, the program
#	will cluster it into the most similar cluster that meet the threshold
#	(accurate but slow mode)
#
#	PSI-CD-HIT
#-G (1/0) use global identity? default 1, sequence identity
#	calculated as total identical residues of local alignments
#	length of shorter sequence
#
#-n 5 for thresholds 0.7 ~ 1.0
#-n 4 for thresholds 0.6 ~ 0.7
#-n 3 for thresholds 0.5 ~ 0.6
#-n 2 for thresholds 0.4 ~ 0.5

#================================================================
# END_OF_HEADER
#================================================================

#SHORT USAGE RULES
#LONG USAGE FUNCTION
usage() {
	cat << EOF

Cdhit_cluster script uses cd-hit/psi-cd-hit to clusterize a FASTA file

usage : $0 <-i inputfile(FASTA)> [-o <directory>] [-n <filename>] [-c <percentage>]
		[-T <threads>] [-g group_name] [-s <int>] [-M <int>][-C <(0|1)>] [-G <(0|1)>] [-b <blast_prog>] [p] [-v] [-h]

	-i input file in sorted BAM format
	-o output directory (optional)
	-c percentage threshold to cluster, default 0.8
	-M max available memory (Mbyte), default 400
	-n file name
	-s length difference cutoff, default 0.8
	-g group name (optional). If unset, samples will be gathered in NO_GROUP group
	-p runs psi-cd-hit instead of cd-hit
	-C psi-cd-hit only: circular sequences, default 1. If set to 0 sequence is assumed lineal
	-G psi-cd-hit only: gobal identity, -G 0 only takes the first local alignment for clustering
	-b psi-cd-hit only: choose blast program, default blastn
	-T number of threads
	-v version
	-h display usage message

example: cdhit_cluster -i ecoli.fasta -c 0.9 -M 50000 -T 0
		 

EOF
}

#cd-hit-est -i $plasmidMappedFasta -o $plasmidMappedFasta"_80" -c $thresholdCdHit -n 4 -d 0 -s 0.8 -B 1 -M 50000 -T 0

#psi-cd-hit.pl -i $sample".ac.covered.fasta" -o $sample".ac.covered.fasta_80" -c $thresholdCdHit -G 1 -g 1 -prog blastn -circle 1


#================================================================
# OPTION_PROCESSING
#================================================================
#Make sure the script is executed with arguments
if [ $# = 0 ] ; then
 usage >&2
 exit 1
fi


#DECLARE FLAGS AND VARIABLES
cwd="$(pwd)"
group="NO_GROUP"
input_file="Input_file"
cluster_percentage=0.8
max_memory=400
length_cutoff=0.8
cd_hit_command=cd-hit-est
is_circle=1
global_psi_cd_hit=1
psi_cd_hit_program=blastn
word_size=0

#PARSE VARIABLE ARGUMENTS WITH getops
#common example with letters, for long options check longopts2getopts.sh
options=":i:o:c:M:n:s:g:C:G:b:T:pvh"
while getopts $options opt; do
	case $opt in
		i )
			input_file=$OPTARG
			;;
		o )
			output_dir=$OPTARG
			;;
		c )
			cluster_percentage=$OPTARG
			;;
		g)
			group=$OPTARG
			;;
		M )
			max_memory=$OPTARG
			;;
		n )
			file_name=$OPTARG
			;;
		s )			
          	length_cutoff=$OPTARG
          	;;
        p )			
          	cd_hit_command=psi-cd-hit.pl
          	;;
        C )			
          	is_circle=$OPTARG
          	;;
        G)			
          	global_psi_cd_hit=$OPTARG
          	;;
        b)			
          	psi_cd_hit_program=$OPTARG
          	;;  	
        h )
		  	usage
		  	exit 1
		  	;;
		v )
		  	echo $VERSION
		  	exit 1
		  	;;
		\?)  
			echo "Invalid Option: -$OPTARG" 1>&2
			usage
			exit 1
			;;
		: )
      		echo "Option -$OPTARG requires an argument." >&2
      		exit 1
      		;;
      	* ) 
			echo "Unimplemented option: -$OPTARG" >&2;
			exit 1
			;;

	esac
done
shift $((OPTIND-1))

#================================================================
# MAIN_BODY
#================================================================
##CHECK DEPENDENCIES, MANDATORY FIELDS, FOLDERS AND ARGUMENTS

#bash check_mandatory_files.sh $input_file

bash check_dependencies.sh cd-hit-est psi-cd-hit.pl
echo "PREV"
echo "word size= " $word_size
echo "cluster%= " $cluster_percentage


if [ "$cluster_percentage" -gt "0.7" ]; then
	echo "si"
else
	echo "NOO"
fi


<<Z
# Set word size (parameter -n for cd-hit) as author recomends
if [[ "$cluster_percentage" -gt 0.7  &&  "$cluster_percentage" -le 1 ]]; then
	word_size=5
elif [[ "$cluster_percentage" -gt "0.6"  &&  "$cluster_percentage" -le "0.7" ]]; then
	word_size=4
elif [[ "$cluster_percentage" -gt "0.5"  &&  "$cluster_percentage" -le "0.6" ]]; then
	word_size=3
elif [[ "$cluster_percentage" -gt "0.4"  &&  "$cluster_percentage" -le "0.5" ]]; then
	word_size=2
else
	echo "please introduce a valid cluster percentage value"
	exit 1
fi
Z
echo "word size= " $word_size
echo "cluster%= " $cluster_percentage
#-n 5 for thresholds 0.7 ~ 1.0
#-n 4 for thresholds 0.6 ~ 0.7
#-n 3 for thresholds 0.5 ~ 0.6
#-n 2 for thresholds 0.4 ~ 0.5
<<C
if [ ! $output_dir ]; then
	output_dir=$(dirname $input_file)
	echo "Default output directory is" $output_dir
	mkdir -p $output_dir
else
	echo "Output directory is" $output_dir
	mkdir -p $output_dir
fi

if [ ! $filename ]; then
	filename=$(basename $input_file | cut -d. -f1)
fi



if [ $positional = true ]; then 
	if [ -f $imageDir/$sample".plasmid.bedgraph" ];then \
		echo "Found a bedgraph file for sample" $sample;
		echo "Omitting bedgraph step"
	else
		echo "$(date)"
		echo "Obtaining coverage coordinates from sequences"

		bedtools genomecov -ibam $input_file -bga -max $max_coverage > $output_dir/$filename".bedgraph"

		echo "$(date)"
		echo "DONE obtaining coverage coordinates from sequences"
	fi
else


	bash check_mandatory_files.sh $database

	if [ -f $database".length" ]; then
		echo "Found length file for" $(basename $database)
		echo "Omitting length calculation"
	else
		echo "$(date)"
		echo "Creating a length file for" $(basename $database)
		bash calculate_seqlen.sh -r -i $database > $database".length"
	fi

	if [ -f $output_dir/$filename".coverage" ];then \
		echo "Found a coverage file for sample" $sample;
		echo "Omitting coverage calculation"
	else
		echo "$(date)"
		echo "Calculating coverage for every position that mapped $filename"

		bedtools genomecov -ibam $input_file -g $database".length" > $output_dir/$filename".coverage"

		echo "$(date)"
		echo "DONE Calculating coverage for every plamid that mapped $sample"
	fi
fi

C



























