#!/bin/bash

# Exit immediately if a pipeline, which may consist of a single simple command, a list, 
#or a compound command returns a non-zero status: If errors are not handled by user
set -e
# Treat unset variables and parameters other than the special parameters ‘@’ or ‘*’ as an error when performing parameter expansion.
# An error message will be written to the standard error, and a non-interactive shell will exit
set -u
#Print everything as if it were executed, after substitution and expansion is applied: Debug|log option
set -x

#=============================================================
# HEADER
#=============================================================

#INSTITUTION:ISCIII
#CENTRE:BU-ISCIII
#AUTHOR: Pedro J. Sola
VERSION=1.0 
#CREATED: 15 March 2018
#REVISION:
#		19 March 2018: Complete usage info
#		19 March 2018: Check mandatory files. folders and variables
#DESCRIPTION:Script that index a database and map a supplied pair-end sequences
#TODO
#	-Handle files extensions for bowtie, now is fastq by default
#================================================================
# END_OF_HEADER
#================================================================

#SHORT USAGE RULES
#LONG USAGE FUNCTION
usage() {
	cat << EOF

Bowtie_mapper script index a database and map a supplied pair-end sequences

usage : $0 [-i <inputfile>] [-o <directory>] <-d database(fasta)> <-s sample_name> <-1 R1> <-2 R2> 
		[-g group_name] [-f <int>] [-T <int>] [-a] [-v] [-h]

	-i input directory (optional)
	-o output directory (optional)
	-d database to map (.fasta)
	-s sample name
	-g group name (optional). If unset, samples will be gathered in NO_GROUP group
	-1 reads corresponding to paired-end R1
	-2 reads corresponding to paired-end R2
	-f offrate index for bowtie_build (optional). Default value 1. for quicker indexing use higher number
	-a use -a mapping (off by default)
	-T number of threads
	-v version
	-h display usage message

example: bowtie_mapper.sh -d database.fasta -s COLI -1 ecoli_1.fastq -2 ecoli_2.fastq -a

EOF
}

#================================================================
# OPTION_PROCESSING
#================================================================
#Make sure the script is executed with arguments
if [ $? != 0 ] ; then
 usage >&2
 exit 1
fi


#DECLARE FLAGS AND VARIABLES
threads=1
offrate=1
cwd="$(pwd)"
a_mapping=""
group="NO_GROUP"
database="Database"
R1="R1"
R2="R2"

#PARSE VARIABLE ARGUMENTS WITH getops
#common example with letters, for long options check longopts2getopts.sh
options=":i:o:s:g:d:1:2:f:avh"
while getopts $options opt; do
	case $opt in
		i )
			input_dir=$OPTARG
			;;
		o )
			output_dir=$OPTARG
			;;
		s )
			sample=$OPTARG
			;;
		g)
			group=$OPTARG
			;;
		d )
			database=$OPTARG
			;;
		1 )
			R1=$OPTARG
			;;
		2 )
			R2=$OPTARG
			;;
		f )			
          	offrate=$OPTARG
      		;;
        T ) 
			threads=$OPTARG
            ;;
        a)
			a_mapping="-a"
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



#CHECK DEPENDENCIES FUNCTION
#This function check all dependencies listed and exits if any is missing
check_dependencies() {
	missing_dependencies=0
	for command in "$@"; do
		if ! [ -x "$(which $command 2> /dev/null)" ]; then
			echo "Error: Please install $command or make sure it is in your path" >&2
			let missing_dependencies++
		else
			echo "$command installed"
		fi
	done

	if [ $missing_dependencies -gt 0 ]; then 
		echo "ERROR: $missing_dependencies missing dependencies, aborting execution" >&2
		exit 1
	fi
}

#CHECK MANDATORY FIELDS
check_mandatory_files() {
	missing_files=0
	for file in "$@"; do
		if [ ! -f $file ]; then
			echo "$(basename $file)" "not supplied, please, introduce a valid file" >&2
			let missing_files++
		fi
	done

	if [ $missing_files -gt 0 ]; then 
		echo "ERROR: $missing_files missing files, aborting execution" >&2
		#exit 1
	fi
}


#================================================================
# MAIN_BODY
#================================================================
##CHECK DEPENDENCIES, MANDATORY FIELDS, FOLDERS AND ARGUMENTS

check_dependencies bowtie2-build bowtie2

check_mandatory_files $database $R1 $R2

if [ ! $sample ]; then
	echo "ERROR: please, provide a sample name"
	usage
	exit 1
fi

if [ ! $output_dir ]; then
	output_dir=$cwd"/$group/$sample/mapping/"
	echo "Default output directory is" $output_dir
	mkdir -p $output_dir
else
	echo "Output directory is" $output_dir
	mkdir -p $output_dir
fi



echo "Mapping $sample in $group Group "
	R1=$(find -L $trimmedDir/ -name "$sample"*"" -name "*_1_paired.fastq.gz" -type f | awk '/'"${sample}"'_/ && NR==1')
	R2=$(find -L $trimmedDir/ -name "$sample"*"" -name "*_2_paired.fastq.gz" -type f | awk '/'"${sample}"'_/ && NR==1')



 bash bowtie_mapper.sh -d ../../../REFERENCES/PLASMIDS/plasmid.all.genomic.dec212017.fasta_psi_90 \
 -g TEST \
 -s ABA622 \
 -1 ../../TRIMMED/ACIN_first/ABA622/ABA622_1_paired.fastq.gz \
 -2 ../../TRIMMED/ACIN_first/ABA622/ABA622_2_paired.fastq.gz

 #ABA622.sam

 bash sam_to_bam.sh -i GROUP/ABA622_2/mapping/ABA622_2.sam

 #ABA622.sorted.bam
 #ABA622.sorted.bam.bai

 bash get_coverage.sh -i TEST/ABA622/mapping/ABA622.sorted.bam -d ../../../REFERENCES/PLASMIDS/plasmid.all.genomic.dec212017.fasta_psi_90

 #ABA622.coverage

 bash adapt_filter_coverage.sh -i TEST/ABA622/mapping/ABA622.coverage -c 70


 #ABA622.coverage_adapted
 #ABA622.coverage_adapted_filtered_70  

bash filter_fasta.sh -i ../../../REFERENCES/PLASMIDS/plasmid.all.genomic.dec212017.fasta_psi_90 -f TEST/ABA622/mapping/ABA622.coverage_adapted_filtered_70

