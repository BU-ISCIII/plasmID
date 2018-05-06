#!/bin/bash

# Exit immediately if a pipeline, which may consist of a single simple command, a list, 
#or a compound command returns a non-zero status: If errors are not handled by user
set -e
#set -x

#=============================================================
# HEADER
#=============================================================

#INSTITUTION:ISCIII
#CENTRE:BU-ISCIII
#AUTHOR: Pedro J. Sola
VERSION=1.0 
#CREATED: 4 May 2018
#REVISION:
#DESCRIPTION:blast_to_bed script obtain a BED file with coordinates of local blast alignments matching some given conditions
#================================================================
# END_OF_HEADER
#================================================================

#SHORT USAGE RULES
#LONG USAGE FUNCTION
usage() {
	cat << EOF

blast_to_bed is a script than obtain a BED file with coordinates of local blast alignments matching some given conditions

usage : $0 <-i inputfile(.blast)> <-b coverage_file> [-o <directory>] [-c <int(0-100)>] [-s <suffix>] [-u] [-v] [-h]

	-i input file 
	-b blast identity cutoff (0 - 100), default 90
    -l blast length percentage cutoff (0 - 100), default 20, use 90 for genes
    -L blast length alignment cutoff, default 0, use 200 or 500 for contigs
	-o output directory (optional). By default the file is replaced in the same location
	-d database chraracter delimiter, default "_"
	-q database chraracter delimiter, default "_"
    -u unique. Outputs only one query entry per database entry
	-v version
	-h display usage message

example: blast_to_bed.sh -i ecoli_prefix.blast -b 80 -l 50 -q "|"

EOF
}

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
input_file="Input_file"
blast_id_cutoff=90
blast_len_percentage=20
blast_len_alignment=0
database_delimiter="_"
query_delimiter="_"
unique=false
suffix=""

#PARSE VARIABLE ARGUMENTS WITH getops
#common example with letters, for long options check longopts2getopts.sh
options=":i:b:q:d:o:l:L:uvh"
while getopts $options opt; do
	case $opt in
		i )
			input_file=$OPTARG
			;;
		b )
			if [ $OPTARG -lt 0 ] || [ $OPTARG -gt 100 ]; then
				echo "please, provide a percentage between 0 and 100"
				usage
				exit 1
			else
				blast_id_cutoff=$OPTARG
			fi
			;;
		o )
			output_dir=$OPTARG
			;;
		l )
			if [ $OPTARG -lt 0 ] || [ $OPTARG -gt 100 ]; then
				echo "please, provide a percentage between 0 and 100"
				usage
				exit 1
			else
				blast_len_percentage=$OPTARG
			fi
			;;
		L )
			blast_len_alignment=$OPTARG
			;;
        d )
			database_delimiter=$OPTARG
			;;
        q )
			query_delimiter=$OPTARG
			;;
        u )
			unique=true
            suffix=".unique.tmp"
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

bash check_mandatory_files.sh $input_file

blast_len_percentage_value=$(echo "($coverage_cutoff_input/100)" | bc -l)

echo "coverage" $coverage_cutoff

if [ ! $output_dir ]; then
	output_dir=$(dirname $input_file)
	#echo "Default output directory is" $output_dir
	mkdir -p $output_dir
else
	#echo "Output directory is" $output_dir
	mkdir -p $output_dir
fi


if [ ! $file_name ]; then
	file_name=$(basename $input_file | cut -d. -f1)
fi

echo "$(date)"
echo "Adapting blast to bed using" $(basename $input_file) "with:"
echo "Blast identity=" $blast_id_cutoff
echo "Min length aligned=" $blast_len_alignment
echo "Min len percentage=" $blast_len_percentage


cat $input_file |\
awk '
    {OFS="\t";
        split($1, database_name, "'"${database_delimiter}"'");
        split($2, query_name, ""'"${query_delimiter}"'"");
    };
    ($3 > '"${blast_id_cutoff}"')&&(($4/$14)>'"${blast_len_percentage}"')&&($4 >'"${blast_len_alignment}"')
    {print database_name[length(database_name)], $7, $8, query_name[1]}' \
> $output_dir/$file_name".bed"$suffix


if [ unique == "true"] then;
    awk '
        (!x[$2$3])
    ' \
    $output_dir/$file_name".bed"$suffix > $output_dir/$file_name".bed"

echo "$(date)"
echo "DONE adapting blast to bed
echo "File can be found at" $output_dir/$file_name".bed"
