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
#CREATED: 21 March 2018
#REVISION:
#		22 March 2018: Handle output directory by default the same as -f file
#DESCRIPTION:Script that extract sequences by term, either by key or file with a list
#AKNOWLEDGE: 
#		-Multiple arguments in one flag: https://stackoverflow.com/questions/7529856/retrieving-multiple-arguments-for-a-single-option-using-getopts-in-bash
#TODO:
#		-Add and remove sequences in the same execution
#================================================================
# END_OF_HEADER
#================================================================

#SHORT USAGE RULES
#LONG USAGE FUNCTION
usage() {
	cat << EOF

Filter_fasta script that extract sequences by term, either by key or file with a list

usage : $0 <-i file.fasta> <(-l term1 -l term2 -l term3 | -f file)> [-n <filename>] [-o <directory>] [N] [-v] [-h]

	-i fasta file to filter
	-o output directory (optional). By default the file is replaced in the same location
	-n file name (optional). By default is the same as -f file with .fasta extension
	-l list of key terms separated by space
	-N Use term to discard sequences with terms (Negative filter)
	-f file with a list of terms to filter
	-v version
	-h display usage message

example: filter_fasta.sh -i ecoli.fasta -l NC00012 -l WC52247 -l hypothetical -l partial -n NAME
		filter_fasta.sh -i ecoli.fasta -l "NC00012 WC52247 hypothetical partial"
		filter_fasta.sh -i ecoli.fasta -f list_with_terms.txt

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
term_option=false
file_option=false
negative_filter=""
cwd="$(pwd)"
input_file="Input_file"

#PARSE VARIABLE ARGUMENTS WITH getops
#common example with letters, for long options check longopts2getopts.sh
options=":i:o:n:l:f:Nvh"
while getopts $options opt; do
	case $opt in
		i )
			input_file=$OPTARG
			;;
		o )
			output_dir=$OPTARG
			;;
		n )
			file_name=$OPTARG
			;;
		N )
			negative_filter="!"
			;;
		l )
			terms_for_filtering+=($OPTARG)
			term_option=true
			;;
		f )
			file_for_filtering=$OPTARG
			file_option=true
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


if [ $file_option = true ]; then
	output_dir=$(dirname $file_for_filtering)
	file_name=$(basename $file_for_filtering | cut -d. -f1)
	mkdir -p $output_dir
else 

#$file_for_filtering

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



if [ $term_option = true ] && [ $file_option = false ]; then

	list_terms_listed=$(for term in "${terms_for_filtering[@]}"; do echo "$term"; done) #process terms into list
	final_list_terms_regexp=$(printf "%s|" $list_terms_listed | sed 's/|$//g') #suitable for regexp

elif [ $term_option = false ] && [ $file_option = true ]; then

	bash check_mandatory_files.sh $file_for_filtering

	final_list_terms_regexp=$(printf "%s|" $(cat $file_for_filtering) | sed 's/|$//g')
else

	bash check_mandatory_files.sh $file_for_filtering

	list_terms_listed=$(for term in "${terms_for_filtering[@]}"; do echo "$term"; done)
	list_terms_regexp_term=$(printf "%s|" $list_terms_listed | sed 's/|$//g')
	list_terms_regexp_file=$(printf "%s|" $(cat $file_for_filtering) | sed 's/|$//g')
	final_list_terms_regexp=$(echo $list_terms_regexp_term"|"$list_terms_regexp_file) #concat all regexp into one
fi


echo "$(date)"
echo "Filtering terms on file" $(basename $input_file)
seq_number_prev=$(cat $input_file | grep ">" | wc -l)

awk '
	BEGIN {RS=">"} 
	'"${negative_filter}"'/'"${final_list_terms_regexp}"'/ {print ">"$0}
	' $input_file > $output_dir/$file_name"_filtered.fasta"

echo "$(date)"
echo "DONE Filtering terms on file" $(basename $input_file)
seq_number_post=$(cat $output_dir/$file_name_filtered.fasta | grep ">" | wc -l)
echo "File with filtered sequences can be found in" $output_dir/$file_name"_filtered.fasta"

echo "Previous number of sequences=" $seq_number_prev
echo "Post number of sequences=" $seq_number_post

echo $file_name