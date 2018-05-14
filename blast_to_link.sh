#!/bin/bash

# Exit immediately if a pipeline, which may consist of a single simple command, a list, 
#or a compound command returns a non-zero status: If errors are not handled by user
#set -e
#set -x

#=============================================================
# HEADER
#=============================================================

#INSTITUTION:ISCIII
#CENTRE:BU-ISCIII
#AUTHOR: Pedro J. Sola
VERSION=1.0 
#CREATED: 14 May 2018
#
#DESCRIPTION:blast_to_link script to obtain a link file that represent duplications between all members of the query 
#================================================================
# END_OF_HEADER
#================================================================

#SHORT USAGE RULES
#LONG USAGE FUNCTION
usage() {
	cat << EOF

blast_to_bed is a script than obtain a BED file with coordinates of local blast alignments matching some given conditions

usage : $0 <-i inputfile(.blast)> <-b id cutoff> [-o <directory>] [-b <int(0-100)>] [-l <int(0-100)>] [-L <int>]
		[-p <prefix>] [-d <delimiter>] [-D (l|r)] [-q <delimiter>] [-Q (l|r)] [-I] [-u] [-v] [-h]

	-i input file 
	-b blast identity cutoff (0 - 100), default 90
	-l blast length percentage cutoff (0 - 100), default 20, use 90 for genes
	-o output directory (optional). By default the file is replaced in the same location
	-q database chraracter delimiter, default "_"
	-Q query field to retrieve (l=left, r=right), default left
	-d database chraracter delimiter, default "_"
	-D database field to retrieve (l=left, r=right), default right
	-I contig mode
	-u unique. Outputs only one query entry per database entry
	-v version
	-h display usage message

example: blast_to_link.sh -i ecoli_prefix.blast -b 80 -l 50
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
blast_len_percentage=10
blast_len_alignment=0
database_delimiter="_"
database_field=r
query_delimiter="_"
query_field=l
unique=false
suffix=""
id_circos=false
id_output=""

#PARSE VARIABLE ARGUMENTS WITH getops
#common example with letters, for long options check longopts2getopts.sh
options=":i:b:q:Q:d:D:o:l:L:Iuvh"
while getopts $options opt; do
	case $opt in
		i )
			input_file=$OPTARG
			;;
		b )
			if [ $OPTARG -lt 0 ] || [ $OPTARG -gt 100 ]; then
				echo "please, provide a percentage between 0 and 100"
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
				exit 1
			else
				blast_len_percentage=$OPTARG
			fi
			;;
		d )
			database_delimiter=$OPTARG
			;;
		D )
			database_field=$OPTARG
			;;
        q )
			query_delimiter=$OPTARG
			;;
		Q )
			query_field=$OPTARG
			;;
        u )
			unique=true
            suffix=".unique.tmp"
			;;
        I)
			id_circos=true
            id_output=",\"id=\"database_name[length(database_name)]"
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


blast_len_percentage_value=$(echo "($blast_len_percentage/100)" | bc -l)
#blast_len_percentage_decimal=$(echo $blast_len_percentage_value | sed 's/0\{1,\}$//')


if [ ! $output_dir ]; then
	output_dir=$(dirname $input_file)
	#echo "Default output directory is" $output_dir
	mkdir -p $output_dir
else
	#echo "Output directory is" $output_dir
	mkdir -p $output_dir
fi


if [ ! $file_name ]; then
	file_name=$(basename $input_file | cut -d. -f1,2)
fi

##CHECK FIELDS TO RETRIEVE

if [ "$database_field" == "l" ] || [ "$database_field" == "r" ]; then

	if [ $database_field == l ]; then
		database_field="1"
	else
		database_field="length(database_name)"
	fi
	
else
	echo "Please introduce r or l for database"
	exit 1
fi

if [ $query_field == "l" ] || [ $query_field == "r" ]; then

	if [ $query_field == l ]; then
		query_field="1"
	else
		query_field="length(query_name)"
	fi
	
else

	echo "Please introduce 0 or 1 for query"
	exit 1
fi

echo "$(date)"
echo "Adapting blast to bed using" $(basename $input_file) "with:"
echo "Blast identity=" $blast_id_cutoff
echo "Min len percentage=" $blast_len_percentage


#cat $input_file |\
#awk '
#	{OFS="\t"
#	split($2, database_name, "'"${database_delimiter}"'")
#	split($1, query_name, "'"${query_delimiter}"'")}
#	(($3 > '"${blast_id_cutoff}"')&&(($4/$14) > '"${blast_len_percentage_value}"')&&($4 > '"${blast_len_alignment}"')) \
#	{print query_name['"$query_field"'], $7, $8, database_name['"$database_field"']'"$id_output"'}
#	' \
#> $output_dir/$file_name".bed"$suffix

awk '(($4/$14)>'"${blast_len_percentage_value}"') && !contigPlasmid[$2$1]++ {print $2$1}' $imageDir/$sample".plasmids.blast" > $output_dir/$file_name".dictionary_link"

awk 'NR==FNR{contigPlasmid[$2]=$2;next} \
{split($2, contigname, "_"); header=$2$1}{if ((header in contigPlasmid) && ($3>90) && (($4/$14)>0.05)) print contigname[length(contigname)], $9,$10,$1,$7,$8, "id=contig_"contigname[length(contigname)]}' \
$imageDir/$sample".dictionary_covered.txt" $imageDir/$sample".plasmids.blast" > $imageDir/$sample".plasmids.blast.links"

#sample.blast.links (file with matching coordinates between coordinates in contigs and plasmids)
#contig_1    915668  915969  NG_048025.1     1178    1476    id=contig_1
#contig_1    123340  123574  NG_048025.1     1361    1127    id=contig_1
#contig_2    1       121     NZ_CP010574.1   66599   66479   id=contig_2
#contig_2    1       121     NZ_CP008930.1   122307  122187  id=contig_2
#contig_2    1       121     NZ_CP011577.1   112628  112748  id=contig_2

##Change coordinates from contig --> plasmid to plasmid-->plasmid in order to represent them in CIRCOSS

awk 'BEGIN{OFS="\t";}{if($1 != savedNode){savedNode= $1; delete chr} else{for(i in chr){print $4" "$5" "$6" "chr[i]" id="savedNode}}chr[$4$5$6] = $4" "$5" "$6}' \
$imageDir/$sample".plasmids.blast.links" > $imageDir/$sample".plasmids.blast.links.sorted"


if [ "$unique" == "true" ]; then
    awk '
        (!x[$1$4]++)
    ' $output_dir/$file_name".bed"$suffix \
> $output_dir/$file_name".bed"
fi

echo "$(date)"
echo "DONE adapting blast to bed"
echo "File can be found at" $output_dir/$file_name".bed"