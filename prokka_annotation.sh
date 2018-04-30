#!/bin/bash

set -e

#=============================================================
# HEADER
#=============================================================

#INSTITUTION:ISCIII
#CENTRE:BU-ISCIII
#AUTHOR: Pedro J. Sola
VERSION=1.0 
#CREATED: 30 May 2018
#REVISION:
#DESCRIPTION:Script that uses prokka to annotate a FASTA file
#
#DOCUMENTATION
#
#Prokka outputs the fasta headers as:
# gnl|center|locustag_01
# gnl|center|locustag_02
#
#================================================================
# END_OF_HEADER
#================================================================

#SHORT USAGE RULES
#LONG USAGE FUNCTION
usage() {
	cat << EOF

Prokka_annotation is a script that uses prokka to annotate a FASTA file

usage : $0 <-i inputfile(FASTA)> <-p prefix> [-o <directory>] [-n <filename>] [-k <kingdom>]
		[-T <threads>] [-g group_name][-G genus] [-S species] [-v] [-h]

	-i input file in FASTA format
	-o output directory
	-p prefix for sample identification (mandatory) and outpùt file name
	-n file name
	-k kingdom (Bacteria by default)
	-g group name (optional). If unset, samples will be gathered in NO_GROUP group
	-G sample genus in case is known by user
	-S sample species in case is known by user
	-T number of threads
	-v version
	-h display usage message


Output directory is the same as input directory by default

example: prokka_annotation -i ecoli.fasta -p ECO -T 5
		 

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
group="NO_GROUP"
input_file="Input_file"
threads=1

#PARSE VARIABLE ARGUMENTS WITH getops
#common example with letters, for long options check longopts2getopts.sh
options=":i:o:p:n:k:g:G:S:T:vh"
while getopts $options opt; do
	case $opt in
		i )
			input_file=$OPTARG
			;;
		
		o )
			output_dir=$OPTARG
			;;
		p)
			prefix=$OPTARG
			file_name=$OPTARG
			;;
		M )
			max_memory=$OPTARG
			;;
		
		k )			
          	kingdom=$OPTARG
          	;;
        g )			
          	group=$OPTARG
          	;;
        S )			
          	species=$OPTARG
          	;;
        G)			
          	genus=$OPTARG
          	;;
        T)			
          	threads=$OPTARG
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

bash check_dependencies.sh prokka



if [ ! $output_dir ]; then
	output_dir=$(dirname $input_file)
	echo "Default output directory is" $output_dir
	mkdir -p $output_dir
else
	echo "Output directory is" $output_dir
	mkdir -p $output_dir
fi

if [ ! $file_name ]; then
	file_name=$(basename $input_file)
	echo "filename is" $file_name
fi


##CD-HIT EXECUTION

echo "$(date)"
echo "Clustering sequences with identity" $cluster_percentage"% or higher"
echo "Using" $cd_hit_command "with file" $input_file
seq_number_prev_clstr=$(cat $input_file | grep ">" | wc -l)

cd $(dirname $input_file)

if [ -f $output_dir/$file_name""_""$cluster_percentage ]; then \
	echo "Found a clustered file for sample" $file_name;
	echo "Omitting clustering process calculation"
	exit 1
else
	if [ $cd_hit_command  == "psi-cd-hit.pl" ]; then 
		$cd_hit_command -i $(basename $input_file) -o $file_name""_""$cluster_percentage -c $cluster_cutoff -G $global_psi_cd_hit -g 1 -prog $psi_cd_hit_program -circle $is_circle
	
	else

		$cd_hit_command -i $(basename $input_file) -o $file_name""_""$cluster_percentage -c $cluster_cutoff -n $word_size -d 0 -s $length_cutoff -B 1 -M $max_memory -T $threads
	fi
fi

seq_number_post_clstr=$(cat $$file_name""_""$cluster_percentage | grep ">" | wc -l)

echo "$(date)"
echo "done Clustering sequences with identity" $cluster_percentage"% or higher"
echo "fasta file can be found in" $output_dir/$file_name""_""$cluster_percentage
echo "Previous number of sequences=" $seq_number_prev_clstr
echo "Number of sequences after clustering=" $seq_number_post_clstr
cd $cwd




#contigFile=$(find $contigDir/ -name "contigs.fasta" -type f | awk '/'"${sample}"'/ ')
contigFile=$(find -L $contigDir/ -name "scaffolds.fasta" -type f 2> /dev/null| awk '/'"${sample}"'/' | awk 'NR==1')
#2> /dev/null avoid permission denied error pront for non contig directories
#contigFile=$nadalesDir/A_257_index.final.scaffolds.fasta
if [ -f $imageDir/$sample".fna" -a -f $imageDir/$sample".gff"  ];then \

	echo "Found an annotation file for sample" $sample;
	echo "Omitting annotation with prokka"

elif [ ! -f $contigFile ]; then \

	echo "ERROR: File with contigs was not found for sample" $sample;
	echo "Exit"
	exit

else

	echo -e "contig file found at:" $contigFile
	echo -e "##### Anotating contigs with Prokka \n"

	prokka --force --outdir $imageDir \
	--prefix $sample \
	--addgenes \
	--kingdom Bacteria \
	--usegenus \
	--centre CNM \
	--locustag $sample \
	--compliant \
	--cpus 5 \
	$contigFile

	echo "##### DONE Anotating contigs with Prokka"

fi

#--genus Klebsiella \
#--species pneumoniae \
