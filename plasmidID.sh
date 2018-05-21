#!/bin/bash

# Exit immediately if a pipeline, which may consist of a single simple command, a list, 
#or a compound command returns a non-zero status: If errors are not handled by user
set -e
# Treat unset variables and parameters other than the special parameters ‘@’ or ‘*’ as an error when performing parameter expansion.

#Print everything as if it were executed, after substitution and expansion is applied: Debug|log option
#set -x

#=============================================================
# HEADER
#=============================================================

#INSTITUTION:ISCIII
#CENTRE:BU-ISCIII
#AUTHOR: Pedro J. Sola
VERSION=Beta 
#CREATED: 15 March 2018
#
#ACKNOLEDGE: longops2getops.sh: https://gist.github.com/adamhotep/895cebf290e95e613c006afbffef09d7
#
#DESCRIPTION: plasmidID is a computational pipeline tha reconstruct and annotate the most likely plasmids present in one sample		
#
#
#================================================================
# END_OF_HEADER
#================================================================

#SHORT USAGE RULES
#LONG USAGE FUNCTION
usage() {
	cat << EOF

plasmidID is a computational pipeline tha reconstruct and annotate the most likely plasmids present in one sample

usage : $0 [-i <inputfile>] [-o <directory>] <-d database(fasta)> <-s sample_name> <-1 R1> <-2 R2> 
		[-g group_name] [-f <int>] [-T <int>] [-a] [-v] [-h]

	-1 reads corresponding to paired-end R1 (mandatory)
	-2 reads corresponding to paired-end R2 (mandatory)
	-d database to map and reconstruct (mandatory)
	-s sample name (mandatory)
	-g group name (optional). If unset, samples will be gathered in NO_GROUP group
	


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



# translate long options to short
reset=true
for arg in "$@"
do
    if [ -n "$reset" ]; then
      unset reset
      set --      # this resets the "$@" array so we can rebuild it
    fi
    case "$arg" in
##MANDATORY MINIMAL OPTIONS
    	--R1)		set -- "$@"	-1 ;;
		--R2)		set -- "$@"	-2 ;;
		--database)		set -- "$@"	-d ;;
		--sample)		set -- "$@"	-s ;;
		--group)		set -- "$@"	-g ;;
##
		--R1)		set -- "$@"	-1 ;;
       	--help)    	set -- "$@" -h ;;
       	--version) 	set -- "$@" -v ;;
       	--config)  	set -- "$@" -c ;;
       # pass through anything else
       *)         set -- "$@" "$arg" ;;
    esac
done

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
		1 )
			r1_file=$OPTARG
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

#check_dependencies bowtie2-build bowtie2

#check_mandatory_files $database $R1 $R2

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





####TRIMMING############
########################
./quality_trim.sh -1 ../../TRIMMED/ACIN_first/ABA622/ABA622_1_paired.fastq.gz -2 ../../TRIMMED/ACIN_first/ABA622/ABA622_2_paired.fastq.gz -s ABA622 -g TEST -T 8

#TEST/ABA622/trimmed/ABA622_1_paired.fastq.gz
#TEST/ABA622/trimmed/ABA622_1_unpaired.fastq.gz
#TEST/ABA622/trimmed/ABA622_2_paired.fastq.gz
#TEST/ABA622/trimmed/ABA622_2_unpaired.fastq.gz

./spades_assembly.sh -q PS_second/PAE1286/trimmed/ -T 8




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

 bash adapt_filter_coverage.sh -i TEST/ABA622/mapping/ABA622.coverage -c 50


 #ABA622.coverage_adapted
 #ABA622.coverage_adapted_filtered_50  

bash filter_fasta.sh -i ../../../REFERENCES/PLASMIDS/plasmid.all.genomic.dec212017.fasta_psi_90 -f TEST/ABA622/mapping/ABA622.coverage_adapted_filtered_50

#ABA622.coverage_adapted_filtered_50_term.fasta

bash cdhit_cluster.sh -i TEST/ABA622/mapping/ABA622.coverage_adapted_filtered_50_term.fasta -c 0.8 -M 5000 -T 0

#ABA622.coverage_adapted_filtered_50_term.fasta_80 >> ARCHIVO FINAL DE PLÁSMIDOS FILTRADO Y CLUSTERIZADO
########################################################################################################
#ABA622.coverage_adapted_filtered_50_term.fasta_80.clstr

bash process_cluster_output.sh -i TEST/ABA622/mapping/ABA622.coverage_adapted_filtered_50_term.fasta_80 -b TEST/ABA622/mapping/ABA622.coverage_adapted -c 50

#ABA622.coverage_adapted_clustered
#ABA622.coverage_adapted_clustered_percentage
#ABA622.coverage_adapted_clustered_ac >> PROBABLY NOT NEEDED
############################################################

###############################################################################################################################################
######################## DONE WITH MAPPING AND CLUSTERING #####################################################################################
###############################################################################################################################################

bash build_karyotype.sh -i TEST/ABA622/mapping/ABA622.coverage_adapted_clustered -K 50 -k 50 -o TEST/ABA622/data

#group/sample/data
#ABA622.karyotype_individual.txt
#ABA622.karyotype_individual.txt

./get_coverage.sh -i TEST/ABA622/mapping/ABA622.sorted.bam -p -o TEST/ABA622/data/

#ABA622.bedgraph

./filter_fasta.sh -i TEST/ABA622/data/ABA622.bedgraph -f TEST/ABA622/mapping/ABA622.coverage_adapted_clustered_ac -G

#ABA622.bedgraph_term

contigFile=$(find -L $contigDir/ -name "scaffolds.fasta" -type f 2> /dev/null| awk '/'"${sample}"'/' | awk 'NR==1')

/prokka_annotation.sh -i ../../ASSEMBLY/ACIN_first/ABA622/scaffolds.fasta -p ABA622 -O TEST/ABA622/data/ -c

#ABA622.fna
#ABA622.gff


./blast_align.sh -i TEST/ABA622/data/ABA622.fna -d TEST/ABA622/mapping/ABA622.coverage_adapted_clustered_ac_term.fasta -o TEST/ABA622/data/ -p plasmids2

#ABA622.plasmids.blast


./blast_to_bed.sh -i TEST/ABA622/data/ABA622.plasmids.blast -l 0 -L 500 -d - -q _ -Q r -I

#ABA622.plasmids.bed

./blast_to_complete.sh -i TEST/ABA622/data/ABA622.plasmids.blast

#ABA622.plasmids.complete

./blast_to_link.sh -i TEST/ABA622/data/ABA622.plasmids.blast -I

#ABA622.plasmids.links
#ABA622.plasmids.blast.links

./gff_to_bed.sh -i TEST/ABA622/data/ABA622.gff -L -u

#ABA622.gff.bed

./coordinate_adapter.sh -i TEST/ABA622/data/ABA622.gff.bed -l TEST/ABA622/data/ABA622.plasmids.blast.links -n 10000

#ABA622.gff.coordinates




######################### ABR _ INCLUDE FILENAME

./blast_align.sh -i /processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PSMP_T/ANALYSIS/PLASMIDID/references/ARGannot.r1.pID.fasta -d TEST/ABA622/data/ABA622.fna -o TEST/ABA622/data -p abr -f ABA622

#ABA622.abr.blast

./blast_to_bed.sh -i TEST/ABA622/data/ABA622.abr.blast -b 95 -l 90 -d _ -D r -q " " -Q r

#ABA622.abr.bed

./coordinate_adapter.sh -i TEST/ABA622/data/ABA622.abr.bed -l TEST/ABA622/data/ABA622.plasmids.blast.links -u

#ABA622.abr.coordinates

###################### INC

./blast_align.sh -i /processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PSMP_T/ANALYSIS/PLASMIDID/references/plasmidFinder_01_26_2018.fsa -d TEST/ABA622/data/ABA622.fna -o TEST/ABA622/data -p inc -f ABA622

#ABA622.inc.blast

/blast_to_bed.sh -i TEST/ABA622/data/ABA622.inc.blast -b 95 -l 80 -d _ -D r -q _ -Q l
#ABA622.inc.bed

./coordinate_adapter.sh -i TEST/ABA622/data/ABA622.inc.bed -l TEST/ABA622/data/ABA622.plasmids.blast.links -u

#ABA622.inc.coordinates