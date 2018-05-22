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

usage : $0 [-1 <R1>] [-2 <R2>] <-d database(fasta)> <-s sample_name> [-g group_name] [options]

	Basic options:
	--default
	--explore

	Mandatory input data:
	-1 | --R1	<filename>	reads corresponding to paired-end R1 (mandatory)
	-2 | --R2	<filename>	reads corresponding to paired-end R2 (mandatory)
	-d | --database	<filename>	database to map and reconstruct (mandatory)
	-s | --sample	<string>	sample name (mandatory)

	Optional input data:
	-g | --group	<string>	group name (optional). If unset, samples will be gathered in NO_GROUP group
	-c | --contig	<filename>	file with contigs. If supplied, plasmidID will not assembly reads
	-a | --annotate <filename>	file with sequences to draw in the final images
	-o 		<output_dir>	output directory, by default is the current directory

	Pipeline options:
	-C | --coverage-cutoff	<int>	minimun coverage percentage to select a plasmid as scafold (0-100), default 80
	-S | --coverage-summary	<int>	minimun coverage percentage to include plasmids in summary image (0-100), default 90
	-f | --cluster		<int>	identity percentage to cluster plasmids into the same representative sequence (0-100), default 80
	-i | --alignment-identity <int>	minimun identity percentage aligned for a contig to annotate
	-l | --alignment-percentage <int>	minimun length percentage aligned for a contig to annotate
	-L | --length-total	<int>	minimun alignment length to filter blast analysis

	--only-reconstruct	Database supplied will not be filtered and all sequences will be used as scaffold 
	
	Additional options:

	-T | --threads	<int>	number of threads
	-v | --version		version
	-h | --help		display usage message

example: ./plasmidID.sh -1 ecoli_R1.fastq.gz -2 ecoli_R2.fastq.gz -d database.fasta -s ECO_553 -G ENTERO

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
		--database)	set -- "$@"	-d ;;
		--sample)	set -- "$@"	-s ;;
##OPTIONAL
		--group)	set -- "$@"	-g ;;
		--contig)	set -- "$@"	-c ;;
		--annotate)	set -- "$@"	-a ;;
##PIPELINE
		--coverage-cutoff)	set -- "$@"	-C ;;
		--coverage-summary)	set -- "$@"	-S ;;
		--cluster)	set -- "$@"	-f ;;
		--alignment-percentage)	set -- "$@"	-l ;;
		--length-total)	set -- "$@"	-L ;;
##ADDITIONAL
		--threads)		set -- "$@"	-T ;;
       	--help)    	set -- "$@" -h ;;
       	--version) 	set -- "$@" -v ;;
       # pass through anything else
       *)         set -- "$@" "$arg" ;;
    esac
done

#DECLARE FLAGS AND VARIABLES
threads=1
cwd="$(pwd)"
group="NO_GROUP"
database="Database"
coverage_cutoff=80
coverage_summary=90
cluster_cutoff=80
alignment_percentage=50
R1="R1_file"
R2="R2_file"

#PARSE VARIABLE ARGUMENTS WITH getops
#common example with letters, for long options check longopts2getopts.sh
options=":1:2:d:s:g:c:a:o:C:S:f:l:L:T:vh"
while getopts $options opt; do
	case $opt in
		1 )
			r1_file=$OPTARG
			;;
		2 )
			r2_file=$OPTARG
			;;
		d )
			database=$OPTARG
			;;
		s )
			sample=$OPTARG
			;;
		g)
			group=$OPTARG
			;;
		c)
			contigs=$OPTARG
			;;
		g)
			group=$OPTARG
			;;
		a )
			annotation+=($OPTARG)
			;;
		C )
			coverage_cutoff=$OPTARG
			;;
		S )
			coverage_summary=$OPTARG
			;;
		f )			
          	cluster_cutoff=$OPTARG
      		;;
      	l )			
          	alignment_percentage=$OPTARG
      		;;
      	L )			
          	alignment_total=$OPTARG
      		;;
        T ) 
			threads=$OPTARG
            ;;
        o)
			output_dir=$OPTARG
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

.lib/check_mandatory_files $r1_file $r2_file $database 

if [ ! $sample ]; then
	echo "ERROR: please, provide a sample name"
	usage
	exit 1
fi

if [ ! $output_dir ]; then
	output_dir=$cwd
	echo "Default output directory is" $output_dir
	mkdir -p $output_dir
else
	echo "Output directory is" $output_dir
	mkdir -p $output_dir
fi


####TRIMMING#################################################################
#############################################################################

.lib/quality_trim.sh -1 $r1_file -2 $r2_file -s $sample -g $group -T $threads

#group/sample/trimmed/sample_1_paired.fastq.gz
#group/sample/trimmed/sample_1_unpaired.fastq.gz
#group/sample/trimmed/sample_2_paired.fastq.gz
#group/sample/trimmed/sample_2_unpaired.fastq.gz


####ASSEMLY##################################################################
#############################################################################

.lib/spades_assembly.sh -q $group/$sample/trimmed/ -T $threads

#group/sample/assembly/scaffolds.fasta

####MAPPING##################################################################
#############################################################################

.lib/bowtie_mapper.sh -d $database \
 -g $group \
 -s $sample \
 -1 $group/$sample/trimmed/$sample"_1_paired.fastq.gz" \
 -2 $group/$sample/trimmed/$sample"_2_paired.fastq.gz"

 #group/sample/mapping/sample.sam

 .lib/sam_to_bam.sh -i $group/$sample/mapping/$sample.sam

 #group/sample/mapping/sample.bam
 #group/sample/mapping/sample.bam.bai


####COVERAGE FILTERING#######################################################
#############################################################################


 .lib/get_coverage.sh -i $group/$sample/mapping/$sample".sorted.bam" -d $database

 #group/sample/mapping/sample.coverage

 .lib/adapt_filter_coverage.sh -i $group/$sample/mapping/$sample".coverage" -c $coverage_cutoff

 #sample.coverage_adapted
 #sample.coverage_adapted_filtered_80










 

bash filter_fasta.sh -i ../../../REFERENCES/PLASMIDS/plasmid.all.genomic.dec212017.fasta_psi_90 -f group/sample/mapping/sample.coverage_adapted_filtered_50

#sample.coverage_adapted_filtered_50_term.fasta

bash cdhit_cluster.sh -i group/sample/mapping/sample.coverage_adapted_filtered_50_term.fasta -c 0.8 -M 5000 -T 0

#sample.coverage_adapted_filtered_50_term.fasta_80 >> ARCHIVO FINAL DE PLÁSMIDOS FILTRADO Y CLUSTERIZADO
########################################################################################################
#sample.coverage_adapted_filtered_50_term.fasta_80.clstr

bash process_cluster_output.sh -i group/sample/mapping/sample.coverage_adapted_filtered_50_term.fasta_80 -b group/sample/mapping/sample.coverage_adapted -c 50

#sample.coverage_adapted_clustered
#sample.coverage_adapted_clustered_percentage
#sample.coverage_adapted_clustered_ac >> PROBABLY NOT NEEDED
############################################################

###############################################################################################################################################
######################## DONE WITH MAPPING AND CLUSTERING #####################################################################################
###############################################################################################################################################

bash build_karyotype.sh -i group/sample/mapping/sample.coverage_adapted_clustered -K 50 -k 50 -o group/sample/data

#group/sample/data
#sample.karyotype_individual.txt
#sample.karyotype_individual.txt

./get_coverage.sh -i group/sample/mapping/sample.sorted.bam -p -o group/sample/data/

#sample.bedgraph

./filter_fasta.sh -i group/sample/data/sample.bedgraph -f group/sample/mapping/sample.coverage_adapted_clustered_ac -G

#sample.bedgraph_term

contigFile=$(find -L $contigDir/ -name "scaffolds.fasta" -type f 2> /dev/null| awk '/'"${sample}"'/' | awk 'NR==1')

/prokka_annotation.sh -i ../../ASSEMBLY/ACIN_first/sample/scaffolds.fasta -p sample -O group/sample/data/ -c

#sample.fna
#sample.gff


./blast_align.sh -i group/sample/data/sample.fna -d group/sample/mapping/sample.coverage_adapted_clustered_ac_term.fasta -o group/sample/data/ -p plasmids2

#sample.plasmids.blast


./blast_to_bed.sh -i group/sample/data/sample.plasmids.blast -l 0 -L 500 -d - -q _ -Q r -I

#sample.plasmids.bed

./blast_to_complete.sh -i group/sample/data/sample.plasmids.blast

#sample.plasmids.complete

./blast_to_link.sh -i group/sample/data/sample.plasmids.blast -I

#sample.plasmids.links
#sample.plasmids.blast.links

./gff_to_bed.sh -i group/sample/data/sample.gff -L -u

#sample.gff.bed

./coordinate_adapter.sh -i group/sample/data/sample.gff.bed -l group/sample/data/sample.plasmids.blast.links -n 10000

#sample.gff.coordinates




######################### ABR _ INCLUDE FILENAME

./blast_align.sh -i /processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PSMP_T/ANALYSIS/PLASMIDID/references/ARGannot.r1.pID.fasta -d group/sample/data/sample.fna -o group/sample/data -p abr -f sample

#sample.abr.blast

./blast_to_bed.sh -i group/sample/data/sample.abr.blast -b 95 -l 90 -d _ -D r -q " " -Q r

#sample.abr.bed

./coordinate_adapter.sh -i group/sample/data/sample.abr.bed -l group/sample/data/sample.plasmids.blast.links -u

#sample.abr.coordinates

###################### INC

./blast_align.sh -i /processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PSMP_T/ANALYSIS/PLASMIDID/references/plasmidFinder_01_26_2018.fsa -d group/sample/data/sample.fna -o group/sample/data -p inc -f sample

#sample.inc.blast

/blast_to_bed.sh -i group/sample/data/sample.inc.blast -b 95 -l 80 -d _ -D r -q _ -Q l
#sample.inc.bed

./coordinate_adapter.sh -i group/sample/data/sample.inc.bed -l group/sample/data/sample.plasmids.blast.links -u

#sample.inc.coordinates