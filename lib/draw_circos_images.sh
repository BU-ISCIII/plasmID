#!/bin/bash



group=$1
sample=$2

cdsDdbbFile="/processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PSMP_T/REFERENCES/PLASMIDS/plasmid.all.genomic.feb212017.bed"



mappedDir=$group/$sample/mapping
imageDir=$group/$sample/data

circos_conf_summary="config_files/circos_summary.conf"
circos_conf_individual="config_files/circos_individual.conf"
circosDir=$group/$sample/images

plasmidMapped=$mappedDir/$sample".coverage_adapted_clustered_ac"

karyotype_file_individual=$imageDir/$sample".karyotype_individual.txt"
karyotype_file_summary=$imageDir/$sample".karyotype_summary.txt"
abr_file=$imageDir/$sample".abr.coordinates"
replisome_file=$imageDir/$sample".inc.coordinates"
additional_file=$imageDir/$sample".additional.coordinates"
coverage_file=$imageDir/$sample".bedgraph_term"
cds_contig_file=$imageDir/$sample".gff.coordinates"

contig_file=$imageDir/$sample".plasmids.bed"
contig_file_complete=$imageDir/$sample".plasmids.complete"
links_file=$imageDir/$sample".plasmids.blast.links"


mkdir -p $circosDir

echo "Creating config file for circos in SAMPLE $sample FILE $circosDir/$sample.circos.conf"

awk '{gsub("PLASMID_KARYOTYPE","'$karyotype_file_summary'"); \
gsub("PLASMID_ANTIBIOTIC_RESISTANCE","'$abr_file'"); \
gsub("PLASMID_REPLISOME_PLASMIDFINDER","'$replisome_file'"); \
gsub("PLASMID_IS_ISFINDER","'$additional_file'"); \
gsub("PLASMID_COVERAGE_GRAPH","'$coverage_file'"); \
gsub("PLASMID_CDS_CONTIG","'$cds_contig_file'"); \
gsub("PLASMID_CDS_DDBB","'$cdsDdbbFile'"); \
gsub("PLASMID_CONTIGS","'$contig_file'"); \
gsub("PLASMID_LINKS","'$links_file'"); \
gsub("OUTPUTDIR","'$circosDir'"); \
gsub("IMAGENAME","'$imageName'"); \
print $0}' $circos_conf_summary > $circosDir/$sample".circos.conf"

echo "DONE Creating config file for circos in SAMPLE $sample"


echo "Creating config file for circos in SAMPLE $sample FILE $circosDir/$sample.circos.conf"

awk '{gsub("PLASMID_KARYOTYPE","'$karyotype_file_individual'"); \
gsub("PLASMID_ANTIBIOTIC_RESISTANCE","'$abrFile'"); \
gsub("PLASMID_REPLISOME_PLASMIDFINDER","'$replisome_file'"); \
gsub("PLASMID_IS_ISFINDER","'$additional_file'"); \
gsub("PLASMID_COVERAGE_GRAPH","'$coverage_file'"); \
gsub("PLASMID_CDS_CONTIG","'$cds_contig_file'"); \
gsub("PLASMID_CDS_DDBB","'$cdsDdbbFile'"); \
gsub("PLASMID_CONTIGS_COMPLETE","'$contig_file_complete'"); \
gsub("PLASMID_CONTIGS","'$contig_file'"); \
gsub("PLASMID_LINKS","'$linksFile'"); \
gsub("OUTPUTDIR","'$circosDir'"); \
print $0}' $circos_conf_individual > $circosDir/$sample"_individual.circos.conf"

echo "DONE Creating config file for circos in SAMPLE $sample"

echo "Executing circos in SAMPLE $sample FILE $circosDir/$sample.circos.conf"


#SAMPLE_SHOWN
#IMAGENAME_SAMPLE

for i in $(cat $plasmidMapped)
do
	awk '{gsub("SAMPLE_SHOWN","'$i'"); \
	gsub("IMAGENAME_SAMPLE_PLASMID","'$sample'_'$i'.png"); \
	print $0}' $circosDir/$sample"_individual.circos.conf" > $circosDir/$sample"_"$i"_individual.circos.conf"
	circos -conf $circosDir/$sample"_"$i"_individual.circos.conf" 
done 

circos -conf $circosDir/$sample".circos.conf"


echo "ALL DONE"
