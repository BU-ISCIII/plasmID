#!/bin/bash

<<Usage
#Copy this command to generate a file (commands_placnet_by_sample.txt) with the order to execute this script for each sample


group=references_simulation
samplesid=/processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PS_MP_T/DOC/SAMPLES/samplesid_$group.txt

cat $samplesid | while read -r line
do echo "qsub -V -b y -j y -l h_vmem=20G -pe openmp 10 -cwd -q all.q -N DRW_${line}_PLSMD bash 3_Draw_plasmid_with_more_coverage.sh $line $group"
done > 3_commands_plasmid_drawing_$group.txt

bash 3_commands_plasmid_drawing_$group.txt
Usage

###############1 Map reads on plasmid database & filter positive matches
###############2 Find plasmid with more coverage
###############3 Draw plasmid with more coverage (>90%)

<<outfmt6
#1	 	Query label.(qseqid)
#2	 	Target or subject(database sequence or cluster centroid) label. (sseqid)
#3	 	Percent identity. (pident)
#4	 	Alignment length. (length)
#5	 	Number of mismatches. (mismatch)
#6	 	Number of gap opens. (gapopen)
#7	 	Start position in query. Query coordinates start with 1 at the first base in the sequence as it appears in the input file. For translated searches (nucleotide queries, protein targets), query start<end for +ve frame and start>end for -ve frame. (qstart)
#8	 	End position in query. (qend)
#9	 	Start position in target. Target coordinates start with 1 at the first base in sequence as it appears in the database. For untranslated nucleotide searches, target start<end for plus strand, start>end for a reverse-complement alignment. (sstart)
#10	 	End position in target. (send)
#11	 	E-value calculated using Karlin-Altschul statistics. (evalue)
#12	 	Bit score calculated using Karlin-Altschul statistics. (bitscore)
#13		Lenght of query (qlen)
#14		Length of target (slen)

outfmt6

#IN: List of ac mapped more than 90% in plasmid DDBB clustered with CCD-HIT (90%)


sample=$1
group=$2
plasmidDdbb="/processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PS_MP_T/REFERENCES/PLASMIDS/plasmid.all.genomic.feb212017.fna"
abReference="/processing_Data/bioinformatics/references/resistance/ARGANNOT/20170213/ARGannot.r1.fasta"
plasmidFinderReference="/processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PS_MP_T/REFERENCES/PLASMIDS/PLASMIDFINDER/plasmidFinder_02_10_2017.fsa"
isFinderReference="/processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PS_MP_T/REFERENCES/PLASMIDS/ISFINDER/ISFinder_2.fasta"
trimmedDir="/processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PS_MP_T/ANALYSIS/TRIMMED"
mappedDir="/processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PS_MP_T/ANALYSIS/MAPPING/PLASMIDS"
imageDir="/processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PS_MP_T/ANALYSIS/IMAGES/PLASMIDS"
contigDir="/processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PS_MP_T/ANALYSIS/ASSEMBLY"
#plasmidLenght="/processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PS_MP_T/REFERENCES/PLASMIDS/NCBI_PLASMID_LENGTH.tsv"
plasmidLenght="/processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PS_MP_T/REFERENCES/PLASMIDS/plasmid.all.genomic.feb212017.length"
plasmidMapped=$mappedDir/$group/$sample/$sample".ac.covered.txt"
plasmidMappedFasta=$mappedDir/$group/$sample/$sample".ac.covered.fasta"
sampleCoverage=$mappedDir/$group/$sample/$sample".coverage.txt"
karyotypeFile=$imageDir/$group/$sample/$sample".karyotype.txt"


mkdir -p $imageDir/$group/$sample


echo "##### Obtain list of cromosomes (idiogram) for CIRCOS karyotype file"

awk '{if ($2 == 0 && $5 < 0.1 && $1 != "genome") {print "chr -", $1, $1, "0", $4, "id="$1}}' $sampleCoverage"_80" > $karyotypeFile

awk '{if ($2 == 0 && $5 < 0.2 && $1 != "genome") {print "chr -", $1, $1, "0", $4, "id="$1}}' $sampleCoverage"_80" > $karyotypeFile"_80"

#chr - NZ_CP015021.1 NZ_CP015021.1 0 81401  id=NZ_CP015021.1                                                           
#chr - NZ_CP015022.1 NZ_CP015022.1 0 95170  id=NZ_CP015022.1
#chr - NZ_CP015856.1 NZ_CP015856.1 0 115432  id=NZ_CP015856.1
#chr - NZ_CP017252.1 NZ_CP017252.1 0 92691  id=NZ_CP017252.1
#chr - NZ_CP017670.1 NZ_CP017670.1 0 92755  id=NZ_CP017670.1
#chr - NZ_CP018238.1 NZ_CP018238.1 0 91420  id=NZ_CP018238.1
#chr - NZ_CP018240.1 NZ_CP018240.1 0 94170  id=NZ_CP018240.1
#chr - NZ_CP018246.1 NZ_CP018246.1 0 95341  id=NZ_CP018246.1
#chr - NZ_CP018244.1 NZ_CP018244.1 0 92522  id=NZ_CP018244.1
#chr - NZ_CP018253.1 NZ_CP018253.1 0 95229  id=NZ_CP018253.1

echo "##### DONE Obtain list of cromosomes (idiogram) for CIRCOS karyotype file"


echo "##### Obtaining coverage coordinates from sequences"

bedtools genomecov -ibam $mappedDir/$group/$sample/$sample".sorted.bam" -bga -max 500 > $imageDir/$group/$sample/$sample".plasmid.bedgraph"

#NZ_CP011540.1   0       44983   0
#NZ_CP011540.1   44983   44985   1
#NZ_CP011540.1   44985   44989   3
#NZ_CP011540.1   44989   44990   5
#NZ_CP011540.1   44990   44991   6
#NZ_CP011540.1   44991   44992   7
#NZ_CP011540.1   44992   44993   15
#NZ_CP011540.1   44993   44994   16
#NZ_CP011540.1   44994   44999   32
#NZ_CP011540.1   44999   45004   33


echo "##### DONE obtaining coverage coordinates from sequences"

echo "##### Filtering coordinate file with matched sequences"
#obtain a list with matched plasmids to filter bedgraph and make drawing faster

listAcRegexp=$(printf "%s|" $(cat $plasmidMapped"_80") | sed 's/|$//g')

awk '{if ($1 ~ /'"${listAcRegexp}"'/) {print $1, $2, $3, $4}}' $imageDir/$group/$sample/$sample".plasmid.bedgraph" > $imageDir/$group/$sample/$sample".plasmid.bedgraph_80"

#rm $imageDir/$group/$sample/$sample"_plasmid.bedgraph"
echo "##### DONE Filtering coordinate file with matched sequences"
<<C
echo "##### Anotating contigs with Prokka"

prokka --force --outdir $imageDir/$group/$sample \
--prefix $sample \
--addgenes \
--kingdom Bacteria \
--usegenus \
--centre CNM \
--locustag $sample \
--compliant \
$contigDir/$group/$sample/contigs.fasta

#-genus Klebsiella \
#--species pneumoniae \
echo "##### DONE Anotating contigs with Prokka"
#From now on fasta sequences from contigs are those reformatted in prokka (program doesn't acept too large sequences names outputed by spades)

C
echo "##### Generatin contigs for  CIRCOS image"

makeblastdb -in $plasmidMappedFasta"_80_80" -out $imageDir/$group/$sample/$sample".blast.tmp" -dbtype nucl

blastn -query $imageDir/$group/$sample/$sample".fna" \
-db $imageDir/$group/$sample/$sample".blast.tmp" \
-out $imageDir/$group/$sample/$sample".plasmids.blast" -evalue 0.0001 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"


#Resulting file from blasi in outputfmt 6:
#gnl|CNM|A14_03_2	NZ_CP012884.1	100.00	122	0	0	1	122	133535	133414	2e-56	  226	484815	194742
#gnl|CNM|A14_03_2	NZ_CP012884.1	80.60	232	41	4	283550	283779	159259	159488	2e-41	  176	484815	194742
#gnl|CNM|A14_03_2	NZ_CP012884.1	87.60	121	15	0	1	121	46556	46436	6e-31	  141	484815	194742
#gnl|CNM|A14_03_2	NZ_CP018365.1	100.00	122	0	0	1	122	83332	83453	2e-56	  226	484815	260772
#gnl|CNM|A14_03_2	NZ_CP018365.1	100.00	121	0	0	1	121	10288	10408	6e-56	  224	484815	260772
#gnl|CNM|A14_03_2	NZ_CP018365.1	100.00	121	0	0	1	121	249067	249187	6e-56	  224	484815	260772
#gnl|CNM|A14_03_2	NZ_CP018365.1	80.60	232	41	4	283550	283779	242632	242861	2e-41	  176	484815	260772
#gnl|CNM|A14_03_2	NZ_CP018365.1	87.60	121	15	0	1	121	123441	123321	6e-31	  141	484815	260772
#gnl|CNM|A14_03_2	NZ_CP011990.1	100.00	121	0	0	1	121	89009	89129	6e-56	  224	484815	162533
#gnl|CNM|A14_03_2	NZ_CP011990.1	80.60	232	41	4	283550	283779	79464	79693	2e-41	  176	484815	162533


#Blast optput format 6 to contigs coordinates obtained with spades (Contig name = NODE_%d+%d+%d*_length%d_coverage_%d)
#awk 'BEGIN {OFS="\t"} {split($1, contigname, "_"); split($2, plasmidname, "|");print plasmidname[4], $9, $10, contigname[1]"_"contigname[2], "id="contigname[1]"_"contigname[2]}' \
#$imageDir/$group/$sample/$sample".plasmids.blast" > $imageDir/$group/$sample/$sample".plasmids.blast.contigs"


#Blast optput format 6 to CONTIGS coordinates obtained with prokka (Contig name = gnl|center(CNM)|locustag(sample)|contig_{%d}6)
#awk 'BEGIN {OFS="\t"} {split($1, contigname, "_");print $2, $9, $10, "contig_"contigname[length(contigname)], "id=contig_"contigname[length(contigname)]}' \

awk 'BEGIN {OFS="\t"} {split($1, contigname, "_"); {if ($3 > 90 && $4 > 500) print $2, $9, $10, "contig_"contigname[length(contigname)], "id=contig_"contigname[length(contigname)]}}' \
$imageDir/$group/$sample/$sample".plasmids.blast" > $imageDir/$group/$sample/$sample".plasmids.blast.contigs"

#NG_048025.1     1178    1476    contig_1    id=contig_1
#NG_048025.1     1361    1127    contig_1    id=contig_1
#NZ_CP010574.1   66599   66479   contig_2    id=contig_2
#NZ_CP008930.1   122307  122187  contig_2    id=contig_2
#NZ_CP011577.1   112628  112748  contig_2    id=contig_2
#NZ_CP006927.1   76749   76629   contig_2    id=contig_2
#NZ_CP008930.1   38761   37665   contig_3    id=contig_3
#NZ_CP008930.1   31098   30226   contig_3    id=contig_3
#NZ_CP010574.1   161444  160572  contig_3    id=contig_3
#NZ_CP006927.1   11494   10622   contig_3    id=contig_3

echo "##### DONE Generatin contigs for CIRCOS image"


echo "##### Generating link track for CIRCOS image"

#Parse blast output and obtain links between both datasets
##Blast optput format 6 to LINKS coordinates obtained with prokka (Contig name = gnl|center(CNM)|locustag|contig_{%d}6)

#awk 'BEGIN {OFS="\t"} {split($1, contigname, "_"); print "contig_"contigname[length(contigname)], $7,$8,$2,$9,$10, "id=contig_"contigname[length(contigname)]}' \

awk 'BEGIN {OFS="\t"} {split($1, contigname, "_"); {if ($3 > 90 && $4 > 500) print "contig_"contigname[length(contigname)], $7,$8,$2,$9,$10, "id=contig_"contigname[length(contigname)]}}' \
$imageDir/$group/$sample/$sample".plasmids.blast" > $imageDir/$group/$sample/$sample".plasmids.blast.links"

#sample.blast.links (file with matching coordinates between coordinates in contigs and plasmids)
#contig_1    915668  915969  NG_048025.1     1178    1476    id=contig_1
#contig_1    123340  123574  NG_048025.1     1361    1127    id=contig_1
#contig_2    1       121     NZ_CP010574.1   66599   66479   id=contig_2
#contig_2    1       121     NZ_CP008930.1   122307  122187  id=contig_2
#contig_2    1       121     NZ_CP011577.1   112628  112748  id=contig_2
#contig_2    1       121     NZ_CP006927.1   76749   76629   id=contig_2
#contig_3    98790   99886   NZ_CP008930.1   38761   37665   id=contig_3
#contig_3    472000  472872  NZ_CP008930.1   31098   30226   id=contig_3
#contig_3    472000  472872  NZ_CP010574.1   161444  160572  id=contig_3
#contig_3    472000  472872  NZ_CP006927.1   11494   10622   id=contig_3

##Change coordinates from contig --> plasmid to plasmid-->plasmid in order to represent them in CIRCOS

awk 'BEGIN{OFS="\t";}{if($1 != savedNode){savedNode= $1; delete chr} else{for(i in chr){print $4" "$5" "$6" "chr[i]" id="savedNode}}chr[$4$5$6] = $4" "$5" "$6}' \
$imageDir/$group/$sample/$sample".plasmids.blast.links" > $imageDir/$group/$sample/$sample".plasmids.blast.links.sorted"

#NG_048025.1 1361 1127 NG_048025.1 1178 1476 id=contig_1
#NZ_CP008930.1 122307 122187 NZ_CP010574.1 66599 66479 id=contig_2
#NZ_CP011577.1 112628 112748 NZ_CP010574.1 66599 66479 id=contig_2
#NZ_CP011577.1 112628 112748 NZ_CP008930.1 122307 122187 id=contig_2
#NZ_CP006927.1 76749 76629 NZ_CP011577.1 112628 112748 id=contig_2
#NZ_CP006927.1 76749 76629 NZ_CP010574.1 66599 66479 id=contig_2
#NZ_CP006927.1 76749 76629 NZ_CP008930.1 122307 122187 id=contig_2
#NZ_CP008930.1 31098 30226 NZ_CP008930.1 38761 37665 id=contig_3
#NZ_CP010574.1 161444 160572 NZ_CP008930.1 38761 37665 id=contig_3
#NZ_CP010574.1 161444 160572 NZ_CP008930.1 31098 30226 id=contig_3

echo "##### DONE Generating link track for CIRCOS image"


echo "##### Creating annotated track from contigs assembled (annotated with Prokka)"


#Filter Gff(3) file from prokka and create a bed file with coordinates with annotated genes (WITH NAMES)
awk 'BEGIN{OFS="\t";}{split($1,query,"_");split($9,description,"Name=");split(description[2],name,";");{if ($3 == "gene" && $1 != "#" && $9 ~ /Name/) {print "contig_"query[length(query)],$4,$5, name[1]}}}' \
$imageDir/$group/$sample/$sample".gff" > $imageDir/$group/$sample/$sample".gff.bed"

#Create a file with matches between file .bed and file .blast.link based on contig name
#$sample.bed (file with annotation from prokka on contigs assembled)
#contig_1    1005    2300    kgtP_1
#contig_1    2669    4024    pssA
#contig_1    7599    8024    trxC
#contig_1    9371    10060   ung
#contig_1    10373   10756   grcA
#contig_1    10802   12133   srmB
#contig_1    12265   13002   yfiC
#contig_1    12987   14606   nadB
#contig_1    15027   15602   rpoE
#contig_1    15635   16285   rseA


echo "##### Adapting coords in annotated contigs"

awk 'NR==FNR{a[NR]=$1;b[NR]=$0;next}{for(i = 1; i <= NR; ++i){if (a[i] == $1) print b[i],"\t", $0}}' \
$imageDir/$group/$sample/$sample".gff.bed" $imageDir/$group/$sample/$sample".plasmids.blast.links" > $imageDir/$group/$sample/$sample".bed.coords.tmp"

#The result file is a combination of annotation and coordinates:
#contig_1    1005    2300    kgtP_1   contig_1   915668  915969  NG_048025.1     1178    1476    id=contig_1
#contig_1    2669    4024    pssA     contig_1   915668  915969  NG_048025.1     1178    1476    id=contig_1
#contig_1    7599    8024    trxC     contig_1   915668  915969  NG_048025.1     1178    1476    id=contig_1
#contig_1    9371    10060   ung      contig_1   915668  915969  NG_048025.1     1178    1476    id=contig_1
#contig_1    10373   10756   grcA     contig_1   915668  915969  NG_048025.1     1178    1476    id=contig_1
#contig_1    10802   12133   srmB     contig_1   915668  915969  NG_048025.1     1178    1476    id=contig_1
#contig_1    12265   13002   yfiC     contig_1   915668  915969  NG_048025.1     1178    1476    id=contig_1
#contig_1    12987   14606   nadB     contig_1   915668  915969  NG_048025.1     1178    1476    id=contig_1
#contig_1    15027   15602   rpoE     contig_1   915668  915969  NG_048025.1     1178    1476    id=contig_1
#contig_1    15635   16285   rseA     contig_1   915668  915969  NG_048025.1     1178    1476    id=contig_1


#If annotation is within blasted fragment of the contig, coordinates are translated from contig to plasmid in order to represent them on CIRCOS

awk '{if (($2 >= $6 - 7000 && $2 <= $7) || ($3 >= $6 && $3 <= $7 + 7000)){{isInverted=($10-$9); \
geneLenght=($3-$2);{if (isInverted < 0) {geneLenght *= -1}}; \
coordChr1=(($2-$6)+$9);coordChr2=((($2-$6)+$9)+geneLenght)}\
{print $8, coordChr1, coordChr2, $4}}}' $imageDir/$group/$sample/$sample".bed.coords.tmp" > $imageDir/$group/$sample/$sample".bed.all.coords"  

#resulting in a bed file with coordinated of plasmid bur refering to contig annotation:
#NZ_CP010574.1 34820 33528 arsB_1
#NZ_CP008930.1 90527 89235 arsB_1
#NZ_CP006927.1 44969 43677 arsB_1
#NZ_CP010574.1 81021 82508 ltrA_1
#NZ_CP008930.1 144220 145707 ltrA_1
#NZ_CP011577.1 91475 89988 ltrA_1
#NG_048025.1 101 3253 bepE_2
#NZ_CP008930.1 146692 147117 klcA_2
#NZ_CP008930.1 147354 147608 holE_2
#NZ_CP008930.1 152085 152627 ssb_2


#Remove duplicate of several matches 

awk '{uniq[$1$4]=$0}END{for (i in uniq) print uniq[i]}' /$imageDir/$group/$sample/$sample".bed.all.coords" \
> $imageDir/$group/$sample/$sample".bed.coords"

echo "##### DONE Adapting coords in annotated contigs"
echo "##### DONE Creating annotated track from contigs assembled (annotated with Prokka)"



echo "##### Creating antibiotic resistance track from contigs assembled"

#blast contigs on ARGannot.fasta DDBB

makeblastdb -in $abReference -out /processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PS_MP_T/REFERENCES/ARGannot -dbtype nucl

blastn -query $imageDir/$group/$sample/$sample".fna" \
-db /processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PS_MP_T/REFERENCES/ARGannot \
-out $imageDir/$group/$sample/$sample".abr.blast" -evalue 0.0001 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"

#/opt/srst2-0.1.8/data/ARGannot.fasta

#gnl|CNM|ECLCNM_contig_1     0__OqxBgb_Flq__OqxBgb__49       73.53   306     70      9       915668  915969  1078    1376    1e-20     106
#gnl|CNM|ECLCNM_contig_1     0__OqxBgb_Flq__OqxBgb__49       73.75   240     53      9       123340  123574  1261    1027    2e-14   86.1
#gnl|CNM|ECLCNM_contig_7     86__AmpH_Bla__AmpH__634 99.22   1161    9       0       62169   63329   1161    1       0.0      2095
#gnl|CNM|ECLCNM_contig_14     371__CatBx_Phe__CatB4__602      100.00  108     0       0       96486   96593   442     549     6e-50     200
#gnl|CNM|ECLCNM_contig_17     0__OqxBgb_Flq__OqxBgb__49       98.86   3153    36      0       31982   35134   1       3153    0.0      5624
#gnl|CNM|ECLCNM_contig_17     77__OqxA_Flq__OqxA__1044        99.23   1176    9       0       30783   31958   1       1176    0.0      2122
#gnl|CNM|ECLCNM_contig_36     193__CTX-M-1_Bla__CTX-M-15__150 100.00  876     0       0       1073    1948    876     1       0.0      1618
#gnl|CNM|ECLCNM_contig_36     193__CTX-M-1_Bla__CTX-M-71__176 99.89   876     1       0       1073    1948    876     1       0.0      1613
#gnl|CNM|ECLCNM_contig_36     193__CTX-M-1_Bla__CTX-M-57__167 99.89   876     1       0       1073    1948    876     1       0.0      1613
#gnl|CNM|ECLCNM_contig_36     193__CTX-M-1_Bla__CTX-M-55__166 99.89   876     1       0       1073    1948    876     1       0.0      1613



#Blast optput format 6 to BED coordinates obtained with prokka (Contig name = gnl|center(CNM)|locustag|contig_{%d}6)
#and matching ARGannot DDBB
#Output = contig name - coordinates on contig - name of ABR gene
#Output is unique since ARGannot DDBB has several serotypes of AR genes
cat $imageDir/$group/$sample/$sample".abr.blast" |sort -k 12 -nr |\
awk '{OFS="\t";split($1, contigname, "_"); split($2, ABR, "__")};($3 > 95)&&(($4/$14)>0.8)&&(!x[ABR[2]contigname[length(contigname)]]++){print "contig_"contigname[length(contigname)], $7, $8, ABR[2]"_"ABR[3]}' \
> $imageDir/$group/$sample/$sample".abr.bed"

#contig_17 31982 35134 OqxBgb_Flq
#contig_55 809 1993 TetD_Tet
#contig_17 30783 31958 OqxA_Flq
#contig_7 62169 63329 AmpH_Bla
#contig_36 1073 1948 CTX-M-1_Bla
#contig_47 1045 1905 SHV-OKP-LEN_Bla
#contig_46 3818 4678 TEM-1D_Bla
#contig_60 1818 2678 Aac3-IIa_AGly
#contig_46 2279 3118 OXA-9_Bla
#contig_61 1658 2455 OXA-48_Bla


#Transfors contigs coord to plasmid coords 

#Create a intermediate file with all matches in contig name
awk 'NR==FNR{a[NR]=$1;b[NR]=$0;next}{for(i = 1; i <= NR; ++i){if (a[i] == $1) print b[i],"\t", $0}}' \
$imageDir/$group/$sample/$sample".abr.bed" \
$imageDir/$group/$sample/$sample".plasmids.blast.links" \
> $imageDir/$group/$sample/$sample".abr.bed.coords.tmp"

#contig_7	168907	170064	AMPH_Ecoli_Bla 	 contig_7	74669	75873	NC_019114.1	21471	20266	id=contig_7
#contig_7	168907	170064	AMPH_Ecoli_Bla 	 contig_7	74676	75870	NZ_CP008906.1		139589	140784	id=contig_7
#contig_7	168907	170064	AMPH_Ecoli_Bla 	 contig_7	172420	172671	CP019443.1		26728	26477	id=contig_7
#contig_7	168907	170064	AMPH_Ecoli_Bla 	 contig_7	172420	172671	NZ_CP011062.1	234981	235232	id=contig_7
#contig_8	45602	47227	Mcr1_Colistin 	 contig_8	101753	157404	NZ_CP015833.1	13262	68917	id=contig_8
#contig_8	45602	47227	Mcr1_Colistin 	 contig_8	47435	89238	NZ_CP015833.1	198847	240662	id=contig_8
#contig_8	45602	47227	Mcr1_Colistin 	 contig_8	1		39667	NZ_CP015833.1	151300	190965	id=contig_8
#contig_8	45602	47227	Mcr1_Colistin 	 contig_8	89239	103015	NZ_CP015833.1	1		13793	id=contig_8
#contig_8	45602	47227	Mcr1_Colistin 	 contig_8	39664	44812	NZ_CP015833.1	192408	197555	id=contig_8
#contig_8	45602	47227	Mcr1_Colistin 	 contig_8	102759	104082	NZ_CP015833.1	12066	13387	id=contig_8


#If annotation is within blasted fragment of the contig, coordinates are translated from contig to plasmid in order to represent them on CIRCOS


awk '{if (($2 >= $6 - 7000 && $2 <= $7) || ($3 >= $6 && $3 <= $7 + 7000)){{isInverted=($10-$9); geneLength=($3-$2);{if (isInverted < 0) {geneLength *= -1}}; coordChr1=(($2-$6)+$9);coordChr2=((($2-$6)+$9)+ geneLength)} {print $8, coordChr1, coordChr2, $4}}}' \
$imageDir/$group/$sample/$sample".abr.bed.coords.tmp" \
> $imageDir/$group/$sample/$sample".abr.bed.all.coords"

#NG_048025.1 1178 1479 OqxBgb_Flq
#NG_048025.1 1361 1127 OqxBgb_Flq
#NZ_CP010574.1 143152 143045 CatBx_Phe
#NZ_CP008930.1 37212 37105 CatBx_Phe
#NZ_CP011577.1 21393 21286 CatBx_Phe
#NG_048025.1 101 3253 OqxBgb_Flq
#NG_048987.1 1177 302 CTX-M-1_Bla
#NG_048987.1 1178 304 CTX-M-1_Bla
#NG_048987.1 1180 308 CTX-M-1_Bla
#NG_048987.1 1187 463 OXY_Bla

#Remove duplicate of several matches 

awk '{uniq[$1$4]=$0}END{for (i in uniq) print uniq[i]}' /$imageDir/$group/$sample/$sample".abr.bed.all.coords" \
> $imageDir/$group/$sample/$sample".abr.bed.coords"

echo "##### DONE Creating antibiotic resistance track from contigs assembled"


echo "##### Creating replisome track from contigs assembled (PlasmidFinder DDBB)"

#blast contigs on PlasmidFinder DDBB

makeblastdb -in $plasmidFinderReference -out /processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PS_MP_T/REFERENCES/PLASMIDS/PLASMIDFINDER -dbtype nucl

blastn -query $imageDir/$group/$sample/$sample".fna" \
-db /processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PS_MP_T/REFERENCES/PLASMIDS/PLASMIDFINDER \
-out $imageDir/$group/$sample/$sample".pfinder.blast" -evalue 0.0001 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"


#Obtain bed coords in contigs
#Filter by identity greater than 95%
#Get the first ocurrence of each contig since only one replisome is expected per contig
cat $imageDir/$group/$sample/$sample".pfinder.blast" | sort -k 12 -rn \
| awk '{OFS="\t"; split($1,contigname,"_"); split($2,repname,"_")};($3 > 95)&&(($4/$14)>0.8)&&(!contig[repname[1]"_"repname[2]contigname[length(contigname)]]++){print "contig_"contigname[length(contigname)], $7, $8, repname[1]"_"repname[2]}' \
> $imageDir/$group/$sample/$sample".pfinder.bed"

#contig_107      4433    5114    IncFIB(AP001918)_1
#contig_14       124847  125476  IncHI2A_1
#contig_57       11563   12205   IncFIB(S)_1
#contig_52       46975   47301   IncHI2_1
#contig_116      122     389     IncFIA_1
#contig_190      1595    1853    IncFII_1

#Transfors contigs coord to plasmid coords 

#Create a intermediate file with all matches in contig name
awk 'NR==FNR{a[NR]=$1;b[NR]=$0;next}{for(i = 1; i <= NR; ++i){if (a[i] == $1) print b[i],"\t", $0}}' \
$imageDir/$group/$sample/$sample".pfinder.bed" \
$imageDir/$group/$sample/$sample".plasmids.blast.links" \
> $imageDir/$group/$sample/$sample".pfinder.bed.coords.tmp"


#contig_25	59506	59653	IncFII(K)_1 	 contig_25	12675	62112	NZ_CP015132.1	92402	42983	id=contig_25
#contig_25	59506	59653	IncFII(K)_1 	 contig_25	6277	11583	NZ_CP015132.1	97699	92393	id=contig_25
#contig_25	59506	59653	IncFII(K)_1 	 contig_25	2947	6212	NZ_CP015132.1	101088	97819	id=contig_25
#contig_25	59506	59653	IncFII(K)_1 	 contig_25	2	2826	NZ_CP015132.1	103977	101153	id=contig_25
#contig_25	59506	59653	IncFII(K)_1 	 contig_25	16775	17556	NZ_CP015132.1	98354	97573	id=contig_25
#contig_25	59506	59653	IncFII(K)_1 	 contig_25	14171	14797	NZ_CP015132.1	98170	97544	id=contig_25
#~contig_25	59506	59653	IncFII(K)_1 	 contig_25	14171	14768	NZ_CP015132.1	88131	87534	id=contig_25
#contig_25	59506	59653	IncFII(K)_1 	 contig_25	16959	17556	NZ_CP015132.1	90912	90322	id=contig_25
#contig_25	59506	59653	IncFII(K)_1 	 contig_25	2407	3069	NZ_CP015132.1	88204	87537	id=contig_25
#contig_25	59506	59653	IncFII(K)_1 	 contig_25	5792	6409	NZ_CP015132.1	101573	100956	id=contig_25


#If annotation is within blasted fragment of the contig, coordinates are translated from contig to plasmid in order to represent them on CIRCOS


awk '{if (($2 >= $6 - 500 && $2 <= $7) || ($3 >= $6 && $3 <= $7 + 500)){{isInverted=($10-$9); geneLenght=($3-$2);{if (isInverted < 0) {geneLenght *= -1}}; \
coordChr1=(($2-$6)+$9);coordChr2=((($2-$6)+$9)+geneLenght)} {print $8, coordChr1, coordChr2, $4}}}' \
$imageDir/$group/$sample/$sample".pfinder.bed.coords.tmp" \
> $imageDir/$group/$sample/$sample".pfinder.bed.coords"

#NZ_CP015132.1 139233 139086 IncFII(K)_1
#NZ_CP012884.1 125294 125441 IncFII(K)_1
#NZ_CP018365.1 207768 207915 IncFII(K)_1
#NZ_CP018430.1 165980 166127 IncFII(K)_1
#NC_009649.1 5612 5759 IncFII(K)_1
#NZ_CP008800.1 194719 194866 IncFII(K)_1
#NZ_CP015386.1 171950 171803 IncFII(K)_1
#NZ_CP018355.1 104714 104861 IncFII(K)_1
#NZ_CP011977.1 218678 218825 IncFII(K)_1
#NZ_CP008930.1 17747 17894 IncFII(K)_1

echo "##### DONE Creating replisome track from contigs assembled (PlasmidFinder DDBB)"



echo "##### Creating IS track from contigs assembled (ISFinder DDBB)"

#blast contigs on ISFinder DDBB

makeblastdb -in $isFinderReference -out /processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PS_MP_T/REFERENCES/PLASMIDS/ISFINDER -dbtype nucl

blastn -query $imageDir/$group/$sample/$sample".fna" \
-db /processing_Data/bioinformatics/research/20160530_ANTIBIOTICS_PS_MP_T/REFERENCES/PLASMIDS/ISFINDER \
-out $imageDir/$group/$sample/$sample".isfinder.blast" -evalue 0.0001 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"


#Obtain bed coords in contigs
#Filter by identity greater than 95%
#Get the first ocurrence of each contig since only one replisome is expected per contig
cat $imageDir/$group/$sample/$sample".isfinder.blast" | sort -k 12 -rn \
| awk '{OFS="\t"; split($1,contigname,"_"); split($2,repname,"_")};($3 > 95)&&($4>500)&&(($4/$13)>0.5)&&(($12)>1000)&&(!contig[isname[1]"_"isname[2]contigname[length(contigname)]]++){print "contig_"contigname[length(contigname)], $7, $8, $2}' \
> $imageDir/$group/$sample/$sample".isfinder.bed"

#| awk '{OFS="\t"; split($1,contigname,"_"); split($2,repname,"_")};($3 > 95)&&(($4)>500)&&(!contig[isname[1]"_"isname[2]contigname[length(contigname)]]++){print "contig_"contigname[length(contigname)], $7, $8, $2}' \


#contig_107      4433    5114    IncFIB(AP001918)_1
#contig_14       124847  125476  IncHI2A_1
#contig_57       11563   12205   IncFIB(S)_1
#contig_52       46975   47301   IncHI2_1
#contig_116      122     389     IncFIA_1
#contig_190      1595    1853    IncFII_1

#Transfors contigs coord to plasmid coords 

#Create a intermediate file with all matches in contig name
awk 'NR==FNR{a[NR]=$1;b[NR]=$0;next}{for(i = 1; i <= NR; ++i){if (a[i] == $1) print b[i],"\t", $0}}' \
$imageDir/$group/$sample/$sample".isfinder.bed" \
$imageDir/$group/$sample/$sample".plasmids.blast.links" \
> $imageDir/$group/$sample/$sample".isfinder.bed.coords.tmp"


#contig_25	59506	59653	IncFII(K)_1 	 contig_25	12675	62112	NZ_CP015132.1	92402	42983	id=contig_25
#contig_25	59506	59653	IncFII(K)_1 	 contig_25	6277	11583	NZ_CP015132.1	97699	92393	id=contig_25
#contig_25	59506	59653	IncFII(K)_1 	 contig_25	2947	6212	NZ_CP015132.1	101088	97819	id=contig_25
#contig_25	59506	59653	IncFII(K)_1 	 contig_25	2	2826	NZ_CP015132.1	103977	101153	id=contig_25
#contig_25	59506	59653	IncFII(K)_1 	 contig_25	16775	17556	NZ_CP015132.1	98354	97573	id=contig_25
#contig_25	59506	59653	IncFII(K)_1 	 contig_25	14171	14797	NZ_CP015132.1	98170	97544	id=contig_25
#~contig_25	59506	59653	IncFII(K)_1 	 contig_25	14171	14768	NZ_CP015132.1	88131	87534	id=contig_25
#contig_25	59506	59653	IncFII(K)_1 	 contig_25	16959	17556	NZ_CP015132.1	90912	90322	id=contig_25
#contig_25	59506	59653	IncFII(K)_1 	 contig_25	2407	3069	NZ_CP015132.1	88204	87537	id=contig_25
#contig_25	59506	59653	IncFII(K)_1 	 contig_25	5792	6409	NZ_CP015132.1	101573	100956	id=contig_25


#If annotation is within blasted fragment of the contig, coordinates are translated from contig to plasmid in order to represent them on CIRCOS


awk '{if (($2 >= $6 - 500 && $2 <= $7) || ($3 >= $6 && $3 <= $7 + 500)){{isInverted=($10-$9); geneLenght=($3-$2);{if (isInverted < 0) {geneLenght *= -1}}; \
coordChr1=(($2-$6)+$9);coordChr2=((($2-$6)+$9)+geneLenght)} {print $8, coordChr1, coordChr2, $4}}}' \
$imageDir/$group/$sample/$sample".isfinder.bed.coords.tmp" \
> $imageDir/$group/$sample/$sample".isfinder.bed.coords"

#NZ_CP015132.1 139233 139086 IncFII(K)_1
#NZ_CP012884.1 125294 125441 IncFII(K)_1
#NZ_CP018365.1 207768 207915 IncFII(K)_1
#NZ_CP018430.1 165980 166127 IncFII(K)_1
#NC_009649.1 5612 5759 IncFII(K)_1
#NZ_CP008800.1 194719 194866 IncFII(K)_1
#NZ_CP015386.1 171950 171803 IncFII(K)_1
#NZ_CP018355.1 104714 104861 IncFII(K)_1
#NZ_CP011977.1 218678 218825 IncFII(K)_1
#NZ_CP008930.1 17747 17894 IncFII(K)_1

echo "##### DONE Creating IS track from contigs assembled (ISFinder DDBB)"















<<Trials


awk'{if($1 != savedNode){savedNode= $1; delete chr} else{for(i in chr){print $4" "$5" "$6" "chr[i]" "savedNode}}chr[$4$5$6] = $4" "$5" "$6}' $sample".plasmids.blast.contigs" >



#awk 'split($1, ac,"_") {p[$1]="chr - "ac[1]"_"ac[2]" "ac[1]"_"ac[2]" 0 "ac[4]" id="ac[1]"_"ac[2]} END {for (i in p){print p[i]}}' $imageDir/$group/$sample/$sample".plasmids.blast" >> $karyotypeFile"contig"



Trials




#cat A14_03_plasmid.links |awk '{vals[$1] =vals[$1]"," $4" "$5" "$6} END {for (i in vals) {print i, vals[i]}}' > A14_03_plasmid.links.sorted

#awk 'BEGIN {RS=","} {print $0}'

#awk '{p[$1]=p[$1]","$4" "$5" "$6} END {for (i in p){ split(p[i], chr, ","); for(j in chr){for(z in chr){ if(z>j && chr[j]!= "")print chr[j]" "chr[z]}}}}' A14_03_plasmid.links 


#awk'{if($1 != savedNode){savedNode= $1; delete chr} else{for(i in chr){print $4" "$5" "$6" "chr[i]" "savedNode}}chr[$4$5$6] = $4" "$5" "$6}'
