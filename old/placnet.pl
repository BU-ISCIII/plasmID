#!/usr/bin/perl  
#v1.03

use strict;
#use warnings;
use Getopt::Long;


### Variables

my @Options;
my $inputFile;
my $prefix;
my @inputText;
my $contigsFile;
my $refDBFile;
my @sams;
my @DBs;
my $MINLENGTH;
my $blastNet;
my $scaffNet ="";
my $fastaProt;
my $fastaNucl;
my $dbRefType;
my $db;

### temp variables
my $line;
my @fields;
my @c;
my $l;
my $s;

@Options = (
	{OPT=>"f=s",	VAR=>\$inputFile,	DESC=>"Input file for PLACNET analysis"},
	{OPT=>"p=s",	VAR=>\$prefix,	DEFAULT => 'placnet', DESC=>"Prefix name for output files"},
	{OPT=>"l=s",	VAR=>\$MINLENGTH,	DEFAULT => '200', DESC=>"Minimun contig length to scaffold process"}
	
);

#Check options and set variables



(@ARGV < 1) && (usage());
if ($ARGV[0] eq "-generate")
{
	generateInputFile();
	exit 0;
}


GetOptions(map {$_->{OPT}, $_->{VAR}} @Options) || usage();

# Now setup default values.
foreach (@Options) {
	if (defined($_->{DEFAULT}) && !defined(${$_->{VAR}})) {
	${$_->{VAR}} = $_->{DEFAULT};
	}
}



unless($inputFile){
	print STDERR "You must specified the input files\n";
	&usage();
}

open(INPUT,$inputFile);
@inputText = <INPUT>;
close INPUT;

print @inputText;



####   PARSING INPUT FILE   #############

foreach $line (@inputText)
{
	if($line =~ /^CONTIGS/)
	{
		chomp $line;
		@fields = split('\t',$line);
		$contigsFile = $fields[1];
	}
	if($line =~ /^REFDB/)
	{
		chomp $line;
		@fields = split('\t',$line);
		$refDBFile = $fields[1];
		$dbRefType = $fields[2];
	}
	if($line =~ /^SAM/)
	{
		chomp $line;
		$line =~ s/SAM:\t//g;
		push(@sams,$line);
	}
	if($line =~ /^DB/)
	{
		chomp $line;
		$line =~ s/DB\d*:\t//g;
		push(@DBs,$line);
	}
	
}

###   CHECKING INPUT FILE   ##########

if($contigsFile eq "")
{
	print "\nError in $inputFile file, Please specify a contig file\n\n";
	exit;
}
else{ 
	
	system("gmhmmp_heuristic.pl -s $contigsFile -out tmpGM -a");
	$fastaProt = gm2fasta("tmpGM.lst");
	system("gmhmmp_heuristic.pl -s $contigsFile -out tmpGM -d");
	$fastaNucl = gm2fasta("tmpGM.lst");
	
	#### CDS and ORF prediction ######
	open(FPROT,">$prefix.gm.faa");
	print FPROT $fastaProt;
	close FPROT;
	open(FNUCL, ">$prefix.gm.cds");
	print FNUCL $fastaNucl;
	close FNUCL;
}
		
#}
if($refDBFile eq "")
{
	print "\nError in $inputFile file, Please specify a reference database file\n\n";
	exit 0;
}
if(scalar(@sams)<1)
{
	print "\nError in $inputFile file, Please specify a SAM file file\n\n";
	exit 0;
}

############## MAIN ROUTINE####################



print "\n\nMaking blast against $refDBFile .........\n\n";
$blastNet = blastRefDB($dbRefType);



foreach $s (@sams)
{
	chomp $s;
	@c = split('\t',$s);
	print "\nMaking scaffold network from $c[0]\n\n";
	
	$scaffNet .= sam2scaffold($s);
}


open(OUT,">$prefix.net.csv");
print OUT $blastNet;
print OUT $scaffNet;
close OUT;

foreach $db (@DBs)
{
		@fields = split('\t',$db);	
		database($fields[0],$fields[1],$fields[2],$fields[3],$fields[4]);
}













################################################


#### BLASTING REFDB   #########

sub blastRefDB     #### blastRefDB(type)
{
	my $node;
	my $value;
	my $n;
	my $l;
	my $out;
	my $type = $_[0];
	
	if($type eq "blast")
	{
		system("blastn -query $contigsFile -db $refDBFile -out tmpMegaBlast.txt -num_alignments 0 -evalue 1e-25");
	}elsif ($type eq "fasta")
	{
		system("makeblastdb -in $refDBFile -out tmpRefDB -dbtype nucl");
		system("blastn -query $contigsFile -db tmpRefDB -out tmpMegaBlast.txt -num_alignments 0 -evalue 1e-25");
	}else{
		print "Error in Reference DB format\n";
		exit 0;
	}
	open(A,"tmpMegaBlast.txt");
	my @mega = <A>;
	close A;
	
	foreach $l (@mega)
	{
		$l =~ s/;//g;
		if ($l =~ /Query= (.*)\n/)
		{
			$node = $1;
			$value =0;
			$n=1;
			#print "$node\n";
		}
		if ($l =~ /gi\|/) 
		{
			#print $l;
			@c = split(' ',$l);
	
			if ($value == 0)
				{
				$value = $c[-2];
			}
			if($c[-2]/$value > 0.85+(0.02*$n))
			{
				$n++;
				$out .= "$node\thit\t$c[0]\t$c[-2]\t$c[-1]\n";
				
				$value = ($value + $c[-2])/2;
			}
		}
	}
	return $out;
}


##### SAM2SCAFFOLD subrutine ###############
sub sam2scaffold    #### attr:	sam2scaffold(SamDefinition)
{
	### attributes
	$s = $_[0];
	chomp $s;
	@c = split('\t',$s);
	
	my $samFile = $c[0];
	my $readLength = $c[1];
	my $insert = $c[2];
	
	print "$c[0]\t$c[1]\t$c[2]\n";
	
	my $FLANK = $insert + 2*$readLength;
	
	
	
	### Variables
	my @txt;
	my $l;
	my $k;
	my $j;
	my $i;
	my %contigs;
	my %length;
	my %scaffolds;
	my $sum;
	my $avg;
	my $out = "";
	my $header;
	my @ks;
	
	open(CONTIGS,$contigsFile);
	@txt = <CONTIGS>;
	close CONTIGS;

	foreach $l (@txt)
	{
		chomp $l;
		if($l =~ /^>(.*)/)
		{
			$header = $1;
			$header =~ s/\s//g;
			
		}else{
			$contigs{$header} .= $l;
		}
	}
	foreach $k (keys(%contigs))
	{
		$contigs{$k} =~ s/\s//g;
		$length{$k} = length($contigs{$k});
				
	}
	
	
	open(SAM,$samFile);
	
	
### SAM Format specifications	
	#$c[0]	1	QNAME	Query template/pair NAME
	#$c[1]	2	FLAG	bitwise FLAG
	#$c[2]	3	RNAME	Reference sequence NAME
	#$c[3]	4	POS	1-based leftmost POSition/coordinate of clipped sequence
	#$c[4]	5	MAPQ	MAPping Quality (Phred-scaled)
	#$c[5]	6	CIAGR	extended CIGAR string
	#$c[6]	7	MRNM	Mate Reference sequence NaMe (‘=’ if same as RNAME)
	#$c[7]	8	MPOS	1-based Mate POSistion
	#$c[8]	9	TLEN	inferred Template LENgth (insert size)
	#$c[9]	10	SEQ	query SEQuence on the same strand as the reference
	#$c[10]	11	QUAL	query QUALity (ASCII-33 gives the Phred base quality)
	#$c[11]	12+	OPT	variable OPTional fields in the format TAG:VTYPE:VALUE
####
	
	while ($l=<SAM>)
	{
		@c = split('\t',$l);
		if ($c[6] !~ /\*|=/ & $length{$c[6]} > $MINLENGTH & $length{$c[2]} > $MINLENGTH ) ### Only contigs longer than MINLENGTH
		{
			
			if($c[3] <= $FLANK | abs($length{$c[2]}-$c[3]) <= $FLANK && ($c[7] <= $FLANK | abs($length{$c[6]}-$c[7]) <= $FLANK))
			{
				
				$scaffolds{$c[2]}{$c[6]}++;
				$scaffolds{$c[6]}{$c[2]}++;
				
			}
		}
		
	}
	close SAM;
	
	$sum =0;
	foreach $k (keys(%length))
	{
		foreach $j (keys(%length))
		{
			$sum += $scaffolds{$k}{$j};
		}
	}
	
	$avg = ($sum/2)/(1+scalar(keys(%scaffolds)));  ##### $sum/2 for simetric matrix
	
	
	@ks = keys(%length);
	for($i=0; $i<scalar(@ks); $i++)
	{
		for($j=$i+1; $j<scalar(@ks);$j++)
		{
				if($scaffolds{$ks[$i]}{$ks[$j]} > $avg*0.30)   #### Only reports scaffold with abundance over 30% of the mean scaffold abundance.
				{
					$out .= "$ks[$i]\tscaff\t$ks[$j]\t\t\t$scaffolds{$ks[$i]}{$ks[$j]}\n";
				}
		}
	}
	return $out;
	
}

sub database   #### Attributes name,fastaFile,type,threshold
{
	my $name = $_[0];
	my $dbFastaFile = $_[1];
	my $dbType = $_[2];
	my $evalue = $_[3];
	my $formatDB =$_[4];
	my @dbText;
	my $line;
	my @fields1;
	my @fields2;
	my $db = "tmpTypeDB";
	
	
	if($formatDB eq "fasta")
	{
		
		system("makeblastdb -in $dbFastaFile -out $db -dbtype $dbType");
		
	}elsif ($formatDB eq "blast")
	{
		$db = $dbFastaFile;
	}else{
		print "Error: unvalid format of $name database.\n";
		exit 0;
	}
	
	if($dbType eq "prot")
	{
		system("blastp -query $prefix.gm.faa -db $db -outfmt 6 -evalue $evalue -out tmp$name.blast -num_alignments 1"); 
	}elsif ($dbType eq "nucl")
	{
		system("blastn -query $contigsFile -db $db -outfmt 6 -evalue $evalue -out tmp$name.blast -num_alignments 1");
	}
	
	open(DB,"tmp$name.blast");
	@dbText = <DB>;
	close DB;

	open(OUT,">$name.annot");
	
	if($dbType eq "prot")
	{
		foreach $line (@dbText)
		{
			@fields1 = split('\t',$line);
			@fields2 = split('\|',$fields1[0]);
			print OUT "$fields2[1]\t$fields1[1]\n";
		}
		close OUT;
	}else{
		foreach $line (@dbText)
		{
			@fields1 = split('\t',$line);
			print OUT "$fields1[0]\t$fields1[1]\n";
		}
		close OUT;	
	}	
}
	
sub gm2fasta ############## attr: (geneMarkOutput.lst) return: fasta
{
	
	open(A,"tmpGM");
	my @txt = <A>;
	close A;


	my $out ="";
	my $cond=0;
	foreach $line (@txt)
	{
		if($line =~ />/)
		{
			$line =~ s/\|GeneMark\.hmm\|\d+_(aa|nt)\|(\-|\+)\|\d+\|\d+\t>/\|/;
		}
		if($line =~ /#===/)
		{
			$cond=0;
		}
		if ($cond==1)
		{
			$out .= $line;
		}
		if($line =~ /Predicted proteins:/ | $line =~ /Nucleotide sequence of predicted genes:/)
		{
		   $cond=1;
		}
	}
	
	
	return $out;
	
}

sub usage
{
	print "Usage:\n\n";
	print "Placnet v1.03 10/06/2015\n";
	print "writen by: Val F. Lanza (valfernandez.vf\@gmail.com) and Maria de Toro (mdtorohernando\@gmail.com\n\n";\
	print "Please cite PLACNET as: \nLanza VF, de Toro M, Garcillán-Barcia MP, Mora A, Blanco J, Coque TM, de la Cruz F: \nPlasmid Flux in Escherichia coli ST131 Sublineages, Analyzed by Plasmid Constellation Network (PLACNET),\na New Method for Plasmid Reconstruction from Whole Genome Sequences. \nPLoS Genet 2014, 10:e1004766\n\n";
	print "Write inputFile Template\n\nplacnet.pl -generate\n\n";
	print "Network process\n\nplacnet.pl -f inputFile.txt -p prefix -l min\n\n\n";

	foreach (@Options) {
		
		printf "  -%-13s %s%s.\n",$_->{OPT},$_->{DESC},
			defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
	}
	print "\n\n\n";
	exit(1);
}

sub generateInputFile
{
	open(OUT,">inputFile.txt");
	my $text = "### Placnet Input File ###



CONTIGS:	multifasta.fasta
SAM:	file1.sam	readLength1	insertSize1
SAM:	file2.sam	readLength2	insertSize2
#...

REFDB:	refdb	type(fasta/blast)


##### Optional Attributes

DB1:	name	file.fasta	type(nucl/prot)	threshold(E value)	format (fasta/blast)
DB2:	name	file.fasta	type(nucl/prot)	threshold(E value)	format (fasta/blast)
DB3:	name	file.fasta	type(nucl/prot)	threshold(E value)	format (fasta/blast)
#....";
	print OUT $text;
	close OUT;
}
