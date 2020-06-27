#!/usr/bin/perl -w
use strict;
die "perl $0 <genomic.gff|genomic.gff.gz> <all.cds|all.cds.gz>\n" unless @ARGV == 2;

##### example inputs: GCA_002891405.2_Csec_1.0_genomic.gff.gz, GCA_002891405.2_Csec_1.0_cds_from_genomic.fna.gz #####

my $in = shift;
my $pool = shift;

##### get the longest transcript of each mRNA #####

if( $in =~ /.*\.gz$/ ){
	open IN, "gzip -dc $in | " || die "$!\n";
}else{
	open IN, $in || die "$!\n";
}

my ($ref, $gene, $trans);
while( <IN> ){
	chomp;
	next if /#/;
	my @tmp = split(/\t/, $_);
	if( $tmp[2] eq "gene" ){
		$gene = $1 if $tmp[8] =~ /ID=(.*);Name.*/;
	}
	if( $tmp[2] eq "mRNA" ){
		$trans = $1 if $tmp[8] =~ /ID=.*\|.*\|(.*);Parent.*/;
	}
	if( $tmp[2] eq "CDS" ){
		my ($cds, $parent) = ($1, $2) if $tmp[8] =~ /ID=cds-(.*);Parent=.*\|.*\|(.*);Dbxref=.*/;
		$ref->{$gene}{$cds} += $tmp[4] - $tmp[3] if $parent eq $trans;
	}
}
close IN;

open O1, ">$in.all.mRNA.len.info.xls" || die "$!\n";
open O2, ">$in.all.longst_mRNA.info.xls" || die "$!\n";
for my $key ( sort keys %$ref ){
	my $gene = $ref->{$key};
	print O1 $key;
	print O2 $key;
	for my $mrna ( sort { $gene->{$a} <=> $gene->{$b} } keys %$gene ){
		print O1 "\t$mrna\t$gene->{$mrna}\n";
	}
	my @sort = sort { $gene->{$a} <=> $gene->{$b} } keys %$gene;
	print O2 "\t$sort[-1]\t$gene->{$sort[-1]}\n";
}
close O1;
close O2;

##### get longest mRNAs from the all mRNAs pool #####

my (%Seq, $id);
if( $pool =~ /.*\.*gz$/ ){
	open CDS, "gzip -dc $pool | " || die "$!\n";
}else{
	open CDS, $pool || die "$!\n";
}

while( <CDS> ){
	chomp;
	if( />.*_cds_(.*)_\d+\s+\[\S+.*/ ){
		$id = $1;
	}else{
		$Seq{$id} .= $_;
	}
}
close CDS;

open IN, "$in.all.longst_mRNA.info.xls" || die "$!\n";
open OUT, ">$in.longest_transcript.cds.fa" || die "$!\n";
while( <IN> ){
	chomp;
	my @tmp = split(/\s+/, $_);
	print OUT ">$tmp[1]\n";
	for(my $i = 0; $i <= length($Seq{$tmp[1]}); $i += 80){
		print OUT substr($Seq{$tmp[1]}, $i, 80), "\n";
	}
}
close IN;
close OUT;
