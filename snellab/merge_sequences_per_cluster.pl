#!/usr/bin/perl
use strict; 
use warnings;
use File::Slurp;

## Script by Jolien

my $clusterfile=$ARGV[0];
my $prokaryotes_dir=$ARGV[1];
my $eukaryotes_dir=$ARGV[2];
my $asgards_dir=$ARGV[3];
my $merged_dir=$ARGV[4];

my %clusters;
open(CLUSTERS, "<".$clusterfile) or die "Cannot open cluster file\n";
while(<CLUSTERS>){
	my $line=$_;
	chomp $line;
	my $cluster=(split("\t", $line))[0];
	$clusters{$cluster}="";
}
close(CLUSTERS);


opendir(PROK, $prokaryotes_dir) or die "Cannot open prokaryotic directory\n";
while(readdir PROK){
	my $file=$_;
	if ($file=~ m/^C\d+/){
		my $text = read_file("$prokaryotes_dir/$file");
		my ($cluster) = $file =~ m/^(C\d+)\./;
		$clusters{$cluster}.=$text;
	}
}
closedir(PROK);


opendir(EUK, $eukaryotes_dir) or die "Cannot open prokaryotic directory\n";
while(readdir EUK){
	my $file=$_;
	if ($file=~ m/^C\d+/){
		my $text = read_file("$eukaryotes_dir/$file");
		my ($cluster) = $file =~ m/^(C\d+)\./;
		$clusters{$cluster}.=$text;
	}
}
closedir(EUK);


opendir(ASG, $asgards_dir) or die "Cannot open prokaryotic directory\n";
while(readdir ASG){
	my $file=$_;
	if ($file=~ m/^C\d+/){
		my $text = read_file("$asgards_dir/$file");
		my ($cluster) = $file =~ m/^(C\d+)\./;
		$clusters{$cluster}.=$text;
	}
}
closedir(ASG);


mkdir($merged_dir) or die "Couldn't create $merged_dir directory\n";
foreach my $cluster (keys %clusters){
	my $fasta = $cluster.".fa";
	my $contents = $clusters{$cluster};
	if ($contents eq ""){
		print "$cluster has no sequence contents\n";
	}
	else{
		open (FASTA, ">$merged_dir/$fasta") or die "Cannot open $fasta\n";
		print FASTA $contents;
		close(FASTA);
	}
}
exit;
