#!/usr/bin/perl
use strict; 
use warnings;

## Original script by Jolien, updated for Eukarya4 dataset by Julian
## This script searches for the eukaryotic sequences of the PANTHER families that are within a cluster;
## As inputs, it uses the file containing the (numbered) clusters, the file containing the sequence members of each PANTHER, and the fasta file containing eukaryotic EggNOG sequences. 
## As output, it generates a tab separated file containing the cluster identifier followed by the identifiers of the sequences, and a directory containing fasta-formatted files with the sequences
##themselves, named by the cluster name. 

my $cluster_file=$ARGV[0];
my $panther_families=$ARGV[1];
my $eukaryotic_sequences=$ARGV[2]; 
my $cluster_contents="BBH_BH_clusters.eukaryotic_PANTHER_sequences.tab";
my $cluster_fasta_directory="BBH_BH_clusters.eukaryotic_PANTHER_sequences_fastas";


##parse clusters, select PANTHERs per cluster
my %clusters;
my %cluster_memberships;
open (CLUSTERS, "<".$cluster_file) or die "Cannot open cluster file\n";
while(<CLUSTERS>){
	my $line=$_;
	chomp $line;
	my @data = split("\t", $line);
	my $cluster=$data[0];
	foreach my $fam(@data[1..$#data]){
		if ($fam =~ m/^PTHR/){ 			
			my ($fam_only)= $fam =~ m/^(PTHR\d{5})_*/;
			push (@{$clusters{$cluster}}, $fam_only);
			$cluster_memberships{$fam_only}=$cluster;
		}
	}
}
close(CLUSTERS);


##select sequences for the panthers, match them to its respective cluster
##panthers should all be in a cluster, otherwise warning
my %proteins_in_cluster;
open (FAMILY, "<".$panther_families) or die "Cannot open PANTHER families file\n";
while(<FAMILY>){
	my $line=$_;
	chomp $line;
	my @data = split("\t", $line);
	my $panther=$data[0];
	if (exists $cluster_memberships{$panther}){
		my $cluster = $cluster_memberships{$panther};
		foreach my $protein (@data[1..$#data]){
			$proteins_in_cluster{$protein}=$cluster; ##protein can only be in a single PANTHER
		}
	}
	else{
		print "$panther is not in a cluster\n";
	}
}
close(FAMILY);


#Save sequences that are in a cluster
my %cluster_sequences;
my %cluster_IDs;
my $current_sequence_ID="";
open(ALL, "<".$eukaryotic_sequences) or die "Cannot open full fasta file of all eukaryotes\n";
while(<ALL>){
	my $line=$_;
	chomp $line;
	if ($line=~ m/>/){
		$line =~ m/^>(.+)$/;
		my $sequence_ID = $1; 
		if (exists $proteins_in_cluster{$sequence_ID}){
			$cluster_sequences{$sequence_ID}="";
			my $cluster=$proteins_in_cluster{$sequence_ID};
			push (@{$cluster_IDs{$cluster}}, $sequence_ID);
			$current_sequence_ID=$sequence_ID;
		}
		else{
			$current_sequence_ID="";
		}
	}
	elsif ($line ne ""){
		if ($current_sequence_ID ne ""){
			$cluster_sequences{$current_sequence_ID}.=$line."\n";
		}
	}
}
close(ALL);

##
open(OUTFILE, ">".$cluster_contents) or die "Cannot open file with the sequence ID contents of the clusters\n";
mkdir($cluster_fasta_directory) or die "Couldn't create $cluster_fasta_directory directory\n";
foreach my $cluster (keys %cluster_IDs){
	my @sequence_IDs=@{$cluster_IDs{$cluster}};
	print OUTFILE $cluster."\t".join("\t",@sequence_IDs)."\n";
	my $fasta = $cluster.".fa";
	open (FASTA, ">$cluster_fasta_directory/$fasta") or die "Cannot open $fasta\n";
	foreach my $sequence (@sequence_IDs){
		print FASTA ">".$sequence."\n".$cluster_sequences{$sequence};
	}
	close(FASTA);
}
close(OUTFILE);
exit;	



