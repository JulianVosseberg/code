#!/usr/bin/perl
use strict; 
use warnings;

## Script by Jolien
## This script searches for the prokaryotic member sequences of the EggNOG (COG/ENOG) families that are within a cluster;
## As inputs, it uses the file containing the (numbered) clusters and the prokaryotic part of EggNOG containing all sequences and their memberships. 
## As output, it generates a tab separated file containing the cluster identifier followed by the identifiers of the EggNOG sequences, and a directory containing fasta-formatted files with the sequences
##themselves, named by the cluster name. 

my $cluster_file=$ARGV[0];
my $prokaryotic_sequences=$ARGV[1]; 
my $cluster_contents="BBH_BH_clusters.prokaryotic_EggNOG_sequences.tab";
my $cluster_fasta_directory="BBH_BH_clusters.prokaryotic_EggNOG_sequences_fastas";

##parse clusters, select COG/ENOGs per cluster
my %clusters;
my %cluster_memberships;
open (CLUSTERS, "<".$cluster_file) or die "Cannot open cluster file\n";
while(<CLUSTERS>){
	my $line=$_;
	chomp $line;
	my @data = split("\t", $line);
	my $cluster=$data[0];
	foreach my $fam(@data[1..$#data]){
		if ($fam =~ m/^COG|^ENOG/){
			push (@{$clusters{$cluster}}, $fam);
			$cluster_memberships{$fam}=$cluster;
		}
	}
}
close(CLUSTERS);

my %cluster_IDs;
my %cluster_sequences;
my $current_sequence_ID="";
open (ALL, "<".$prokaryotic_sequences) or die "Cannot open sequence file\n";
while(<ALL>){
	my $line=$_;
	chomp $line;
	if ($line=~ m/>/){
		$current_sequence_ID="";			##empty, only substance if the sequence is a member of a COG/ENOG that is in a cluster
		my ($membership)= $line =~ m/^>\d+\..+\.(.*)$/;
		my ($sequence_ID)= $line =~ m/^>(\d+\..+)\..*$/; ##check format in which ID can best be saved. What about multiple transcripts?
		if ($membership eq ""){
			next;
		}
		my @memberships=($membership); ## if the sequence is assigned to multiple COGs, it is saved as part of the individual COGs, and of the combined. 
		if ($membership=~ m/_/){
			my @memberships_split=split("_", $membership);
			push (@memberships, @memberships_split);
		}
		foreach my $fam (@memberships){		
			if (exists $cluster_memberships{$fam}){
				my $cluster=$cluster_memberships{$fam};
				push (@{$cluster_IDs{$cluster}}, $sequence_ID);
				$current_sequence_ID=$sequence_ID;
			}	
		}
	}
	elsif ($line ne ""){
		if ($current_sequence_ID ne ""){
			$cluster_sequences{$current_sequence_ID}.=$line."\n";
		}
	}
}
close(ALL);

open(OUTFILE, ">".$cluster_contents) or die "Cannot open file with the sequence ID contents of the clusters\n";
mkdir($cluster_fasta_directory) or die "Couldn't create $cluster_fasta_directory directory\n";
foreach my $cluster (keys %cluster_IDs){
	my @sequence_IDs=@{$cluster_IDs{$cluster}};
	print OUTFILE $cluster."\t".join("\t",@sequence_IDs)."\n";
	my $fasta = $cluster.".fa";
	open (FASTA, ">$cluster_fasta_directory/$fasta") or die "Cannot open $fasta\n";
	foreach my $sequence (@sequence_IDs){
		if ($cluster_sequences{$sequence} ne ""){
			print FASTA ">".$sequence."\n".$cluster_sequences{$sequence};		
		}
		else{
			print "No sequence found for $sequence\n";
		}
	}
	close(FASTA);
}
close(OUTFILE);
exit;
	
	
		
