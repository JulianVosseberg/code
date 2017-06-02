#!/usr/bin/perl
use strict; 
use warnings;

## Original script by Jolien, updated for Pfam by Julian
## This file reads in the domain tabular output format from HMMsearch and based on this assigns asgard protein sequences to Pfam fams
## Per HMMsearch, it saves all hits of the Pfam HMM (general parsing)
## For the hit sequences, it identifies the best Pfam hit, and furthermore a list of additional Pfams by which it is hit 
## The output is a list of the Pfams and their protein sequence contents, referred to by their IDs


##files
my $hmm_output=$ARGV[0];
my $Pfam_protein_contents=$hmm_output.".parsed.asgard_memberships";


##variables
my $hit_length_cutoff=40;		##min hit length	


##select all Pfams significantly hit by a sequence, save these hits
##for each of the hit sequences, save the Pfams by which it is hit, including the i-value.
my %pfam_all_hits;
my %asgard_sequences; 
open(HMM, "<".$hmm_output) or die "Cannot open hmmsearch infile\n";
while(<HMM>){
	my $line=$_;
	chomp $line;
	if ($line=~ m/^#/){
		next;
	}
	my @data=split(/\s+/, $line);
	my $sequence_ID = $data[0];
	my ($pfam) = $data[4] =~ m/^(PF\d{5})\./;
	my $i_value = $data[12];
	my $hmm_start = $data[15];
	my $hmm_stop = $data[16];
	my $seq_start =$data[17]; 
	my $seq_stop = $data[18];
	if (($hmm_stop - $hmm_start) < $hit_length_cutoff){	## only consider hits spanning at least X positions in the hmm profile
		next;
	}
	if (exists $pfam_all_hits{$pfam}{$sequence_ID}){
		if ($i_value < $pfam_all_hits{$pfam}{$sequence_ID}){
			$pfam_all_hits{$pfam}{$sequence_ID}=[$i_value, $seq_start, $seq_stop];
			$asgard_sequences{$sequence_ID}{$pfam}=[$i_value, $seq_start, $seq_stop];
		}
		else{
			next;
		}
	}
	else{
		$pfam_all_hits{$pfam}{$sequence_ID}=[$i_value, $seq_start, $seq_stop];
		$asgard_sequences{$sequence_ID}{$pfam}=[$i_value, $seq_start, $seq_stop];
	}
}
close(HMM);


##find the contents of each Pfam based on the best matches of the sequences
my %pfam_contents;
my %sequence_membership;
foreach my $seq (keys %asgard_sequences){
	my %pfam_hit_scores=();
	foreach my $pfam_hit (keys %{$asgard_sequences{$seq}}){
		my $i_value=@{$asgard_sequences{$seq}{$pfam_hit}}[0];
		$pfam_hit_scores{$pfam_hit}=$i_value;
	}
	my @pfams_ordered=sort {$pfam_hit_scores{$a} <=> $pfam_hit_scores{$b}} keys(%pfam_hit_scores);
	my $best_hit_pfam=$pfams_ordered[0];
	$pfam_contents{$best_hit_pfam}{$seq}=$asgard_sequences{$seq}{$best_hit_pfam};
	$sequence_membership{$seq}=$best_hit_pfam;
}


##print the Pfam contents
open (OUT, ">".$Pfam_protein_contents) or die "Cannot open Pfam protein contents outfile\n";
foreach my $pfam (keys %pfam_contents){
	print OUT $pfam;
	foreach my $member (keys %{$pfam_contents{$pfam}}){
		print OUT "\t$member";
	}
	print OUT "\n";
}
close (OUT);

exit;

