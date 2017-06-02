#!/usr/bin/perl
use strict; 
use warnings;

##Original script by Jolien, updated for Pfam by Julian
##This script identifies BBH pairs consisting of Pfam - COG/ENOG, based on the hmmsearch of Pfam vs. prokaryotic EggNOG. 
##It adds to the BBH clusters the Pfam and EggNOG family that have, as their best hit, a family in a BBH pair. 
##It uses the parsed output files. 

my $parsed_Pfam=$ARGV[0];
my $parsed_EggNOG=$ARGV[1];
my $clusterlist="BBH_BH_clusters.".$parsed_Pfam;

my $overlap_cutoff=40; ## minimal number of overlapping positions in the Pfam profile

my %combined_families; ##contains the currently unassigned families

##For each Pfam, find the best EggNOG hit.
my %Pfam_best_hit;
open (IN1, "<".$parsed_Pfam) or die "Cannot open parsed Pfam file\n";
while(<IN1>){
	my $line=$_;
	chomp $line;
	my @data=split "\t", $line;
	my $Pfam=$data[0];
	if ($data[3]){
		my $best_hit=$data[3];
		$Pfam_best_hit{$Pfam}=$best_hit;	
		$combined_families{$Pfam}="";
	}
}
close(IN1);

##For each EggNOG, find the best Pfam hit.
##Save the region in which it hits the Pfam
my %EggNOG_best_hit;
my %EggNOG_best_hit_region;
my %EggNOG_hits;
open (IN2, "<".$parsed_EggNOG) or die "Cannot open parsed Pfam file\n";
while(<IN2>){
	my $line=$_;
	chomp $line;
	my @data=split "\t", $line;
	my $EggNOG=$data[0];
	if ($data[1]){
		$data[1] =~ m/^(\w+)\s\((\d+..\d+)\)$/;
		my $best_hit = $1; my $region=$2;	
		$EggNOG_best_hit{$EggNOG}=$best_hit;			##also include region
		$EggNOG_best_hit_region{$EggNOG}{$best_hit}=$region;
		$combined_families{$EggNOG}="";
	}
	if (scalar @data > 2){
		foreach my $other (@data[2..$#data]){
			$other =~ m/^(\w+)\s\((\d+..\d+)\)$/;
			my $hit = $1; my $region=$2;				
			$EggNOG_hits{$EggNOG}{$hit}=$region;			
		}
	} 
}
close(IN2);


##Find BBH pairs
my %clusters; 
my %cluster_assignment_Pfam;
my %cluster_assignment_Pfam_region;
my %cluster_assignment_EggNOG;
my $cluster_counter=0;
foreach my $Pfam (keys %Pfam_best_hit){
	my $best_hit=$Pfam_best_hit{$Pfam};		## best_hit is EggNOG fam
	if ($EggNOG_best_hit{$best_hit} eq $Pfam){	##check
		$cluster_counter++;
		my $cluster_name="C".$cluster_counter;	
		$clusters{$cluster_name}=[$Pfam,$best_hit];
		my $region=$EggNOG_best_hit_region{$best_hit}{$Pfam};
		$cluster_assignment_Pfam{$Pfam}=$cluster_name;
		$cluster_assignment_Pfam_region{$Pfam}{$cluster_name}=$region;
		$cluster_assignment_EggNOG{$best_hit}=$cluster_name;
		delete $combined_families{$Pfam};
		delete $combined_families{$best_hit};
	}
}


##Iterate and assign to cluster containing the best hit
##Region check: only if a EggNOG family hits a Pfam in the same region as through which it came in a cluster, it is added. 
##Otherwise the EggNOG is removed, and the chain is broken
my $number_unassigned=scalar(keys %combined_families);
my %removed_EggNOG;
my %removed_Pfam;
while ($number_unassigned > 0){
	foreach my $fam (keys %combined_families){
		if ($fam=~ m/^PF/){
			my $best_hit=$Pfam_best_hit{$fam};
			if (exists $cluster_assignment_EggNOG{$best_hit}){
				my $cluster = $cluster_assignment_EggNOG{$best_hit};
				my $region =$EggNOG_hits{$best_hit}{$fam};
				push (@{$clusters{$cluster}}, $fam);
				$cluster_assignment_Pfam{$fam}=$cluster;
				$cluster_assignment_Pfam_region{$fam}{$cluster}=$region;
				delete $combined_families{$fam};
			}
			if (exists $removed_EggNOG{$best_hit}){
				delete $combined_families{$fam};
				$removed_Pfam{$fam}="";
			}
		}
		elsif ($fam=~ m/^COG|^ENOG/){
			my $best_hit=$EggNOG_best_hit{$fam};	
			my @span1=split (/\.\./, $EggNOG_best_hit_region{$fam}{$best_hit});
			my @best_hit_region=($span1[0]..$span1[1]);
			if (exists $cluster_assignment_Pfam{$best_hit}){
				my $cluster = $cluster_assignment_Pfam{$best_hit};
				my @span2=split (/\.\./,$cluster_assignment_Pfam_region{$best_hit}{$cluster});
				my @cluster_region=($span2[0]..$span2[1]);
				my $overlap=overlap(\@best_hit_region, \@cluster_region);
				if ($overlap >= $overlap_cutoff){
					push (@{$clusters{$cluster}}, $fam);
					$cluster_assignment_EggNOG{$fam}=$cluster;
					delete $combined_families{$fam};
				}
				else{
					delete $combined_families{$fam}; 
					$removed_EggNOG{$fam}="";
				}			
			}
			if (exists $removed_Pfam{$best_hit}){
				delete $combined_families{$fam};
				$removed_EggNOG{$fam}="";
			}
		}
	}
	$number_unassigned=scalar(keys %combined_families);
}
	
	
open(OUT, ">".$clusterlist) or die "Cannot open $clusterlist\n";
foreach my $group (keys %clusters){
	my @group_members=@{$clusters{$group}};
	print OUT $group;
	foreach my $member (@group_members) {
		if ($member =~ m/^PF/){
			print OUT "\t".$member." (".$cluster_assignment_Pfam_region{$member}{$group}.")";
		}
		else{
			print OUT "\t".$member;
		}
	}
	print OUT "\n";
}
close(OUT);


sub overlap {
	my ($region1, $region2) = @_;
	my $overlap_counter=0;
	foreach my $position1 (@$region1){
		foreach my $position2 (@$region2){
			if ($position1 == $position2){
				$overlap_counter++;
			}
		}		
	} 
	return $overlap_counter;
}

exit;

