#!/usr/bin/perl
use strict; 
use warnings;

## Original script by Jolien, updated for Pfam by Julian
## This file reads in the domain tabular output format from HMMsearch
## The parsing procedure consists of two parts: assigning protein sequences to Pfam fams (part 1), and potential merging of Pfam fams into clusters (part 2);
## Per HMMsearch, it saves all significants hit of the Pfam HMM (general parsing)
## For the hit sequences, it identifies the best Pfam hit, and furthermore a list of additional Pfams by which it is hit (part 1)
## The output is a list of the Pfams and their protein sequence contents, referred to by their IDs (part 1)
## If Pfams share hits, and if they are not in separate clusters in the Pfam-EggNOG merge, they will be merged (part 2)
## The output is a new list of clusters, and the protein contents of the individual Pfam fams. (part 2)


##files
my $hmm_output=$ARGV[0];
my $Pfam_EggNOG_clusters=$ARGV[1];
my $Pfam_protein_contents=$hmm_output.".parsed.eukaryotic_memberships";
my $Pfam_EggNOG_clusters_update=$Pfam_EggNOG_clusters.".eukaryotic_update";


##variables
my $hit_length_cutoff=40;		##min hit length	
my $overlap_cutoff=40;
my $shared_sequences_cutoff=5;		##min number of hits that Pfams should share before considering merging. 


##select all Pfams hit by a sequence, save these hits
##for each of the hit sequences, save the Pfams by which it is hit, including the e-value.
##general
my %pfam_all_hits;
my %eukaryotic_sequences; ##saving only the ID itself, and not the memberships, will suffice.
open(HMM, "<".$hmm_output) or die "Cannot open hmmsearch infile\n";
while(<HMM>){
	my $line=$_;
	chomp $line;
	if ($line=~ m/^#/){
		next;
	}
	my @data=split(/\s+/, $line);
	my $hit = $data[0];
	my ($sequence_ID)=$hit; ##only one transcript per sequence, so no need to change something
	my ($pfam) = $data[4] =~ m/^(PF\d{5})\./;
	my $i_value = $data[12];
	my $hmm_start = $data[15];
	my $hmm_stop = $data[16];
	my $seq_start =$data[17]; 
	my $seq_stop = $data[18];
	if (($hmm_stop - $hmm_start) < $hit_length_cutoff){	## only consider hits spanning at least 30 positions
		next;
	}
	if (exists $pfam_all_hits{$pfam}{$sequence_ID}){
		if ($i_value < $pfam_all_hits{$pfam}{$sequence_ID}){
			$pfam_all_hits{$pfam}{$sequence_ID}=[$i_value, $hmm_start, $hmm_stop, $seq_start, $seq_stop];
			$eukaryotic_sequences{$sequence_ID}{$pfam}=[$i_value, $hmm_start, $hmm_stop, $seq_start, $seq_stop];
		}
		else{
			next;
		}
	}
	else{
		$pfam_all_hits{$pfam}{$sequence_ID}=[$i_value, $hmm_start, $hmm_stop, $seq_start, $seq_stop];
		$eukaryotic_sequences{$sequence_ID}{$pfam}=[$i_value, $hmm_start, $hmm_stop, $seq_start, $seq_stop];
	}
}
close(HMM);


##find the contents of each Pfam based on the best matches of the sequences
##find also the affiliated sequences of each Pfams (sequences that are significantly hit by the Pfam, but that do not have that Pfam as their best match)
##part 1
my %pfam_contents;
my %pfam_affiliates;
my %sequence_membership;
foreach my $seq (keys %eukaryotic_sequences){
	my %pfam_hit_scores=();
	foreach my $pfam_hit (keys %{$eukaryotic_sequences{$seq}}){
		my $i_value=@{$eukaryotic_sequences{$seq}{$pfam_hit}}[0];
		$pfam_hit_scores{$pfam_hit}=$i_value;
	}
	my @pfams_ordered=sort {$pfam_hit_scores{$a} <=> $pfam_hit_scores{$b}} keys(%pfam_hit_scores);
	my $best_hit_pfam=$pfams_ordered[0];
	$pfam_contents{$best_hit_pfam}{$seq}=$eukaryotic_sequences{$seq}{$best_hit_pfam};
	$sequence_membership{$seq}=$best_hit_pfam;
	if (scalar (keys %pfam_hit_scores) > 1){
		foreach my $other (@pfams_ordered[1..$#pfams_ordered]){
			$pfam_affiliates{$other}{$seq}=$eukaryotic_sequences{$seq}{$other};	##save sequences as 'affiliated sequence' to the Pfams it is hit by, but not the best
		}
	}
}


##print the Pfam contents
##part 1
open (OUT1, ">".$Pfam_protein_contents) or die "Cannot open Pfam protein contents outfile\n";
foreach my $pfam (keys %pfam_contents){
	print OUT1 $pfam;
	foreach my $member (keys %{$pfam_contents{$pfam}}){
		print OUT1 "\t$member";
	}
	print OUT1 "\n";
}
close (OUT1);


##save clusters to check if Pfams with shared hits are in one and/or the same;
##part 2
my %pfam_clusters;
my %pfam_clusters_split;
my %pfam_clusters_regions;
my %cluster_contents;
my $cluster_counter=0;
open (CLUSTERS, "<".$Pfam_EggNOG_clusters) or die "Cannot open cluster file\n";
while(<CLUSTERS>){
	my $line=$_;
	chomp $line;
	my @data=split ("\t", $line);
	my $cluster=$data[0];
	$cluster_counter++;
	my @families=@data[1..$#data];
	foreach my $fam (@families){
		if ($fam=~ m/^PF/){
			$fam =~ m/^(PF\S+)\s\((\d+)\.\.(\d+)\)$/;
			my $fam_name=$1; my $start=$2; my $stop=$3;
			if ($fam =~ m/^PF\d+_\d/){
				my ($pfam_only) = $fam_name =~ m/^(PF\d+)_\d/;
				push (@{$pfam_clusters_split{$pfam_only}}, $fam_name); 
				$pfam_clusters{$pfam_only}="";
			}
			$pfam_clusters{$fam_name}=$cluster;
			$pfam_clusters_regions{$fam_name}{$cluster}=[($start..$stop)];
			push(@{$cluster_contents{$cluster}}, $fam_name);
		}
		else{
			push(@{$cluster_contents{$cluster}}, $fam);
		}
	}
}
close(CLUSTERS);


#find the Pfam fams that have contents, but are not in a cluster
my %pfams_unclustered;
foreach my $pfam (keys %pfam_contents){
	if (exists $pfam_clusters{$pfam}){
		next;
	}
	else{
		$pfams_unclustered{$pfam}="";
	}
}


##For each Pfam, check if it hits sequences from another Pfam; 
##Save the best hit 'other Pfam' if they hit sequences in the same region and if they share at least ... sequence hits
##Save the region by which the best hit 'other Pfam' hits the shared best sequence. 
##needed to save best sequence for each affiliated Pfam. 
##part 2
my %best_links;
my %best_links_hit_region;
my %linked_hit_region;
my %pfams_unlinked;		##unclustered pfams that cannot be linked (either because they have no affl sequences or because they have too few/no overlap);
foreach my $pfam_unclustered (keys %pfams_unclustered){
	if (exists $pfam_affiliates{$pfam_unclustered}){
		my %affiliated_pfam_shared_seq_number=();
		my %affiliated_pfam_sum_i_value=(); 	
		my %affiliated_pfam_seq_i_value=();	## sequence with the lowest i-value, used for determining which region is hit in the affiliated Pfam. 
		foreach my $affiliated_sequence (keys %{$pfam_affiliates{$pfam_unclustered}}){
			my $affiliated_pfam=$sequence_membership{$affiliated_sequence};
			my @region1=($pfam_affiliates{$pfam_unclustered}{$affiliated_sequence}[3]..$pfam_affiliates{$pfam_unclustered}{$affiliated_sequence}[4]);			
			my @region2=($pfam_contents{$affiliated_pfam}{$affiliated_sequence}[3]..$pfam_contents{$affiliated_pfam}{$affiliated_sequence}[4]);
			my $overlap = overlap(\@region1, \@region2);
			if ($overlap >= $hit_length_cutoff){
				if (exists $affiliated_pfam_shared_seq_number{$affiliated_pfam}){
					$affiliated_pfam_shared_seq_number{$affiliated_pfam}++;
					my $i_value=$pfam_affiliates{$pfam_unclustered}{$affiliated_sequence}[0];
					if ($affiliated_pfam_shared_seq_number{$affiliated_pfam} <= $shared_sequences_cutoff){
						$affiliated_pfam_sum_i_value{$affiliated_pfam}+=$i_value;
						$affiliated_pfam_seq_i_value{$affiliated_pfam}{$affiliated_sequence}=$i_value;
					}
				}
				else{
					$affiliated_pfam_shared_seq_number{$affiliated_pfam}=1;
					my $i_value=$pfam_affiliates{$pfam_unclustered}{$affiliated_sequence}[0];
					$affiliated_pfam_sum_i_value{$affiliated_pfam}=$i_value;
					$affiliated_pfam_seq_i_value{$affiliated_pfam}{$affiliated_sequence}=$i_value;
				}
			}
		}
		if (%affiliated_pfam_shared_seq_number){ 			## not empty
			foreach my $affiliated_pfam (keys %affiliated_pfam_shared_seq_number){
				if ($affiliated_pfam_shared_seq_number{$affiliated_pfam} < $shared_sequences_cutoff){
					delete $affiliated_pfam_shared_seq_number{$affiliated_pfam};
					delete $affiliated_pfam_sum_i_value{$affiliated_pfam};
					delete $affiliated_pfam_seq_i_value{$affiliated_pfam};
				}
			}
			if (scalar (keys %affiliated_pfam_shared_seq_number) > 0){		
				my @pfam_links=sort {$affiliated_pfam_sum_i_value{$a} <=> $affiliated_pfam_sum_i_value{$b}} keys (%affiliated_pfam_sum_i_value);
				my $best_link = $pfam_links[0];
				$best_links{$pfam_unclustered}=$best_link;
				my @seq = sort {$affiliated_pfam_seq_i_value{$best_link}{$a}<=>$affiliated_pfam_seq_i_value{$best_link}{$b}} keys(%{$affiliated_pfam_seq_i_value{$best_link}});
				my $best_seq = $seq[0];
				my @region_best_link = ($pfam_contents{$best_link}{$best_seq}[1]..$pfam_contents{$best_link}{$best_seq}[2]);
				$best_links_hit_region{$pfam_unclustered}{$best_link}=[@region_best_link];
				my @region_linked = ($pfam_affiliates{$pfam_unclustered}{$best_seq}[1]..$pfam_affiliates{$pfam_unclustered}{$best_seq}[2]);
				$linked_hit_region{$pfam_unclustered}=[@region_linked];
			}
			else{
				$pfams_unlinked{$pfam_unclustered}=""; ## although there are affiliated sequences, there are no affiliated Pfams found due to lack of shared sequences
			}
		}
		else{
			$pfams_unlinked{$pfam_unclustered}="";	## although there are affiliated sequences, there are no affiliated Pfams found (likely due to no overlap)
		}
	}
	else{
		$pfams_unlinked{$pfam_unclustered}=""; ## no affiliated sequences, so no affiliated Pfams. 
	}
}


## assign unclustered, but linked Pfams to existing clusters, based on best match;
## check if regions overlap 
## if best match is a split Pfam, compare the different overlapping regions in the best match, assign cluster to that of split pfam with highest overlap
## do this iteratively
## part 2
my %updated_clusters=%cluster_contents;
my %updated_pfam_clusters=%pfam_clusters;
my %updated_pfam_clusters_regions=%pfam_clusters_regions;
my $best_links_in_cluster=0;
foreach my $unassigned (keys %best_links){
	my $best_link=$best_links{$unassigned};
	if (exists $pfam_clusters{$best_link}){
		$best_links_in_cluster++;
	}
}
while($best_links_in_cluster > 0){
	foreach my $unassigned (keys %best_links){
		my $best_link = $best_links{$unassigned};
		if (exists $updated_pfam_clusters{$best_link}){
			if (exists $pfam_clusters_split{$best_link}){
				my @region_in_link = @{$best_links_hit_region{$unassigned}{$best_link}};
				my @split_pfams=@{$pfam_clusters_split{$best_link}};
				my $best_overlap=0;
				my $best_split="";
				foreach my $split_pfam (@split_pfams){
					my $cluster2=$updated_pfam_clusters{$split_pfam};
					my @region_in_cluster = @{$updated_pfam_clusters_regions{$split_pfam}{$cluster2}};
					my $overlap = overlap(\@region_in_cluster, \@region_in_link);
					if ($overlap > $best_overlap){
						$best_overlap = $overlap;
						$best_split=$split_pfam;
					}
				}
				if ($best_overlap >= $overlap_cutoff){
					my $cluster2=$updated_pfam_clusters{$best_split};
					push (@{$updated_clusters{$cluster2}}, $unassigned);
					$updated_pfam_clusters{$unassigned}=$cluster2;
					$updated_pfam_clusters_regions{$unassigned}{$cluster2}=[@{$linked_hit_region{$unassigned}}];
					print "$unassigned and $best_split share hits and $unassigned is added to $cluster2\n";	
				}
				else{
					$pfams_unlinked{$unassigned}="";
				}
				delete $best_links{$unassigned};
			}
			else{		
				my $cluster2=$updated_pfam_clusters{$best_link};
				my @region_in_cluster = @{$updated_pfam_clusters_regions{$best_link}{$cluster2}};
				my @region_in_link = @{$best_links_hit_region{$unassigned}{$best_link}};
				my $overlap = overlap(\@region_in_cluster, \@region_in_link);
				if ($overlap >= $overlap_cutoff){
					push (@{$updated_clusters{$cluster2}}, $unassigned);
					$updated_pfam_clusters{$unassigned}=$cluster2;
					$updated_pfam_clusters_regions{$unassigned}{$cluster2}=[@{$linked_hit_region{$unassigned}}];
					print "$unassigned and $best_link share hits and $unassigned is added to $cluster2\n";	
				}
				else{
					$pfams_unlinked{$unassigned}="";
				}
				delete $best_links{$unassigned};	## best link is used for merging OR best link is not qualified for merging due to limited region overlap
			}
		}
	}
	$best_links_in_cluster=0;	
	foreach my $unassigned (keys %best_links){
		my $best_link=$best_links{$unassigned};
		if (exists $updated_pfam_clusters{$best_link}){
			$best_links_in_cluster++;
		}
	}
}


##'BBH' pairs of still unassigned --> new clusters
foreach my $leftover1 (keys %best_links){
	my $best_link=$best_links{$leftover1};
	unless($best_link){
		next;
	}
	if (exists $best_links{$best_link}){
		if ($best_links{$best_link} eq $leftover1){	
			$cluster_counter++;
			my $new_cluster="C".$cluster_counter;
			$updated_clusters{$new_cluster}=[$leftover1, $best_link]; 
			$updated_pfam_clusters{$leftover1}=$new_cluster; 
			$updated_pfam_clusters{$best_link}=$new_cluster;
			$updated_pfam_clusters_regions{$leftover1}{$new_cluster}=[@{$linked_hit_region{$leftover1}}];
			$updated_pfam_clusters_regions{$best_link}{$new_cluster}=[@{$linked_hit_region{$best_link}}];
			print "$leftover1 and $best_link share hits and together form a new, Pfam only cluster: $new_cluster\n";
			delete $best_links{$leftover1};delete $best_links{$best_link};
		}
	}
}	


## Iterate, add to BBH pairs
my $removed=scalar(keys(%best_links));		
while ($removed > 0){	
	my $number_leftover2 = scalar(keys(%best_links));
	foreach my $leftover2 (keys %best_links){
		my $best_link = $best_links{$leftover2};
		if (exists $updated_pfam_clusters{$best_link}){
			my $add_cluster=$updated_pfam_clusters{$best_link};
			my @best_hit_region=@{$best_links_hit_region{$leftover2}{$best_link}};
			my @cluster_region=@{$updated_pfam_clusters_regions{$best_link}{$add_cluster}};
			my $overlap =overlap(\@best_hit_region, \@cluster_region);
			if ($overlap >= $overlap_cutoff){
				push(@{$updated_clusters{$add_cluster}}, $leftover2);
				$updated_pfam_clusters{$leftover2}=$add_cluster;
				$updated_pfam_clusters_regions{$leftover2}{$add_cluster}=[@{$linked_hit_region{$leftover2}}];
				print "$leftover2 is added to the Pfam only cluster $add_cluster\n";	
			}	
			else{
				$pfams_unlinked{$leftover2}="";
			}
			delete $best_links{$leftover2};			## best link is used for merging OR best link is not qualified for merging due to limited region overlap
		}
	}
	$removed = $number_leftover2 -  scalar(keys(%best_links));
}


##if a Pfam is not yet in a cluster and also cannot be assigned to one, let it be a cluster on its own;
##2 situations: has a best match, but could nevertheless not be assigned to a (existing/new) cluster OR has no best match
##part 2
foreach my $pfam_leftover (keys %best_links){
	$cluster_counter++;
	my $cluster_new="C".$cluster_counter;
	$updated_clusters{$cluster_new}=[$pfam_leftover];
	$updated_pfam_clusters{$pfam_leftover}=$cluster_new;
	print "$pfam_leftover could not be assigned to a cluster and now forms a cluster on its own\n";
}
foreach my $pfam_unlinked (keys %pfams_unlinked){
	$cluster_counter++;
	my $cluster_new="C".$cluster_counter;
	$updated_clusters{$cluster_new}=[$pfam_unlinked];
	$updated_pfam_clusters{$pfam_unlinked}=$cluster_new;
	print "$pfam_unlinked could not be assigned to a cluster and now forms a cluster on its own\n";
}


##print the new clusters
##part 2
open (OUT2, ">".$Pfam_EggNOG_clusters_update) or die "Cannot open the update of the cluster file\n";
foreach my $cluster_final (keys %updated_clusters){
	print OUT2 "$cluster_final\t".join("\t", @{$updated_clusters{$cluster_final}})."\n";
}
close(OUT2);


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
