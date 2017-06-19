#!/usr/bin/perl
use strict; 
use warnings;

## Script by Jolien, updated for Eukarya4 by Julian
## This file reads in the domain tabular output format from HMMsearch
## The parsing procedure consists of two parts: assigning protein sequences to PANTHER fams (part 1), and potential merging of PANTHER fams into clusters (part 2);
## Per HMMsearch, it saves all significants hit of the PANTHER HMM (general parsing)
## For the hit sequences, it identifies the best PANTHER hit, and furthermore a list of additional PANTHERs by which it is hit (part 1)
## The output is a list of the PANTHERs and their protein sequence contents, referred to by their IDs (part 1)
## If PANTHERs share hits, and if they are not in separate clusters in the PANTHER-EggNOG merge, they will be merged (part 2)
## The output is a new list of clusters, and the protein contents of the individual PANTHER fams. (part 2)


##files
my $hmm_output=$ARGV[0];
my $PANTHER_EggNOG_clusters=$ARGV[1];
my $PANTHER_protein_contents=$hmm_output.".parsed.eukaryotic_memberships";
my $PANTHER_EggNOG_clusters_update=$PANTHER_EggNOG_clusters.".eukaryotic_update";


##variables
my $i_value_cutoff=10 ** -6;		##max i_value
my $hit_length_cutoff=40;		##min hit length	
my $overlap_cutoff=40;
my $shared_sequences_cutoff=5;		##min number of hits that PANTHERs should share before considering merging. 


##select all PANTHERs significantly hit by a sequence, save these hits
##for each of the hit sequences, save the PANTHERs by which it is hit, including the e-value.
##general
my %panther_all_hits;
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
	my ($sequence_ID)=$hit;
	my ($panther) = $data[3] =~ m/^(PTHR\d{5})\./;
	my $i_value = $data[12];
	my $hmm_start = $data[15];
	my $hmm_stop = $data[16];
	my $seq_start =$data[17]; 
	my $seq_stop = $data[18];
	if ($i_value > $i_value_cutoff){			## only consider hits with e-value < 1e-5	
		next;
	}
	if (($hmm_stop - $hmm_start) < $hit_length_cutoff){	## only consider hits spanning at least 30 positions
		next;
	}
	if (exists $panther_all_hits{$panther}{$sequence_ID}){
		if ($i_value < $panther_all_hits{$panther}{$sequence_ID}){
			$panther_all_hits{$panther}{$sequence_ID}=[$i_value, $hmm_start, $hmm_stop, $seq_start, $seq_stop];
			$eukaryotic_sequences{$sequence_ID}{$panther}=[$i_value, $hmm_start, $hmm_stop, $seq_start, $seq_stop];
		}
		else{
			next;
		}
	}
	else{
		$panther_all_hits{$panther}{$sequence_ID}=[$i_value, $hmm_start, $hmm_stop, $seq_start, $seq_stop];
		$eukaryotic_sequences{$sequence_ID}{$panther}=[$i_value, $hmm_start, $hmm_stop, $seq_start, $seq_stop];
	}
}
close(HMM);


##find the contents of each PANTHER based on the best matches of the sequences
##find also the affiliated sequences of each PANTHERs (sequences that are significantly hit by the PANTHER, but that do not have that PANTHER as their best match)
##part 1
my %panther_contents;
my %panther_affiliates;
my %sequence_membership;
foreach my $seq (keys %eukaryotic_sequences){
	my %panther_hit_scores=();
	foreach my $panther_hit (keys %{$eukaryotic_sequences{$seq}}){
		my $i_value=@{$eukaryotic_sequences{$seq}{$panther_hit}}[0];
		$panther_hit_scores{$panther_hit}=$i_value;
	}
	my @panthers_ordered=sort {$panther_hit_scores{$a} <=> $panther_hit_scores{$b}} keys(%panther_hit_scores);
	my $best_hit_panther=$panthers_ordered[0];
	$panther_contents{$best_hit_panther}{$seq}=$eukaryotic_sequences{$seq}{$best_hit_panther};
	$sequence_membership{$seq}=$best_hit_panther;
	if (scalar (keys %panther_hit_scores) > 1){
		foreach my $other (@panthers_ordered[1..$#panthers_ordered]){
			$panther_affiliates{$other}{$seq}=$eukaryotic_sequences{$seq}{$other};	##save sequences as 'affiliated sequence' to the PANTHERs it is hit by, but not the best
		}
	}
}


##print the PANTHER contents
##part 1
open (OUT1, ">".$PANTHER_protein_contents) or die "Cannot open PANTHER protein contents outfile\n";
foreach my $panther (keys %panther_contents){
	print OUT1 $panther;
	foreach my $member (keys %{$panther_contents{$panther}}){
		print OUT1 "\t$member";
	}
	print OUT1 "\n";
}
close (OUT1);


##save clusters to check if PANTHERs with shared hits are in one and/or the same;
##part 2
my %panther_clusters;
my %panther_clusters_split;
my %panther_clusters_regions;
my %cluster_contents;
my $cluster_counter=0;
open (CLUSTERS, "<".$PANTHER_EggNOG_clusters) or die "Cannot open cluster file\n";
while(<CLUSTERS>){
	my $line=$_;
	chomp $line;
	my @data=split ("\t", $line);
	my $cluster=$data[0];
	$cluster_counter++;
	my @families=@data[1..$#data];
	foreach my $fam (@families){
		if ($fam=~ m/^PTHR/){
			$fam =~ m/^(PTHR\S+)\s\((\d+)\.\.(\d+)\)$/;
			my $fam_name=$1; my $start=$2; my $stop=$3;
			if ($fam =~ m/^PTHR\d+_\d/){
				my ($panther_only) = $fam_name =~ m/^(PTHR\d+)_\d/;
				push (@{$panther_clusters_split{$panther_only}}, $fam_name); 
				$panther_clusters{$panther_only}="";
			}
			$panther_clusters{$fam_name}=$cluster;
			$panther_clusters_regions{$fam_name}{$cluster}=[($start..$stop)];
			push(@{$cluster_contents{$cluster}}, $fam_name);
		}
		else{
			push(@{$cluster_contents{$cluster}}, $fam);
		}
	}
}
close(CLUSTERS);


#find the PANTHER fams that have contents, but are not in a cluster
my %panthers_unclustered;
foreach my $panther (keys %panther_contents){
	if (exists $panther_clusters{$panther}){
		next;
	}
	else{
		$panthers_unclustered{$panther}="";
	}
}


##For each PANTHER, check if it hits sequences from another PANTHER; 
##Save the best hit 'other PANTHER' if they hit sequences in the same region and if they share at least ... sequence hits
##Save the region by which the best hit 'other PANTHER' hits the shared best sequence. 
##needed to save best sequence for each affiliated PANTHER. 
##part 2
my %best_links;
my %best_links_hit_region;
my %linked_hit_region;
my %panthers_unlinked;		##unclustered panthers that cannot be linked (either because they have no affl sequences or because they have too few/no overlap);
foreach my $panther_unclustered (keys %panthers_unclustered){
	if (exists $panther_affiliates{$panther_unclustered}){
		my %affiliated_panther_shared_seq_number=();
		my %affiliated_panther_sum_i_value=(); 	
		my %affiliated_panther_seq_i_value=();	## sequence with the lowest i-value, used for determining which region is hit in the affiliated PANTHER. 
		foreach my $affiliated_sequence (keys %{$panther_affiliates{$panther_unclustered}}){
			my $affiliated_panther=$sequence_membership{$affiliated_sequence};
			my @region1=($panther_affiliates{$panther_unclustered}{$affiliated_sequence}[3]..$panther_affiliates{$panther_unclustered}{$affiliated_sequence}[4]);			
			my @region2=($panther_contents{$affiliated_panther}{$affiliated_sequence}[3]..$panther_contents{$affiliated_panther}{$affiliated_sequence}[4]);
			my $overlap = overlap(\@region1, \@region2);
			if ($overlap >= $hit_length_cutoff){
				if (exists $affiliated_panther_shared_seq_number{$affiliated_panther}){
					$affiliated_panther_shared_seq_number{$affiliated_panther}++;
					my $i_value=$panther_affiliates{$panther_unclustered}{$affiliated_sequence}[0];
					if ($affiliated_panther_shared_seq_number{$affiliated_panther} <= $shared_sequences_cutoff){
						$affiliated_panther_sum_i_value{$affiliated_panther}+=$i_value;
						$affiliated_panther_seq_i_value{$affiliated_panther}{$affiliated_sequence}=$i_value;
					}
				}
				else{
					$affiliated_panther_shared_seq_number{$affiliated_panther}=1;
					my $i_value=$panther_affiliates{$panther_unclustered}{$affiliated_sequence}[0];
					$affiliated_panther_sum_i_value{$affiliated_panther}=$i_value;
					$affiliated_panther_seq_i_value{$affiliated_panther}{$affiliated_sequence}=$i_value;
				}
			}
		}
		if (%affiliated_panther_shared_seq_number){ 			## not empty
			foreach my $affiliated_panther (keys %affiliated_panther_shared_seq_number){
				if ($affiliated_panther_shared_seq_number{$affiliated_panther} < $shared_sequences_cutoff){
					delete $affiliated_panther_shared_seq_number{$affiliated_panther};
					delete $affiliated_panther_sum_i_value{$affiliated_panther};
					delete $affiliated_panther_seq_i_value{$affiliated_panther};
				}
			}
			if (scalar (keys %affiliated_panther_shared_seq_number) > 0){		
				my @panther_links=sort {$affiliated_panther_sum_i_value{$a} <=> $affiliated_panther_sum_i_value{$b}} keys (%affiliated_panther_sum_i_value);
				my $best_link = $panther_links[0];
				$best_links{$panther_unclustered}=$best_link;
				my @seq = sort {$affiliated_panther_seq_i_value{$best_link}{$a}<=>$affiliated_panther_seq_i_value{$best_link}{$b}} keys(%{$affiliated_panther_seq_i_value{$best_link}});
				my $best_seq = $seq[0];
				my @region_best_link = ($panther_contents{$best_link}{$best_seq}[1]..$panther_contents{$best_link}{$best_seq}[2]);
				$best_links_hit_region{$panther_unclustered}{$best_link}=[@region_best_link];
				my @region_linked = ($panther_affiliates{$panther_unclustered}{$best_seq}[1]..$panther_affiliates{$panther_unclustered}{$best_seq}[2]);
				$linked_hit_region{$panther_unclustered}=[@region_linked];
			}
			else{
				$panthers_unlinked{$panther_unclustered}=""; ## although there are affiliated sequences, there are no affiliated PANTHERs found due to lack of shared sequences
			}
		}
		else{
			$panthers_unlinked{$panther_unclustered}="";	## although there are affiliated sequences, there are no affiliated PANTHERs found (likely due to no overlap)
		}
	}
	else{
		$panthers_unlinked{$panther_unclustered}=""; ## no affiliated sequences, so no affiliated PANTHERs. 
	}
}


## assign unclustered, but linked PANTHERs to existing clusters, based on best match;
## check if regions overlap 
## if best match is a split PANTHER, compare the different overlapping regions in the best match, assign cluster to that of split panther with highest overlap
## do this iteratively
## part 2
my %updated_clusters=%cluster_contents;
my %updated_panther_clusters=%panther_clusters;
my %updated_panther_clusters_regions=%panther_clusters_regions;
my $best_links_in_cluster=0;
foreach my $unassigned (keys %best_links){
	my $best_link=$best_links{$unassigned};
	if (exists $panther_clusters{$best_link}){
		$best_links_in_cluster++;
	}
}
while($best_links_in_cluster > 0){
	foreach my $unassigned (keys %best_links){
		my $best_link = $best_links{$unassigned};
		if (exists $updated_panther_clusters{$best_link}){
			if (exists $panther_clusters_split{$best_link}){
				my @region_in_link = @{$best_links_hit_region{$unassigned}{$best_link}};
				my @split_panthers=@{$panther_clusters_split{$best_link}};
				my $best_overlap=0;
				my $best_split="";
				foreach my $split_panther (@split_panthers){
					my $cluster2=$updated_panther_clusters{$split_panther};
					my @region_in_cluster = @{$updated_panther_clusters_regions{$split_panther}{$cluster2}};
					my $overlap = overlap(\@region_in_cluster, \@region_in_link);
					if ($overlap > $best_overlap){
						$best_overlap = $overlap;
						$best_split=$split_panther;
					}
				}
				if ($best_overlap >= $overlap_cutoff){
					my $cluster2=$updated_panther_clusters{$best_split};
					push (@{$updated_clusters{$cluster2}}, $unassigned);
					$updated_panther_clusters{$unassigned}=$cluster2;
					$updated_panther_clusters_regions{$unassigned}{$cluster2}=[@{$linked_hit_region{$unassigned}}];
					print "$unassigned and $best_split share hits and $unassigned is added to $cluster2\n";	
				}
				else{
					$panthers_unlinked{$unassigned}="";
				}
				delete $best_links{$unassigned};
			}
			else{		
				my $cluster2=$updated_panther_clusters{$best_link};
				my @region_in_cluster = @{$updated_panther_clusters_regions{$best_link}{$cluster2}};
				my @region_in_link = @{$best_links_hit_region{$unassigned}{$best_link}};
				my $overlap = overlap(\@region_in_cluster, \@region_in_link);
				if ($overlap >= $overlap_cutoff){
					push (@{$updated_clusters{$cluster2}}, $unassigned);
					$updated_panther_clusters{$unassigned}=$cluster2;
					$updated_panther_clusters_regions{$unassigned}{$cluster2}=[@{$linked_hit_region{$unassigned}}];
					print "$unassigned and $best_link share hits and $unassigned is added to $cluster2\n";	
				}
				else{
					$panthers_unlinked{$unassigned}="";
				}
				delete $best_links{$unassigned};	## best link is used for merging OR best link is not qualified for merging due to limited region overlap
			}
		}
	}
	$best_links_in_cluster=0;	
	foreach my $unassigned (keys %best_links){
		my $best_link=$best_links{$unassigned};
		if (exists $updated_panther_clusters{$best_link}){
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
			$updated_panther_clusters{$leftover1}=$new_cluster; 
			$updated_panther_clusters{$best_link}=$new_cluster;
			$updated_panther_clusters_regions{$leftover1}{$new_cluster}=[@{$linked_hit_region{$leftover1}}];
			$updated_panther_clusters_regions{$best_link}{$new_cluster}=[@{$linked_hit_region{$best_link}}];
			print "$leftover1 and $best_link share hits and together form a new, PANTHER only cluster: $new_cluster\n";
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
		if (exists $updated_panther_clusters{$best_link}){
			my $add_cluster=$updated_panther_clusters{$best_link};
			my @best_hit_region=@{$best_links_hit_region{$leftover2}{$best_link}};
			my @cluster_region=@{$updated_panther_clusters_regions{$best_link}{$add_cluster}};
			my $overlap =overlap(\@best_hit_region, \@cluster_region);
			if ($overlap >= $overlap_cutoff){
				push(@{$updated_clusters{$add_cluster}}, $leftover2);
				$updated_panther_clusters{$leftover2}=$add_cluster;
				$updated_panther_clusters_regions{$leftover2}{$add_cluster}=[@{$linked_hit_region{$leftover2}}];
				print "$leftover2 is added to the PANTHER only cluster $add_cluster\n";	
			}	
			else{
				$panthers_unlinked{$leftover2}="";
			}
			delete $best_links{$leftover2};			## best link is used for merging OR best link is not qualified for merging due to limited region overlap
		}
	}
	$removed = $number_leftover2 -  scalar(keys(%best_links));
}


##if a PANTHER is not yet in a cluster and also cannot be assigned to one, let it be a cluster on its own;
##2 situations: has a best match, but could nevertheless not be assigned to a (existing/new) cluster OR has no best match
##part 2
foreach my $panther_leftover (keys %best_links){
	$cluster_counter++;
	my $cluster_new="C".$cluster_counter;
	$updated_clusters{$cluster_new}=[$panther_leftover];
	$updated_panther_clusters{$panther_leftover}=$cluster_new;
	print "$panther_leftover could not be assigned to a cluster and now forms a cluster on its own\n";
}
foreach my $panther_unlinked (keys %panthers_unlinked){
	$cluster_counter++;
	my $cluster_new="C".$cluster_counter;
	$updated_clusters{$cluster_new}=[$panther_unlinked];
	$updated_panther_clusters{$panther_unlinked}=$cluster_new;
	print "$panther_unlinked could not be assigned to a cluster and now forms a cluster on its own\n";
}


##print the new clusters
##part 2
open (OUT2, ">".$PANTHER_EggNOG_clusters_update) or die "Cannot open the update of the cluster file\n";
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
