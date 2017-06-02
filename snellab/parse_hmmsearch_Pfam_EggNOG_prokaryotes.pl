#!/usr/bin/perl
use strict; 
use warnings;

## Original script written by Jolien, updated for Pfam by Julian
## This file reads in the domain tabular output format from HMMsearch
## Per HMMsearch, it saves all hits of the Pfam HMM


my $hmm_output=$ARGV[0];
my $Pfam_EggNOG_hits=$hmm_output."_parsed_Pfam";
my $EggNOG_Pfam_hits=$hmm_output."_parsed_EggNOG";


## variables
my $hit_number_cutoff=5; ##minimun
my $hit_length_cutoff=40; ## minimum
my $overlap_cutoff=40; ## minimum


##for all pfams, collect the hits with COG/ENOG membership;
##COG, save the e-value and the region of the Pfam HMM by which it is hit; cog_data=array(counter, best_sequence, i-value, Pfam_start, Pfam_stop)
my %pfams;
my %cogs;
my $current_pfam="";
my %current_cogs=();
my %current_hit_sequences=();
my $flag=0; 
open(HMM, "<".$hmm_output) or die "Cannot open hmmsearch infile\n";
while(<HMM>){
	my $line=$_;
	chomp $line;
	if ($line=~ m/^#/){
		next;
	}
	my @data=split(/\s+/, $line);
	my $hit = $data[0];
	my ($cog_hit)=$hit=~ m/\.(COG.+|ENOG.+)$/;
	my ($pfam) = $data[4] =~ m/^(PF\d{5})\./;
	my $i_value = $data[12];
	my $hmm_start = $data[15];
	my $hmm_stop = $data[16];
	if ($current_pfam ne $pfam){	
		$current_pfam=$pfam;
		%current_cogs=();
		%current_hit_sequences=();
		$pfams{$current_pfam}="";		## will be replaced if COGs are found
	}
	unless ($cog_hit){		## only consider hits with COG/ENOG membership	
		next;
	}
	if (exists $current_hit_sequences{$hit}){
		if ($i_value > $current_hit_sequences{$hit}){
			next;
		}
		else{
			$flag=1;
			$current_hit_sequences{$hit}=$i_value;
		}
	}
	else{
		$current_hit_sequences{$hit}=$i_value;
		$flag=0;
	}
	if (exists $current_cogs{$cog_hit}){
		my @cog_data=@{$current_cogs{$cog_hit}};
		my $counter=$cog_data[0];
		if ($flag == 0){				## only add up to the number of hits if the sequence was not included yet by other domains that were significantly hit
			$counter++;
		} 			
		if ($i_value < $cog_data[2]){			## check if domain e-value is lower
			@cog_data=($counter,$hit,$i_value,$hmm_start,$hmm_stop);
			$current_cogs{$cog_hit}=[@cog_data];
			$cogs{$cog_hit}{$current_pfam}=[@cog_data];
		}
		else{
			shift(@cog_data);unshift(@cog_data,$counter);
			$current_cogs{$cog_hit}=[@cog_data];
			$cogs{$cog_hit}{$current_pfam}=[@cog_data];
		}
	}
	else{
		my $counter=1;		
		my @cog_data=($counter,$hit,$i_value,$hmm_start,$hmm_stop);
		$current_cogs{$cog_hit}=[@cog_data];
		$cogs{$cog_hit}{$current_pfam}=[@cog_data];
	}	
	$pfams{$current_pfam}={%current_cogs};		##update the hashes of the current pfam
}
close(HMM);


##Filter COG/ENOGs that are hit rarely ($hit_number_cutoff)
##Filter COG/ENOGs that are hit only by a small region of the Pfam profile ($hit_length_cutoff)
##Filter potential COG/ENOG fusions (while the Pfam only matches one of them). 
##If one of the other cog hits equals one of the fusion elements, remove this fusion;
##If none of the fused COGs/ENOGs occurs without the other, it is (in combination) kept as a COG/ENOG hit to the Pfam.
##Since the hit is completely removed if one (or more) of the COGs/ENOGs occur without the other, it is possible that the highest hit of the independent COG is not the factual highest (because this was the fused one. The fused COG hit might however span multiple regions (e.g. due to backtransfer), due to which this high score does not represent the Pfam - (single) COG/ENOG relation.
foreach my $pfam_family (keys %pfams){
	if ($pfams{$pfam_family} eq ""){
		next;
	}
	my %cog_hits = %{$pfams{$pfam_family}};
	foreach my $cog_hit(keys %cog_hits){
		my @info = @{$pfams{$pfam_family}{$cog_hit}};
		my $hit_number=$info[0];
		my $hit_length=$info[4] - $info[3];
		if (($hit_number < $hit_number_cutoff) or ($hit_length < $hit_length_cutoff)){			 
			delete $pfams{$pfam_family}{$cog_hit};			
			delete $cogs{$cog_hit};
		}
	}
	foreach my $cog_hit (keys %cog_hits){
		if ($cog_hit=~ m/_/){
			my @cogs_fused=split('_', $cog_hit);
			foreach my $cog_fused(@cogs_fused){
				if (exists $pfams{$pfam_family}{$cog_fused}){
					delete $pfams{$pfam_family}{$cog_hit};			
					delete $cogs{$cog_hit};
					last;
				}
			}
		}
	}
	if (!%{$pfams{$pfam_family}}){
		$pfams{$pfam_family}="";
	}
}



##For each Pfam family, iterate over the hits
##Remove COGs that have only one member sequence hit by the Pfam 
##Check if the different hits hit the Pfam in the same or in different regions. If different, split Pfam.  	

my %new_pfams;
my %cogs_new_assigned;
my %split_pfams;
foreach my $pfam_family (keys %pfams){
	if ($pfams{$pfam_family} eq ""){
		print "$pfam_family has no significant hit to a sequence that is a member of a COG/ENOG\n";
		next;
	}
	my %cog_hit_length;	## sort on the longest region
	foreach my $cog_hit (keys %{$pfams{$pfam_family}}){
		my @info = @{$pfams{$pfam_family}{$cog_hit}};
		my @range=($info[3]..$info[4]);
		$cog_hit_length{$cog_hit}=scalar(@range);
	}		
	my %current_regions=();
	my %current_subs=();
	my $sub_counter=0;
	foreach my $cog_hit (sort {$cog_hit_length{$b} <=> $cog_hit_length{$a}} keys(%cog_hit_length)){
		my @info = @{$pfams{$pfam_family}{$cog_hit}};
		my @range=($info[3]..$info[4]);
		my $new_domain_flag=1;
	 	if (%current_regions){
			foreach my $cog (keys %current_regions){
				my $overlap=overlap(\@{$current_regions{$cog}}, \@range);
				if ($overlap >= $overlap_cutoff){
					$current_regions{$cog_hit}=[@range];
					my $affiliated_sub_counter=$current_subs{$cog};
					$current_subs{$cog_hit}=$affiliated_sub_counter;
					$new_domain_flag=0;
					last;
				}
			}
			if ($new_domain_flag == 1){
				$current_regions{$cog_hit}=[@range];
				$sub_counter++;
				$current_subs{$cog_hit}=$sub_counter;
			}
		}
		else{
			$current_regions{$cog_hit}=[@range];
			$sub_counter++;
			$current_subs{$cog_hit}=$sub_counter;
		}
	}
	if ($sub_counter > 1){
		print "$pfam_family is regionally split into different families\n";
		my $pfam_split=0;
		my @split_pfam_names=();
		while ($sub_counter > 0){
			$pfam_split++;
			my $pfam_name=$pfam_family."_".$pfam_split;
			foreach my $cog_assign (keys %current_subs){
				if ($sub_counter == $current_subs{$cog_assign}){
					$new_pfams{$pfam_name}{$cog_assign}=[@{$pfams{$pfam_family}{$cog_assign}}];
					$cogs_new_assigned{$cog_assign}{$pfam_name}=[@{$pfams{$pfam_family}{$cog_assign}}];
				}
			}
			push(@split_pfam_names, $pfam_name);
			$sub_counter--;
		}
		$split_pfams{$pfam_family}=[@split_pfam_names];
	}
	elsif ($sub_counter == 1){		
		foreach my $cog_assign (keys %current_subs){
			$new_pfams{$pfam_family}{$cog_assign}=[@{$pfams{$pfam_family}{$cog_assign}}];
			$cogs_new_assigned{$cog_assign}{$pfam_family}=[@{$pfams{$pfam_family}{$cog_assign}}];
		}
	}
	else{
		print "$pfam_family cannot be assigned to a COG/ENOG\n";
	}
}


##iterate over new_pfams, print best hit (based on lowest i-value) and other (overlapping) hits
##print the broadest region hit by this set of COGs. 
open(OUT1, ">".$Pfam_EggNOG_hits) or die "Cannot open outfile Pfam-EggNOG\n";
foreach my $pfam (keys %new_pfams){
	print OUT1 $pfam;
	my %best_COG=();
	my $start="";
	my $stop="";
	foreach my $COG (keys %{$new_pfams{$pfam}}){
		my @COG_data=@{$new_pfams{$pfam}{$COG}};
		my $COG_i_value=$COG_data[2];
		my $COG_start=$COG_data[3];
		my $COG_stop=$COG_data[4];
		if ($start eq ""){
			$start=$COG_start;
		}
		elsif($COG_start < $start){
			$start=$COG_start;
		}
		if ($stop eq ""){
			$stop=$COG_stop;
		}
		elsif($COG_stop > $stop){
			$stop=$COG_stop;
		}	
		$best_COG{$COG}=$COG_i_value;
	}
	print OUT1 "\t$start\t$stop";
	foreach my $sorted_COG (sort {$best_COG{$a} <=> $best_COG{$b}} keys(%best_COG)){
		print OUT1 "\t$sorted_COG";
	}
	print OUT1 "\n";
}
close(OUT1);

##iterate over the COGs, and print to which Pfams they are assigned
##print in which region they hit each Pfam's profile
open(OUT2, ">".$EggNOG_Pfam_hits) or die "Cannot open outfile EggNOG-Pfam\n";
foreach my $collected_cog (keys %cogs_new_assigned){
	print OUT2 $collected_cog;
	my %best_pfam=();
	foreach my $pfam (keys %{$cogs_new_assigned{$collected_cog}}){
		my @pfam_data=@{$cogs_new_assigned{$collected_cog}{$pfam}};
		my $i_value=$pfam_data[2];
		$best_pfam{$pfam}=$i_value;
	}
	foreach my $sorted_pfam (sort {$best_pfam{$a} <=> $best_pfam{$b}} keys(%best_pfam)){
		my @info = @{$cogs_new_assigned{$collected_cog}{$sorted_pfam}};
		print OUT2 "\t".$sorted_pfam." (".$info[3]."..".$info[4].")";
	}
	print OUT2 "\n";	
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
	

