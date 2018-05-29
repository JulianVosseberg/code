#!/usr/bin/perl

use strict;
use warnings;
use List::MoreUtils qw(uniq);

# Written by Leny; modified by Julian (usage and determine first which species are present and checks for empty files)
# This script determines best blast hits (BBHs) for species versus species blast results and uses a species-supergroup list to determine at what level BBHs are determined.
# The output is a list with BBHs.

# State usage
sub usage {
   print STDERR "\nUsage:\n\tselect_bbhs.pl <supergroups.tsv> <species.list> <prefix> <output file>\n";
   exit;
}

# Check number of arguments
if (scalar(@ARGV) != 4){
   usage();
}

# Get arguments
my $speciesSupergroupList = $ARGV[0];
my $speciesList = $ARGV[1];
my $protein_type = $ARGV[2];
my $bbhlist = $ARGV[3]; # output file

# Checks which species are present
open FH, "<$speciesList" or
    die "Cannot open $speciesList\n";

chomp(my @speciesarray = <FH>);
#print "@speciesarray\n";

close FH;

# List with species and their supergroup is parsed

open FH, "<$speciesSupergroupList" or 
   die "Cannot open $speciesSupergroupList\n";

my %SupergroupperSpecies;

#print "#These are the species and their supergroup in $speciesSupergroupList\n";
while( my $line = <FH>){
      chomp($line);
      my @fields = split(/\t/, $line);
      my $species = $fields[0];
      my $supergroup = $fields [1];
#     print "$species\t$supergroup\n";
      if ($species ~~ @speciesarray){
          $SupergroupperSpecies{$species}=$supergroup;
#          print "$species found\n";
      }
}
     
close FH;

# Blastfiles are parsed, only the ones that are the result of a blast between two supergroups. BBHs between supergroups are collected in %besthit based on highest bitscore.
# If bitscores are exactly the same multiple hits are collected. 

my %besthit;

print "#Best hits between supergroups are collected\n\n";
print "Type\tQuery\tHit\tBitscore\tSupergroup_query\tSupergroup_hit\n";
open FH, "<${protein_type}_blastp.txt" or
    die "Cannot open ${protein_type}_blastp.txt\n";
while( my $line = <FH> ){
    chomp $line;
    my @fields = split(/\t/, $line);
    my $query = $fields[0];
    my $species1 = substr($query, 0, 4);
    my $hit = $fields[1];
    my $species2 = substr($hit, 0, 4);
    if ($SupergroupperSpecies{$species1} ne $SupergroupperSpecies{$species2}){
	my $bitscore = $fields[11];
        if (not exists $besthit{$SupergroupperSpecies{$species1}}{$SupergroupperSpecies{$species2}}{$query}){
	    $besthit{$SupergroupperSpecies{$species1}}{$SupergroupperSpecies{$species2}}{$query}=[$bitscore,$hit];
            print "Initial\t$query\t$hit\t$bitscore\t$SupergroupperSpecies{$species1}\t$SupergroupperSpecies{$species2}\n";
         } elsif ($bitscore > $besthit{$SupergroupperSpecies{$species1}}{$SupergroupperSpecies{$species2}}{$query}[0]) {
            $besthit{$SupergroupperSpecies{$species1}}{$SupergroupperSpecies{$species2}}{$query}=[$bitscore,$hit];
            print "Replaced\t$query\t$hit\t$bitscore\t$SupergroupperSpecies{$species1}\t$SupergroupperSpecies{$species2}\n";
         } elsif ($bitscore == $besthit{$SupergroupperSpecies{$species1}}{$SupergroupperSpecies{$species2}}{$query}[0]) {
            push(@{$besthit{$SupergroupperSpecies{$species1}}{$SupergroupperSpecies{$species2}}{$query}},$hit);
            print "Expanded\t$query\t$hit\t$bitscore\t$SupergroupperSpecies{$species1}\t$SupergroupperSpecies{$species2}\n";                             
      }
   }  
}
close FH; 

# Check if not all files were empty
my $size = keys %besthit;
if ($size == 0) {
	print "\n#All supergroup-supergroup BLAST output files are empty\n";
	print STDERR "\nError: All supergroup-supergroup BLAST output files are empty. Script aborted.\n";
	exit;
}

# For queries with multiple best hits with equal bitscore, if these best hits form BBHs with the query, only the hits with the highest reciprocal bitscore are selected.
# This is done by removing the query from the values of the best hit in %besthit and the hit from the values of the query, if the reciprocal bitscore of this best hit is not the maximum reciprocal 
# bitscore. Before removing it is checked whether a query-best hit pair actually forms a BBH, because if only the hits with maximum reciprocal bitscore are kept while a hit with a lower reciprocal 
# bitscore forms an BBH and the tophit does not, BBHs are lost that should have been kept

print "\n#BBH removal if reciprocal bitscores are below maximum:\n";
foreach my $supergroup1 (keys %besthit) {
   foreach my $supergroup2 (keys %besthit){
      my $besthit_size = keys %{ $besthit{$supergroup1}{$supergroup2} };
      # Check if there are any hits between these supergroups
      if ($besthit_size == 0 and $supergroup1 ne $supergroup2){
         print "\n#$supergroup1 versus $supergroup2 BLAST did not result in any significant hits\n";
         print STDERR "\n$supergroup1 versus $supergroup2 BLAST did not result in any significant hits\n";
         next;
      }
      foreach my $query (keys  %{ $besthit{$supergroup1}{$supergroup2} }) {
         my $print=0;
         my @hits_1_2 = @{$besthit{$supergroup1}{$supergroup2}{$query}};
         my $bitscore_query = shift @hits_1_2;
         # Hit names and bitscores are sorted on bitscore, large to small
         my @bitscore_hits;
	 my @hits_1_2_mutual; #will replace @hits_1_2 because @hits_1_2 can be larger than @bitscore_hits and then sorting with @idx removes hits from @hits_1_2 but not necessily the ones that are not in %besthit
         foreach (@hits_1_2) {
           if (exists $besthit{$supergroup2}{$supergroup1}{$_}){ # Otherwise if $_ did not yet exist in $besthit, it is added...
               push @bitscore_hits, @{$besthit{$supergroup2}{$supergroup1}{$_}}[0];               
               push @hits_1_2_mutual, $_;
            }
         }
         my @idx = sort { $bitscore_hits[$b] <=> $bitscore_hits[$a] } 0 .. $#bitscore_hits;
         @bitscore_hits = @bitscore_hits[@idx];
         @hits_1_2_mutual = @hits_1_2_mutual[@idx];
         # Hit names and bitscores are collected if they are part of an BBH with the query
         my @bitscore_hits_bbh;
         my @hits_1_2_bbh;
         my $switch=0;
         foreach my $i (0..$#hits_1_2_mutual) {  
            my @hits_2_1 = @{$besthit{$supergroup2}{$supergroup1}{$hits_1_2_mutual[$i]}};
            shift @hits_2_1;
            if($query ~~ @hits_2_1){
               push @bitscore_hits_bbh, $bitscore_hits[$i];
               push @hits_1_2_bbh, $hits_1_2_mutual[$i]; 
            }                       
         } 
         #Hits with a bitscore lower than the maximum are removed
         foreach my $j (1..$#bitscore_hits_bbh) {   
            if ($bitscore_hits_bbh[$j] < $bitscore_hits_bbh[0]){
               $switch++;
               $print = 1;
               if ($switch==1){
                  print "Query:$query\tBBHs:";
                  foreach my $k (0..$#hits_1_2_bbh){
                     print "$hits_1_2_bbh[$k]:$bitscore_hits_bbh[$k]";
                  }
                  print "\tRemoved:";
               }
               print "$hits_1_2_bbh[$j]\t";
               # Removal occurs by removing the query from the hit and the hit from the query in %besthit 
               @{$besthit{$supergroup2}{$supergroup1}{$hits_1_2_bbh[$j]}} = grep {!/$query/} @{$besthit{$supergroup2}{$supergroup1}{$hits_1_2_bbh[$j]}};                   
               @{$besthit{$supergroup1}{$supergroup2}{$query}} = grep {!/$hits_1_2_bbh[$j]/} @{$besthit{$supergroup1}{$supergroup2}{$query}};
            }
         } 
         if ($print == 1){
            print "\n";
         }
      }
   }
}


#BBHs are collected

my @bbhs;

print "\n#Best hits per query are:\n";
foreach my $supergroup1 (keys %besthit) {
   foreach my $supergroup2 (keys %besthit){
      foreach my $query (keys  %{ $besthit{$supergroup1}{$supergroup2} }) {
         my @hits_1_2 = @{$besthit{$supergroup1}{$supergroup2}{$query}};
         my $bitscore_query = shift @hits_1_2;
         my $hit_total = scalar(@hits_1_2);
         my @bbh_hits;
         print "Query:$query\tHits:@hits_1_2\tHit_total:$hit_total\tBitscore_query:$bitscore_query\tSupergroups:$supergroup1-$supergroup2\t";
         foreach (@hits_1_2) {
            if (exists $besthit{$supergroup2}{$supergroup1}{$_}){ # It is possible that the hit does not exist...
               my @hits_2_1 = @{$besthit{$supergroup2}{$supergroup1}{$_}};
               shift @hits_2_1;
               if ($query ~~ @hits_2_1){
                  print "BBH:$query $_\t"; 
                  print "Bitscore_hit:@{$besthit{$supergroup2}{$supergroup1}{$_}}[0]\t";  
                  push @bbh_hits, $_;                                                   
               }
            }
         }
         my $bbh_total = scalar(@bbh_hits);
         print "BBH_total:$bbh_total\t";                         
         if ($bbh_total > 0){
            push @bbhs, $query;
         }               
         print "\n"; 
      }
   }
}

@bbhs = sort @bbhs;
my $total = scalar(@bbhs);
print "\nThe number of BBHs is $total\n";


# BBH IDs are printed to file

open FH_OUT, ">$bbhlist" or
    die "Cannot open $bbhlist\n";

foreach(@bbhs) {
   print FH_OUT "$_\n";
}

close FH_OUT;
