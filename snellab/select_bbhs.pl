#!/usr/bin/perl

use strict;
use warnings;

# Written by Leny
# This script determines best blast hits for species versus species blast results and uses a species-supergroup list to determine at what level bbh's are determined.
# The output is a list with bbh's.

my $speciesSupergroupList = $ARGV[0];
my $protein_type = $ARGV[1];
my $supergroup_number = $ARGV[2];
my $path = $ARGV[3];

# List with species and their supergroup is parsed

open FH, "<$speciesSupergroupList" or 
   die "Cannot open $speciesSupergroupList\n";




#for my $i (0 ..$#cmd){
#      print $print[$i];
#      shell_command($cmd[$i]);
#    }
 
my %SupergroupperSpecies;

print "#This are the species and their supergroup in $speciesSupergroupList\n";
while( my $line = <FH>){
      chomp($line);
      my @fields = split(/\t/, $line);
      my $species = $fields[0];
      my $supergroup = $fields [1];
      print "$species\t$supergroup\n";
      $SupergroupperSpecies{$species}=$supergroup;
}
      
close FH;

# Blastfiles are parsed, only the ones that are the result of a blast between two supergroups. BBH's between supergroups are collected in %besthit based on highest bitscore.
# If bitscores are exactly the same they are collected in %exacthit. Later on it is checked whether besthits are in %exacthit, if that's the case both entries in %exacthit should be
# taken into account. Otherwise the problem is solved by a hit with a higher bitscore than the ones that have exactly the same bitscore.

my %besthit;
my %exacthit;




open FH_OUT, ">$path/exact.hits" or die
   "Cannot open exact.hits";

print "#These are best hits between supergroups\n";
foreach my $species1 (keys %SupergroupperSpecies) {
   foreach my $species2 (keys %SupergroupperSpecies){ 
      if ($SupergroupperSpecies{$species1} ne $SupergroupperSpecies{$species2}){
         open FH, "<Blast.$protein_type.$species1.$species2.txt";
         while( my $line = <FH> ){
            chomp $line;
            my @fields = split(/\t/, $line);
            my $query = $fields[0];
            my $subject = $fields[1];
            my $bitscore = $fields[11];
            if (not exists $besthit{$SupergroupperSpecies{$species1}}{$SupergroupperSpecies{$species2}}{$query}){
               $besthit{$SupergroupperSpecies{$species1}}{$SupergroupperSpecies{$species2}}{$query}=[$subject,$bitscore];
               print "$SupergroupperSpecies{$species1}\t$SupergroupperSpecies{$species2}\t$query\t$subject\t$bitscore\n";
            } elsif ($bitscore > $besthit{$SupergroupperSpecies{$species1}}{$SupergroupperSpecies{$species2}}{$query}[1]) {
               print "\n$query best hit changed from $besthit{$SupergroupperSpecies{$species1}}{$SupergroupperSpecies{$species2}}{$query}[0] into:\n";
               $besthit{$SupergroupperSpecies{$species1}}{$SupergroupperSpecies{$species2}}{$query}=[$subject,$bitscore];
               print "$SupergroupperSpecies{$species1}\t$SupergroupperSpecies{$species2}\t$query\t$subject\t$bitscore\n\n";
            } elsif ($bitscore == $besthit{$SupergroupperSpecies{$species1}}{$SupergroupperSpecies{$species2}}{$query}[1]) {
               print "\nFor query $query bitscore of $besthit{$SupergroupperSpecies{$species1}}{$SupergroupperSpecies{$species2}}{$query}[0] is exact that of $subject\n\n";                 
               $exacthit{$besthit{$SupergroupperSpecies{$species1}}{$SupergroupperSpecies{$species2}}{$query}[0]}{$query}=$bitscore;
               $exacthit{$query}{$besthit{$SupergroupperSpecies{$species1}}{$SupergroupperSpecies{$species2}}{$query}[0]}=$bitscore;
               $exacthit{$query}{$subject}=$bitscore;
               $exacthit{$subject}{$query}=$bitscore;
               print FH_OUT "$besthit{$SupergroupperSpecies{$species1}}{$SupergroupperSpecies{$species2}}{$query}[0]\t$query\t$bitscore\n";
               print FH_OUT "$subject\t$query\t$bitscore\n";
               }
         }  
         close FH;
      }
   }  
} 

close FH_OUT;

print "#These are the sequences in %exacthit\n";
foreach my $query (keys %exacthit) {
   foreach my $subject (keys %{ $exacthit{$query} }) {
      print "$query\t$subject\t$exacthit{$query}{$subject}\n";
   }
}

# All best hits are printed and bbh's are also printed to bbhlist.

open FH_OUT, ">$path/bbhlist.txt" or
  die "Cannot open bbhlist.txt";

my %bbhs;
my $counter=0;

print "#Best hits per query are:\n";
foreach my $supergroup1 (keys %besthit) {
   foreach my $supergroup2 (keys %besthit){
      foreach my $query (keys  %{ $besthit{$supergroup1}{$supergroup2} }) {
         print "$supergroup1\t$supergroup2\t$query\t$besthit{$supergroup1}{$supergroup2}{$query}[0]\n";
         if ($besthit{$supergroup2}{$supergroup1}{$besthit{$supergroup1}{$supergroup2}{$query}[0]}[0] eq $query){
            print "\nYes! this is a BBH: $query\t$besthit{$supergroup1}{$supergroup2}{$query}[0]\t$besthit{$supergroup1}{$supergroup2}{$query}[1]\n";
            print "The bitscores are $besthit{$supergroup1}{$supergroup2}{$query}[1] and $besthit{$supergroup2}{$supergroup1}{$besthit{$supergroup1}{$supergroup2}{$query}[0]}[1]\n";
            my $highest_bitscore; 
            if ($besthit{$supergroup1}{$supergroup2}{$query}[1] >  $besthit{$supergroup2}{$supergroup1}{$besthit{$supergroup1}{$supergroup2}{$query}[0]}[1]){
               $highest_bitscore = $besthit{$supergroup1}{$supergroup2}{$query}[1];
            } else {
               $highest_bitscore = $besthit{$supergroup2}{$supergroup1}{$besthit{$supergroup1}{$supergroup2}{$query}[0]}[1];
            }
            print "Highest bitscore is $highest_bitscore\n";    
            if (exists $exacthit{$query}{$besthit{$supergroup1}{$supergroup2}{$query}[0]}){
               print "\nWait! for this pair $query\t$besthit{$supergroup1}{$supergroup2}{$query}[0] there is another best hit with exactly the same bitscore which is:\n";
               foreach my $hit (keys %{ $exacthit{$query}} ) {
                  if (($exacthit{$query}{$hit} == $exacthit{$query}{$besthit{$supergroup1}{$supergroup2}{$query}[0]}) and ($besthit{$supergroup1}{$supergroup2}{$query}[0] ne $hit) ){
                     print "$hit\n";
                     if (exists $besthit{$supergroup2}{$supergroup1}{$hit}) {
                        print "Its best hit is $besthit{$supergroup2}{$supergroup1}{$hit}[0]\n";
                        if ($besthit{$supergroup2}{$supergroup1}{$hit}[0] eq $query){
                           print "So bbh is actually triangle\n";
                           print "The other bitscore is $besthit{$supergroup2}{$supergroup1}{$hit}[1]\n";
                           if ($besthit{$supergroup2}{$supergroup1}{$hit}[1] > $highest_bitscore){ 
                              print "bbh should be different\n";
                              print FH_OUT "$query\n";
            		      print FH_OUT "$hit\n";
                              $bbhs{$query}=[$hit,$exacthit{$query}{$hit}];
                              $counter++; 
                            } elsif ($besthit{$supergroup2}{$supergroup1}{$hit}[1] == $highest_bitscore){
                              print "Two sets of bbh's should be added\n";
                               print FH_OUT "$query\n";
                               print FH_OUT "$besthit{$supergroup1}{$supergroup2}{$query}[0]\n";
                               $bbhs{$query}=[$besthit{$supergroup1}{$supergroup2}{$query}];
                               $counter++;
                               print FH_OUT "$query\n";
                               print FH_OUT "$hit\n";
                               $counter++;
                             } else {
                             print FH_OUT "$query\n";
                             print FH_OUT "$besthit{$supergroup1}{$supergroup2}{$query}[0]\n";
                             $bbhs{$query}=[$besthit{$supergroup1}{$supergroup2}{$query}];
                             $counter++;
                          }

                     } else {  
                       print FH_OUT "$query\n";
                       print FH_OUT "$besthit{$supergroup1}{$supergroup2}{$query}[0]\n";
                       $bbhs{$query}=[$besthit{$supergroup1}{$supergroup2}{$query}];
                       $counter++; 
                     }
                    }
                  }
               }
               } else {
               print FH_OUT "$query\n";
               print FH_OUT "$besthit{$supergroup1}{$supergroup2}{$query}[0]\n";
               $bbhs{$query}=[$besthit{$supergroup1}{$supergroup2}{$query}];
               $counter++;
            }
         }
      }
   }
}

print "#The total number of BBHs is $counter\n"; 

# This subroutine handles every command that is executed by the shell. Sprintf, chomp and printing are used to check whether the command is really doing what
# you think it is doing. The result is also printed, chomped and returned for further use.

sub shell_command{
  my $cmd = shift @_;
  $cmd = sprintf($cmd);
  chomp $cmd;
  print "\tShell command: $cmd\n\n";
  my $result = `$cmd`;
  print "\t$result\n";
  print "\n";
  chomp $result;
  return $result;
}

