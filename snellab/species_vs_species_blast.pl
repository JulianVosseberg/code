#!/usr/bin/perl

#use strict; (with strict on, the species filehandles cannot be opened)
use warnings;
use List::MoreUtils qw(uniq);

# Original script written by Leny; modified for fasta files containing eukaryotic sequences starting with their abbreviation by Julian
# This script redistributes sequences in a fasta file to per species fasta files. Each species fasta file is converted to a blast database. Then all species versus all species blastp is executed and only the best hit for each search is reported. 
# I have to rewrite the script, the species is list is not really needed, I could get the species from the fastaheader.

my $speciesList = $ARGV[0];
my $fastaFile = $ARGV[1];
my $protein_type = $ARGV[2]; # For example: Kinases
my $cpu_number = $ARGV[3];
my $path = $ARGV[4];

# $speciesList is parsed

open FH, "<$speciesList" or 
   die "Cannot open $speciesList\n";

my $linecounter; 

my @Species;
print STDERR "The species in $speciesList are:\n";
while( my $line = <FH>){
      chomp($line);
      my $species = substr($line, 0, 4);
      print STDERR "$species\n";
      push @Species, $species;
}

close FH;
      
print STDERR "The number of species in $speciesList is: ";
print STDERR scalar(@Species);
print STDERR "\n";

# For each species a filehandle is opened for distribution of fasta's in $fastaFile over different species

foreach my $species (@Species) {
   my $filehandle = "FH.$species";
   open $filehandle, ">$path/$protein_type.$species.fa" or 
   die "Cannot open $path/$protein_type.$species.fa";
} 	

# Sequences in $fastaFile are distributed over new files based on species 

open FH, "<$fastaFile" or
  die "Cannot open $fastaFile\n";

my $species;
my $seq = "";
my $description = "";
my $counter=0;

while( my $line = <FH> ){
    if ($line =~ /^>/){
	# Write previous sequence to file 
        if ( $seq ) { 
           my $filehandle = "FH.$species";
    	   print $filehandle "$description\n$seq";
           $counter++;
        }
        chomp $line;
        # Determine species of sequence and collect fastaheader
        $species = substr($line,4, 4);
        print "$species\n";
        $description = $line;          	
        # Clear variable seq
        $seq = "";
     # Collect squence
     } else {
	$seq .= $line;
    }
}

# Write the last sequence

if ( $seq ) {
   my $filehandle = "FH.$species";
   print $filehandle "$description\n$seq";
   $counter++;
}

# Print statistics
print STDERR "Number of distributed sequences: $counter\n";

# Make Blast databases
foreach my $species (@Species) {
   close "FH.$species";
   my $fasta = "$path/$protein_type.$species.fa";
   my $cmd = "~/.local/makeblastdb -in $fasta -dbtype prot";
   shell_command($cmd);
}

# Run all versus all blast, excluding your own species. Only one target sequence is reported. 
foreach my $species1 (@Species) {
   foreach my $species2 (@Species){ 
      if ($species1 ne $species2) {
      my $cmd= "~/.local/blastp -query $path/$protein_type.$species1.fa -db $path/$protein_type.$species2.fa -max_target_seqs 1 -outfmt 6 -num_threads $cpu_number >$path/Blast.$protein_type.$species1.$species2.txt";
      shell_command($cmd);
      }  
   }
}

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





