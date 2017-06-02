#!/usr/bin/perl
use strict; 
use warnings;

## Written by Jolien

my $fasta_directory = $ARGV[0];
my $fasta_output_directory = $ARGV[1];
my $fasta_output_directory_kclust_files = "$fasta_output_directory/kClust_files";
my $score_per_column=$ARGV[2];

mkdir($fasta_output_directory) or die "Couldn't create $fasta_output_directory directory\n";
mkdir($fasta_output_directory_kclust_files) or die "Couldn't create $fasta_output_directory_kclust_files directory\n";

opendir(FASTA, $fasta_directory) or die "Cannot open fasta input directory\n";
while(readdir FASTA){
	my $file=$_;
	if ($file=~ m/\.fa$/){
		my ($cluster) = $file =~ m /^(.+)\.fa$/;
		system("kClust -i $fasta_directory/$file -d $fasta_output_directory -s $score_per_column");
		rename("$fasta_output_directory/representatives.fas", "$fasta_output_directory/$file");
		rename("$fasta_output_directory/db_sorted.fas", "$fasta_output_directory_kclust_files/$cluster.db_sorted.fas");
		rename("$fasta_output_directory/clusters.dmp", "$fasta_output_directory_kclust_files/$cluster.clusters.dmp");
		rename("$fasta_output_directory/headers.dmp", "$fasta_output_directory_kclust_files/$cluster.headers.dmp");
	}
}
closedir(FASTA);
exit;

