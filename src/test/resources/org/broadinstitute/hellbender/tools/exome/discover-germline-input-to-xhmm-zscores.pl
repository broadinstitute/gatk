#!/usr/bin/env perl
use strict;
use warnings;

# script use to transform a read-counts file into the format accepted by XHMM --discover program (basically the transpose and chaning target names to its coordinates as XHMM use them that way)

my $header_line = <STDIN>;

chomp $header_line;
my ($x1, $x2, $x3, $x4, @sample_names) = split(/\t/, $header_line);

my @target_names = ();
my @sample_value_lists = (map { [] } @sample_names); 
while (my $line = <STDIN>) {
  chomp $line;
  my ($chr,$start,$end,$name,@sample_values) = split(/\t/, $line);
  for (my $i = 0; $i < scalar(@sample_names); $i++) {
    push @{$sample_value_lists[$i]}, $sample_values[$i];
  } 
  push @target_names, "$chr:$start-$end";
}

print STDOUT join("\t", "Matrix", @target_names), "\n";
for (my $i = 0; $i < scalar(@sample_names); $i++) {
   print STDOUT join("\t", $sample_names[$i], @{$sample_value_lists[$i]}), "\n";
}
