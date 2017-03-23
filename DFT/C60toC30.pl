#!/usr/bin/perl

use strict;
use warnings;

my ($inp_filename) = @ARGV;
open(INP, $inp_filename) or die "\nCould not open $inp_filename !!\n";


# Now retrieve lines from .xyz file and print in .bas format
my @x;
my @y;
my @z;
my $med_x = 0;
my $med_y = 0;
my $med_z = 0;
my $counter = 0;
foreach my $line (<INP>){
  if ($line =~m/\s\s6\s+(-?\d.\d+)\s+(-?\d.\d+)\s+(-?\d.\d+)/){
    # print "$1\n"; # $1 contains the y coordinate of an atom
    push(@x,$1);
    $med_x = $med_x + $1;
    push(@y,$2);
    $med_y = $med_y + $2;
    push(@z,$3);
    $med_z = $med_z + $3;
    $counter = $counter + 1;
  }
}

close INP;

$med_x = $med_x / $counter;
$med_y = $med_y / $counter;
$med_z = $med_z / $counter;
# print "-----------------------------------\n";
print "The mean value of x coord is $med_x\n";
print "The mean value of y coord is $med_y\n";
print "The mean value of z coord is $med_z\n";



my ($out_filename) = "C30_generated.bas";
open(OUT, ">", $out_filename) or die "\nCould not open $out_filename !!\n";

# Print header of .bas file
print OUT "30\n\n";
# print "30\n\n";

my $c;
open(INP, $inp_filename) or die "\nCould not open $inp_filename !!\n";
foreach my $line (<INP>){
  if ($line =~m/\s\s6\s+(-?\d.\d+)\s+(-?\d.\d+)\s+(-?\d.\d+)/){
    if ($3 > $med_z ){
      $c = 1;
    } else {
      print OUT "$line";
    }
  }
}

close INP;
close OUT;
