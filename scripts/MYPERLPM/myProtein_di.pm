#!/usr/bin/env perl

# utilities related to handling proteins, amino acids

package myProtein;
use strict;
use warnings;
use Exporter;

our @ISA= qw( Exporter );
our @EXPORT = qw(_mapAAletters);

my %AA = (ALA=>'A',TYR=>'Y',MET=>'M',LEU=>'L',CYS=>'C',GLY=>'G',
         ARG=>'R',ASN=>'N',ASP=>'D',GLN=>'Q',GLU=>'E',HIS=>'H',TRP=>'W',
         LYS=>'K',PHE=>'F',PRO=>'P',SER=>'S',THR=>'T',ILE=>'I',VAL=>'V');
 
my %aa = reverse %AA;

sub _mapAAletters {
  
  my $input = shift;
  my $len = length($input);
   ( $len == 1 || $len == 3 ) || die "unrecognized aa code: $input\n";

  if ( $len == 1 ) { 
    (exists $aa{$input} ) || die "unrecognized aa code: $input\n";
    return $aa{$input};
  }
  else {
    (exists $AA{$input} ) || die "unrecognized aa code: $input\n";
    return $AA{$input};
  }
}

1;
