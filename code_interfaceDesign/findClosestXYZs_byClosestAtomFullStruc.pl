#!/usr/bin/env perl

# script by enghui yap - Jan 2013
# takes a full xyz file, a full structure file (reformated: noH,CYS/HIS), 
# and a interface atom list (full and interface index must match!)
# for each xyz, find the nearest heavy atom to it
# if this atom is part of the interface list, print the xyz mesh
# also print out the nearest interface atom to each xyz
# 
# !!!! IMPORTANT !!!!
# the PDB file must be generated from amber, such that the atom numbers match that used in 
# simulation

use strict;

# main program
(@ARGV==3) || die("Usage: [fullxyzfile] [full_reformated pdbfile][interface atomlist] (full & int pdb files in amber index!)\n");
 
my $fullxyzfile = $ARGV[0];
my $fullpdbfile = $ARGV[1];
my $interfaceatomlist = $ARGV[2];

#============================================================================
# load the fullpdb file
my (@atomN,@residues, @residueN, @atomnames, @x,@y,@z,@r,@atomTag);
my $n=0;
open(PDBFILE, $fullpdbfile) or die ("Can't open $fullpdbfile\n");
while(my $line=<PDBFILE>){
  my @buf = split(" ",$line);
  if($buf[0] !~/^ATOM/) {next; } 
  $atomN[$n]    = substr($line,6,5);
  $atomnames[$n]= substr($line,12,4);
  $residues[$n] = substr($line,17,3);
  $residueN[$n] = substr($line,22,4);
  $x[$n] = substr($line,30,8); 
  $y[$n] = substr($line,38,8); 
  $z[$n] = substr($line,46,8); 
  $r[$n] = substr($line,62,7);
  
  (my $residno = $residueN[$n] ) =~ s/^\s+|\s+$//g; 
  (my $atom = $atomnames[$n]) =~ s/^\s+|\s+$//g;  
  $atomTag[$n]=$residno.".".$atom;
  #print "$x[$n] $y[$n] $z[$n] $atomTag[$n]\n";
  $n++;
}
my $natoms=$n;
close PDBFILE;
#============================================================================
# load the int file (residue,residno, atom ==> for atomTag only)
my %hashIntAtomTag;

open(INTFILE, $interfaceatomlist) or die ("Can't open $interfaceatomlist\n");
while(my $line=<INTFILE>){
  my ($residue,$residno,$atom) = split(" ",$line);
  my $key = $residno.".".$atom;
  #printf "X%sX\n",$key;
  ( ! exists $hashIntAtomTag{$key}) || die "error: interface atomTag $key already exist";
  $hashIntAtomTag{$key} = 1; 
}
close INTFILE;
#============================================================================
# load the xyz points
open(XYZFILE, $fullxyzfile) or die ("Can't open $fullxyzfile\n");
my (@xm,@ym,@zm);
my $m = 0;
while (<XYZFILE>) {
  ($xm[$m],$ym[$m],$zm[$m]) = split(" ");
  $m++;
}
my $nxyzpoint=$m;
close XYZFILE;
#============================================================================
# for each xyz point, find the closest atom in the full structure 

for(my $m=0;$m<$nxyzpoint;$m++) {
  
  my $mindist = 1000;
  my $minatomid = -1;
  for(my $n=0;$n<$natoms;$n++) {
    my $distsq = ($xm[$m]-$x[$n])**2 + ($ym[$m]-$y[$n])**2 + ($zm[$m]-$z[$n])**2;    
    if ($distsq < $mindist) { $mindist = $distsq; $minatomid = $n; }
  }
  my $closestAtomTag = $atomTag[$minatomid];

  if ( exists $hashIntAtomTag{$closestAtomTag} ) {
      printf "%8.3f %8.3f %8.3f %8.3f -> %4d %-4s %3s %4d %6.4f\n",$xm[$m],$ym[$m],$zm[$m],sqrt($mindist),$atomN[$minatomid],$atomnames[$minatomid], $residues[$minatomid],$residueN[$minatomid],$r[$minatomid];
  }
#  else { printf "X%sX %8.3f %8.3f %8.3f %8.3f -> %4d %-4s %3s %4d %6.4f\n", $closestAtomTag, $xm[$m],$ym[$m],$zm[$m],sqrt($mindist),$atomN[$minatomid],$atomnames[$minatomid], $residues[$minatomid],$residueN[$minatomid],$r[$minatomid]; }

}
