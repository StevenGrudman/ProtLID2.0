#!/usr/bin/env perl

# script by enghui yap - Jan 2013
# score candidate pdb file against a interface template file
# interface template file has the following format: 
# x y z [aa type] [functional atom; RC=center of benzene ring]
 
use strict;

# main program
(@ARGV==2) || die("Usage: [surface xyz meshpoint file] [pdbfile]");
 
my $xyzfile = $ARGV[0];
my $pdbfile = $ARGV[1];

# load the pqr file
my (@atomN,@residues, @residueN, @atomnames, @x,@y,@z,@r);

my $n=0;
open(PDBFILE, $pdbfile) or die ("Can't open $pdbfile\n");
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

  $n++;
}
my $natoms=$n;

close PDBFILE;

# for each xyz point, find the closest atom
open(XYZFILE, $xyzfile) or die ("Can't open $xyzfile\n");
my $m=1;
while (<XYZFILE>) {
    my ($xm,$ym,$zm) = split(" ");

    my $mindist=1000;
    my $minatomid = -1;
    for(my $n=0;$n<$natoms;$n++) {
      my $dist = &rmsd( $x[$n],$y[$n],$z[$n],$xm,$ym,$zm);
      if ($dist < $mindist) { $mindist = $dist; $minatomid = $n; }
    }
    printf "%8.3f %8.3f %8.3f %8.3f -> %4d %-4s %3s %4d %6.4f\n",$xm,$ym,$zm,sqrt($mindist),$atomN[$minatomid],$atomnames[$minatomid], $residues[$minatomid],$residueN[$minatomid],$r[$minatomid];

    $m++;
}
#============================================================================

sub rmsd {
    my $x1 = $_[0];
    my $y1 = $_[1];
    my $z1 = $_[2];
    my $x2 = $_[3];
    my $y2 = $_[4];
    my $z2 = $_[5];
    return sqrt( ($x1-$x2)*($x1-$x2) +  ($y1-$y2)*($y1-$y2)  + ($z1-$z2)*($z1-$z2)  );
}

