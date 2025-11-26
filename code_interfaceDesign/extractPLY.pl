#!/usr/bin/perl
# written by Enghui Yap - Oct 2012
# extracts vertices from PLY file

(@ARGV == 1 ) || die "usage: plyfile\n"; 
 
$plyfile = $ARGV[0];
#$outfile = "out.xyz";

open(FILE, $plyfile) or die ("Can't open $plyfile\n");
#open(OUT, ">$outfile") or die ("Can't open $outfile\n");

my $bStart = 0; 
my $bEnd = 0; 
$n=0;

while( defined($line=<FILE>) && !$bEnd ){ 

  if( $line =~/^element vertex/ ) { @linebuf = split(" ", $line); $nVert = $linebuf[2];}
  elsif ($line =~ /^end_header/) { $bStart = 1; }
  elsif($bStart) {
    @linebuf = split(" ",$line); 
    my $x =  $linebuf[0];
    my $y =  $linebuf[1];
    my $z =  $linebuf[2];
    print "$x $y $z \n"; 

    $n++;
    if($n == $nVert) {$bEnd = 1; }
  }
}
close FILE;
#close OUT;
