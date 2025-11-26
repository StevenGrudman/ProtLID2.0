#!/usr/bin/env perl

# script by enghui yap
# part of PROLID suite of programs
# - generate pir file for MODELLER to perform probe placement
# - assume there is a file (tleap.rec.in) in the current directory 
#   containing the receptor sequence (in 3-letter code)
# - input argument: probe 3-letter sequence, e.g. ALA-ARG-GLU

use strict;
use lib $ENV{'PROTLIDv1HOME'}."/scripts/MYPERLPM";
use myProtein;

# main program
(@ARGV==1) || die("Usage: [probe's aaTag e.g. ALA-ARG-GLU]\n") ;


my @probeSeq_threeLetterArray = split("-",$ARGV[0]);
my $probeSeq;
foreach my $aa3 (@probeSeq_threeLetterArray) { $probeSeq .= _mapAAletters($aa3);}


# load tleap.rec.in to get one-letter receptor sequence
open(FILE, "tleap.rec.in") || die "cannot open tleap.rec.in\n";
my @oneLetterRecSeqArray;
my $bStart;
while (<FILE>) {
    chomp;
    if (/recseq/) { $bStart = 1; next;}
    elsif (/}/) { last; }
    if ($bStart) {
	if (/CYX/) { $_ = "CYS";}
	elsif (/HIE/) { $_ = "HIS";}
	push(@oneLetterRecSeqArray,_mapAAletters($_));
    }
}
close FILE;
my $oneLetterRecSeq=join("",@oneLetterRecSeqArray);

# writing out PIR file
 
## sequence = receptor+probe
print ">P1;PROBE_PLUS_RECEPTOR\n";
print "sequence:none::.::.::::\n";
print $oneLetterRecSeq."/".$probeSeq."*\n";

## structure1 = receptor + blankprobe
(my $blankProbeSeq = $probeSeq) =~ s/\w/-/g;
print ">P1;REC\n";
print "structureX:rec.reformat.pdb:FIRST: :LAST: ::::\n";
print $oneLetterRecSeq."/".$blankProbeSeq."*\n";

## structure2 = blankreceptor + centerProbeSeq (center res: 1 for odd probe,2 for even probe, blank otherwise)
(my $blankRecSeq = $oneLetterRecSeq) =~ s/\w/-/g;

my $centerProbeSeq = $blankProbeSeq;
my $probeLen = length($probeSeq); 
if ( $probeLen % 2 == 0 ) { 
    my $i = $probeLen/2;

    substr($centerProbeSeq,$i-1,1) = substr($probeSeq,$i-1,1);
    substr($centerProbeSeq,$i  ,1) = substr($probeSeq,$i  ,1);
}
else {
    my $i = ($probeLen-1)/2;
    substr($centerProbeSeq,$i,1) = substr($probeSeq,$i,1);
}
print ">P1;CENTER_RES\n";
print "structureX:temp_centerRes_full.pdb:FIRST: :LAST: ::::\n";
print $blankRecSeq."/".$centerProbeSeq."*\n";
exit;
