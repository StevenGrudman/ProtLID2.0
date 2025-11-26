#!/usr/bin/env perl

# script by enghui yap - june 2008
# takes an input pdb file, checks that residue numbers are continuous
# if chain is given, just check that chain
# if there are alternative location use 'A'

(@ARGV>0) || die("Usage: [pdbfile]");
 
$pdbfile = $ARGV[0];
if(@ARGV>1) { $bChain=1; $neededChainID=$ARGV[1];} else {$bChain=0;} 
if(@ARGV>2) { $minRes=$ARGV[2];} else {$minRes=1;}

open(PDBFILE, $pdbfile) or die ("Can't open $pdbfile\n");

# read in the pdb file line by line                                                                  
$prevChainID="0";

while($line=<PDBFILE>){

    @buf = split(" ",$line);
    
    if($buf[0] =~/ATOM/) 
    {

	$altloc =substr($line,16,1);
	$chainid=substr($line, 21,1);
	
	if ($bChain) {
	    if ($chainid !~ $neededChainID ) { 
		if ( !$bChainFound) {next;}
		else {last;}
	    }
	    else {
		$bChainFound = 1;
	    }
	}


	if($altloc=~m/(A| )/)
	{
	    $resno   = substr($line,22,4); 

# check that residues are continuous
	    if ($resno < $minRes ) {next; }
	    
	    if($chainid=~$prevChainID) {
		if($resno-$prevResNo>1) {
		    print "$_\n$pdbfile chain $chainid: Residue Number not continuous: prev: $prevResNo; curr: $resno\n";

		}
	    }
	    else {
		$prevChainID=$chainid;
	    }

	    $prevResNo = $resno;
	}
	   
    }
    
}

