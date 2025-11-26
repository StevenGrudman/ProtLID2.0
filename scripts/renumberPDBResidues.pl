#!/usr/bin/perl -w
use strict;

die "\nUSAGE: renumberPDBResidues [option] <pdb>\n
  [option]
        -i: Residues independently numbered within chains
        -n: Residues numbered from 1 to N in the whole file (ignores chains)\n\n"
unless scalar@ARGV == 2;
my $option = $ARGV[0];
my $file   = $ARGV[1];

# Open pdb file
open(IN, $file ) or die "ERROR: Could not open file $file\n";
my @data = <IN>;
close IN;

# Renumber residues on specified pdb file
( my $pdb = $file ) =~ s/.*\///g;
$pdb =~ s/tmp_//;
open( OUT, ">tmp_$pdb" ) or die "ERROR: Could not create file tmp_$pdb\n";
my( $presSeq, $pchainID );
my $i = 0;
my $j = 1;
foreach my $line ( @data ) {
	next unless $line =~ /^ATOM|^HETATM/;

	my $name    = substr($line, 13, 3);
	my $resSeq  = substr($line, 22, 4);
	my $chainID = substr($line, 21, 1);
	# Initialize values
	if( $i == 0 ) {
		$presSeq  = $resSeq;
		$pchainID = $chainID;
	}
	$j = sprintf( "%4u", $j );

	# Renumber residues independently for each chain in the pdb (start with 1 at the begining of each chain)
	if( $option eq '-i' ) {
		if( $pchainID eq $chainID ) { # same chainID
			if( $presSeq ne $resSeq ) { # resSeq transition
				$j++;
				$j = sprintf( "%4u", $j );
				substr($line, 22, 4, $j);
				print OUT $line;
			}
			else { # same resSeq
				substr($line, 22, 4, $j);
				print OUT $line;
			}
		}
		else { # chainID transition (implies resSeq transition)
			$j = 1; # start over
			$j = sprintf( "%4u", $j );
			substr($line, 22, 4, $j );
			print OUT $line;
		}
	}
	# Renumber residues from 1 to N in the whole PDB file (chains ignored)
	elsif( $option eq '-n' ) {
		if( $presSeq ne $resSeq ) { # resSeq transition
			$j++;
			$j = sprintf( "%4u", $j );
			substr($line, 22, 4, $j);
			print OUT $line;
		}
		else { # same resSeq
			substr($line, 22, 4, $j );
			print OUT $line;
		}
		
	}
	else {
		print STDERR "ERROR: [option] should be either -i or -n\n";
	}
	$pchainID = $chainID;
	$presSeq  = $resSeq;
	$i++;
}
close OUT;
print STDERR "INFO: Residues renumbered in $pdb\n";

exit;
