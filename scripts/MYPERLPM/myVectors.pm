#!/usr/bin/perl -w

# vector class to perform basic functions
# 
package myVectors;
use strict;
use warnings;
use Exporter;

our @ISA= qw( Exporter );
our @EXPORT = qw(_norm _normalize _crossProduct _rotateByAxisAngle _multiply_matvec);

sub _norm {
  my $ref_v=shift;
  my $norm = sqrt ( $ref_v->[0] ** 2 + $ref_v->[1] ** 2 + $ref_v->[2] ** 2 );
  return $norm;
}

sub _normalize {
  my $ref_v=shift;
  my $norm = _norm($ref_v);
  $ref_v->[0] /= $norm ;
  $ref_v->[1] /= $norm ;
  $ref_v->[2] /= $norm ;
}
sub _crossProduct {
  
  my $ref_v1=shift;
  my $ref_v2=shift;
  my $ref_crossprod = shift;
  $ref_crossprod->[0] =    $ref_v1->[1]*$ref_v2->[2] - $ref_v2->[1]*$ref_v1->[2];
  $ref_crossprod->[1] = -( $ref_v1->[0]*$ref_v2->[2] - $ref_v2->[0]*$ref_v1->[2] );
  $ref_crossprod->[2] =    $ref_v1->[0]*$ref_v2->[1] - $ref_v2->[0]*$ref_v1->[1];
        
}

sub _rotateByAxisAngle {

    my $ref_axis = shift;
    my $angle = shift; # in degrees
    my $ref_invec = shift;
    my $ref_outvec = shift;

    my $PI = 3.14159265359;
    my $angle_rad = $angle * $PI / 180;

    my $x = $ref_axis->[0];
    my $y = $ref_axis->[1];
    my $z = $ref_axis->[2];
    my $cosTheta=cos($angle_rad); 
    if (abs($cosTheta) < 1e-6) {$cosTheta=0;} 
    if (abs($cosTheta-1) < 1e-6) {$cosTheta=1;}
    my $sinTheta=sin($angle_rad); 
    if (abs($sinTheta) < 1e-6) {$sinTheta=0;}
    if (abs($sinTheta-1) < 1e-6) {$sinTheta=1;}

    my $oneminusCosTheta = (1-$cosTheta);

    my $r11 = $x*$x*$oneminusCosTheta +    $cosTheta;
    my $r12 = $x*$y*$oneminusCosTheta - $z*$sinTheta; 
    my $r13 = $x*$z*$oneminusCosTheta + $y*$sinTheta; 

    my $r21 = $y*$x*$oneminusCosTheta + $z*$sinTheta;
    my $r22 = $y*$y*$oneminusCosTheta +    $cosTheta;
    my $r23 = $y*$z*$oneminusCosTheta - $x*$sinTheta;

    my $r31 = $z*$x*$oneminusCosTheta - $y*$sinTheta;
    my $r32 = $z*$y*$oneminusCosTheta + $x*$sinTheta;
    my $r33 = $z*$z*$oneminusCosTheta +    $cosTheta;

    my @rotmat = ([$r11,$r12,$r13],[$r21,$r22,$r23],[$r31,$r32,$r33]);
     _multiply_matvec ( \@rotmat, $ref_invec, $ref_outvec );
}
	
sub _multiply_matvec {

    my $ref_mat = shift;
    my $ref_invec = shift;
    my $ref_outvec= shift;

    for(my $i=0;$i<3;$i++) {

        my $t = 0;
	
        for(my $j=0;$j<3;$j++) {
	    $t += $ref_mat->[$i][$j] * $ref_invec->[$j];
	      
	}   
	$ref_outvec->[$i] = $t;
    }
}

1;
