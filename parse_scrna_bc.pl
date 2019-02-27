### Song Cao ####

### Extract barcodes for each reads ### 

## last updated: Jan 30, 2019 ###

#!/usr/bin/perl
use strict;
use warnings;
my $usage = ' 
perl parse_scrna_bc.pl $f_in $f_out
';

die $usage unless @ARGV == 2;
my ($f_in,$f_out)=@ARGV;

open(IN,"<$f_in"); 
open(OUT,">$f_out"); 

my %rc_ref=(); 
my %rc_var=(); 
my %varlist=();
my @temp; 
my $count; 
 
while(<IN>)
{
	my $l=$_; 
	chomp($l); 
	@temp=split("\t",$l); 
	if($l=~/Support/) {
 	print OUT $temp[0],"\t",$temp[1],"\t",$temp[2],"\t",$temp[3],"\t",$temp[5],"\n";  
	}

	else 
	{
	my $find_b=0;  
	my $barcodes; 
	for(my $i=0;$i<scalar @temp; $i++) { 
	if($temp[$i]=~/^CB:Z:/) {
	$find_b=1; 
	$barcodes=$temp[$i]; $barcodes=~s/^^CB:Z://g; 
	last; 	
	}
	}	
	if($find_b==1) { print OUT $temp[0],"\t",$barcodes,"\n";} 
	}
 
}
close IN; 
close OUT; 
