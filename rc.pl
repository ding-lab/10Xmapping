### Song Cao ####

### Extract supporting reads for ref and var for somatic variants from bam ### 

## last updated: Jan 30, 2019 ###

#!/usr/bin/perl
use strict;
use warnings;
my $usage = ' 
perl rc.pl $f_in $f_out
';

die $usage unless @ARGV == 2;
my ($f_in,$f_out)=@ARGV;
open(IN,"<$f_in"); 
open(OUT,">$f_out"); 

my %rc_ref=(); 
my %rc_var=(); 
my %varlist=();
my @temp; 
my $var; 
my @vars;
 
while(<IN>)
{
	my $l=$_; 
	chomp($l); 
	@temp=split("\t",$l); 


	if($l=~/Support/) { 
    $var=$temp[0]."_".$temp[1]."_".$temp[2]."_".$temp[3];

    if(!defined $varlist{$var}) { push @vars, $var; $varlist{$var}=1; }

	if($l=~/Ref\-Support/) { $rc_ref{$var}=$temp[4]; }
    if($l=~/Var\-Support/) { $rc_var{$var}=$temp[4]; }	
	}
 
}

foreach my $v (@vars) 
 {
	my $vtr=$v; 
	chomp($vtr); 
	@temp=split(/\_/,$vtr);
	#print $v,"\n"; 
	if($rc_var{$v}+$rc_ref{$v}==0)
	{
	print OUT $temp[0],"\t",$temp[1],"\t",$temp[2],"\t",$temp[3],"\t",$rc_ref{$v},"\t",$rc_var{$v},"\t","0","\n"; 	
	}
	else 
	{
	print OUT $temp[0],"\t",$temp[1],"\t",$temp[2],"\t",$temp[3],"\t",$rc_ref{$v},"\t",$rc_var{$v},"\t",$rc_var{$v}/($rc_var{$v}+$rc_ref{$v}),"\n";
	}
 }

