### Song Cao ####

### add cell type  ### 

## last updated: Feb 8, 2019 ###

#!/usr/bin/perl
use strict;
use warnings;
my $usage = ' 
perl add_celltype_bc.pl f_ctype f_in f_out
';

### f_ctype: cell type and barcode table ##

### f_in output file from parse_scrna_bc.pl ##


die $usage unless @ARGV == 3;

my ($f_ctype,$f_in,$f_out)=@ARGV;

my $f_count=$f_out.".rc"; 

open(TYPE,"<$f_ctype"); 
open(IN,"<$f_in"); 
open(OUT,">$f_out"); 
open(OUT1,">$f_count"); 

my %ctype=(); 
my %rc_ref=(); 
my %rc_var=(); 
my %varlist=();
my %exist=();
my @temp; 
my $count; 
my $dir_l=(split(/\//,$f_in))[-1];
my $sn=(split(/\./,$dir_l))[0];
print $sn,"\n";

while(<TYPE>) 
{
	my $l=$_; 
	chomp($l); 
	@temp=split("\t",$l); 
 	if($l=~/^sample_id/) { next; } 
	else {
	$ctype{$temp[0]}{$temp[1]}=$temp[2];  
	}		
} 

my $var_id; 
my $var_s=0;
 
while(<IN>)
{
	my $l=$_; 
	chomp($l); 
	@temp=split("\t",$l); 
	#my $var;	
	if($l=~/Support/) {
 	print OUT $sn,"\t",$temp[0],"\t",$temp[1],"\t",$temp[2],"\t",$temp[3],"\t",$temp[4],"\n";  
	$var_id=$temp[0]."_".$temp[1]."_".$temp[2]."_".$temp[3]; 
	if($temp[4]=~/Var/) {  $var_s=1; }
	if($temp[4]=~/Ref/) { $var_s=0; }
	}
	else 
	{
	if(defined $ctype{$sn} && defined $ctype{$sn}{$temp[1]}) 
	{
	my $ct=$ctype{$sn}{$temp[1]}; 
	print OUT $temp[0],"\t",$temp[1],"\t",$ctype{$sn}{$temp[1]},"\n";				
	if($var_s==1) {  $rc_var{$var_id}{$ct}++; $exist{$var_id}{$ct}=1; }
	else { $rc_ref{$var_id}{$ct}++; $exist{$var_id}{$ct}=1; }
	}
	} 
}

foreach my $v (sort keys %exist) 
{
	foreach my $c (sort keys %{$exist{$v}}) 
	{

	print OUT1 $sn,"\t",$v,"\t",$c; 

	if(defined $rc_ref{$v} && defined $rc_ref{$v}{$c}) 
	{
	print OUT1 "\t",$rc_ref{$v}{$c}; 
	}
	else {	print OUT1 "\t","0"; }

	if(defined $rc_var{$v} && defined $rc_var{$v}{$c})
        {
        print OUT1 "\t",$rc_var{$v}{$c},"\n";
        }       
        else {  print OUT1 "\t","0","\n"; }		
	}
}

close IN; 
close OUT; 

