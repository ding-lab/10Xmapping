### Song Cao ####

### generate number of reads supporting var and ref ### 

## last updated: Feb 8, 2019 ###

#!/usr/bin/perl
use strict;
use warnings;
my $usage = ' 
perl var_heatmap_cell.pl $ctype $f_in $f_maf $f_out
';

die $usage unless @ARGV == 4;

## f_in: add_celltype_bc output file ##
## ctype: you would like to show, like Plasma ##

my ($ctype,$f_in,$f_maf,$f_out)=@ARGV;

open(IN,"<$f_in"); 
open(MAF,"<$f_maf"); 
open(OUT,">$f_out"); 

#my %vaf=(); 
my %mutlist=(); 
my %bclist=(); 
my %var2pr=(); 
my @temp; 
my $var_id; 
my $var_s;
my %rc_var=(); 
my %rc_ref=(); 
my %exist=(); 
#my %mutlist=();
#my %bclist=();

while(<MAF>)
{
        my $l=$_;
        chomp($l);
        @temp=split("\t",$l);
        my $pr=$temp[36];
        $pr=~s/^p\.//g;
        my $gn=$temp[0];
        my $cb=$gn."-".$pr;
        my $chr=$temp[4];
        $chr=~s/chr//g;
        my $var=$chr."_".$temp[5]."_".$temp[10]."_".$temp[12];
        if($temp[0]=~/^IGH/) { next; }
        $var2pr{$var}=$cb;
}
 
while(<IN>)
{
        my $l=$_;
        chomp($l);
        @temp=split("\t",$l);
        if($l=~/Support/) {

        #print OUT $sn,"\t",$temp[0],"\t",$temp[1],"\t",$temp[2],"\t",$temp[3],"\t",$temp[4],"\n";

        $var_id=$temp[1]."_".$temp[2]."_".$temp[3]."_".$temp[4];

        if($temp[5]=~/Var/) {  $var_s=1; }
        if($temp[5]=~/Ref/) {  $var_s=0; }

        }

        else
        {
	if(defined $var2pr{$var_id})
	{ 
	#print $var_s,"\n"; <STDIN>;
       # if(defined $ctype{$sn} && defined $ctype{$sn}{$temp[1]})
       # {
       # my $ct=$ctype{$sn}{$temp[1]};
       # print OUT $temp[0],"\t",$temp[1],"\t",$ctype{$sn}{$temp[1]},"\n";

	if($temp[2] eq $ctype || lc($ctype)=~/all/) 
	{
        if($var_s==1) {  $rc_var{$var_id}{$temp[1]}++; $exist{$var_id}{$temp[1]}=1; $mutlist{$var_id}=1; $bclist{$temp[1]}=1;  }
        else { $rc_ref{$var_id}{$temp[1]}++; $exist{$var_id}{$temp[1]}=1; $mutlist{$var_id}=1; $bclist{$temp[1]}=1; }
	}

       }
       }
       # }
       # }
}
 

print OUT "Mutatation"; 
foreach my $bc (sort keys %bclist) 
{
	print OUT ",",$bc; 
}
print OUT "\n";

foreach my $m (sort keys %mutlist) 
{

 if(defined $var2pr{$m}) 
 {
  print OUT $var2pr{$m}."-Ref"; 
  
  foreach my $bc (sort keys %bclist)
  {
        if(defined $rc_ref{$m}{$bc}) { print OUT ",",$rc_ref{$m}{$bc}; }
        else { if(defined $rc_var{$m}{$bc}) { print OUT ",","0"; }  else { print OUT ","," "; } }
  }

  print OUT "\n";

   print OUT $var2pr{$m}."-Var";  

  foreach my $bc (sort keys %bclist)  
  {
	if(defined $rc_var{$m}{$bc}) { print OUT ",",$rc_var{$m}{$bc}; }
	else { if(defined $rc_ref{$m}{$bc}) { print OUT ",","0"; }  else { print OUT ","," "; } } 
  }	
 
  print OUT "\n"; 			
 }

} 

close IN;
close MAF;  
close OUT; 
