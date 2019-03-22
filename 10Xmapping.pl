
### Song Cao ####

### Extract supporting reads for ref and var for somatic variants from bam ### 

## last updated: Jan 30, 2019 ###

#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

(my $usage = <<OUT) =~ s/\t+//g;
perl 10Xmapping.pl --mapq --bam --maf --out

mapq: mapping quality: default, 0 
bam: input bam 
maf: input maf
out: output

OUT

my $p_mapq=0; 
my $f_bam=""; 
my $f_maf=""; 
my $f_out=""; 

my $status = &GetOptions (
      "mapq=i" => \$p_mapq,
      "bam=s" => \$f_bam,
      "maf=s"  => \$f_maf,
      "out=s"  => \$f_out,
    );

if ($f_bam eq "" || $f_maf eq "" || $f_out eq "" ) {
      print $usage;
      exit;
   }

open(OUT,">$f_out");


foreach my $l (`cat $f_maf`)
	{
		my $ltr=$l; chomp($ltr); 
		my @temp=split("\t",$ltr);
	#	if(length($temp[15])<15) { $temp[12].="-01"; }
		my $sn=$temp[15];
		$sn=~s/\_T//g;  
		#print $sn,"\n";
		my $chr=$temp[4];
		$chr=~s/chr//g;  

		my $pos=$temp[5];
		my $ref=$temp[10];
		my $var=$temp[12];
		my $dellen=0;
 		my $inslen=0; 	

		if($var eq "-") { $dellen=length($ref); } 
		if($ref eq "-") { $inslen=length($var); }

		#print $chr,"\t",$pos,"\n";
		my $chr_s=1; 
		foreach my $h (`samtools view -H $f_bam`) 
		{
			my $htr=$h;
			chomp($htr); 
			if($htr=~/SN:1/) { $chr_s=0; last; }
			if($htr=~/SN:chr1/) { $chr_s=1; last; }	
		}

		if($chr_s==1) { $chr="chr".$chr; }

		my $chr_pos=$chr.":".$pos."-".$pos; 

		#print $chr_pos,"\t",$ref,"\t",$var,"\n";
 		#<STDIN>;		
	
		if($f_bam ne "NULL" && (-e $f_bam)) 
		{

			my %count_read_ref=();
			my %count_read_var=();
 
			foreach my $t (`samtools view $f_bam \"$chr_pos\"`)
			{
				chomp($t);	
				#print $t,"\n";
				#<STDIN>;	
				my @temp2=split("\t",$t);
				my $cigar=$temp2[5]; 
				my @intstr=split /(\D+)/, $cigar;
				my $find_var=0; 
				my $find_ref=0;
				my $seq=$temp2[9]; 
	 			my $startp=$temp2[3];
				my $index=0; 
				my $mapq=$temp2[4];
				my $id=$temp2[0]; 
				my $flag=$temp2[1];
				my $rid; 
				if($mapq<$p_mapq) { next; } 
                                #print $t,"\n";
                                #<STDIN>;
 				if($dellen==0 && $inslen==0)
				{					
					for(my $i=0;$i<scalar @intstr;$i+=2)
					{
					### junction ##
				#	print $intstr[$i],"\t",$intstr[$i+1],"\n";

				 	if($intstr[$i+1] eq "N") 
					{ 
					  $startp+=$intstr[$i];  
					}
 
					### deletion and insertion, then not a snv, over ##
				 	if($intstr[$i+1] eq "D" || $intstr[$i+1] eq "I") { $find_ref=0; $find_var=0; last; } 
					## soft clipping ##
					if($intstr[$i+1] eq "S") { $index+=$intstr[$i]; }
				   	if($intstr[$i+1] eq "M") { 	
					for(my $j=0; $j<$intstr[$i];$j++) {
					my $nt=substr($seq,$index,1); 
					#print $startp,"\t",$index,"\t",$nt,"\n"; 
					#<STDIN>;
					if($startp eq $pos && $nt eq $var) 
					{
						$find_var=1;
						#print "Finding Variant","\t",$find_var,"\t",$startp,"\t",$nt,"\n"; 
						#print $id,"\n"; 
						#<STDIN>;
					} 
					if($startp eq $pos && $nt eq $ref)  
                                        {
                                                $find_ref=1;
						#print "Finding Reference","\t",$find_var,"\t",$startp,"\t",$nt,"\n";
                                                #print $id,"\n"; 
                                        }
					$startp++; $index++; 		
					}
					}					
					}
				 } ## end if snv 
				### insertion 

				elsif ($dellen==0 && $inslen>0) 
				{
					$find_ref=1;
				    for(my $i=0;$i<scalar @intstr;$i+=2)
                    {
				    if($intstr[$i+1] eq "N")
                    {
                      $startp+=$intstr[$i];
                    }
					if($intstr[$i+1] eq "S") { $index+=$intstr[$i]; }
					if($intstr[$i+1] eq "M") {
					#for(my $j=0; $j<$intstr[$i];$j++) {
					
                    #my $nt=substr($seq,$index,1);
					#if($startp eq $pos && $nt eq $ref) 
					#{
					#  $find_ref=1;
					#  print $t,"\n"; <STDIN>;	  		
					#}	
					$startp+=$intstr[$i]; $index+=$intstr[$i];	
					}
					
					if($intstr[$i+1] eq "I")  {
					 #$find_ref=0;	
					 if($startp eq $pos+1 && $intstr[$i] eq $inslen) 
					 {
						$find_var=1; $find_ref=0; 				
					 }
					}
					}
					for(my $i=0;$i<scalar @intstr;$i+=2)
                    {
					 if($intstr[$i+1] eq "D") { $find_ref=0; $find_var=0; last; }		
					 if($intstr[$i+1] eq "I") { $find_ref=0; last; }  	
					}
				} ## elseif insertion 	

		
				### deletion ##

			    elsif ($dellen>0 && $inslen==0)
                {
                    $find_ref=1;
                    for(my $i=0;$i<scalar @intstr;$i+=2)
                    {
                    if($intstr[$i+1] eq "N")
                    {
                      $startp+=$intstr[$i];
                    }
                    if($intstr[$i+1] eq "S") { $index+=$intstr[$i]; }
                    if($intstr[$i+1] eq "M") {
                    #for(my $j=0; $j<$intstr[$i];$j++) {

                    #my $nt=substr($seq,$index,1);
                    #if($startp eq $pos && $nt eq $ref) 
                    #{
                    #  $find_ref=1;
                    #  print $t,"\n"; <STDIN>;          
                    #}  
                    $startp+=$intstr[$i]; $index+=$intstr[$i];
                    }

                    if($intstr[$i+1] eq "D")  {
                     #$find_ref=0;  
                     if($startp eq $pos && $intstr[$i] eq $dellen)
                     {
                        $find_var=1; $find_ref=0;
                     }
                    }
                    }
                    for(my $i=0;$i<scalar @intstr;$i+=2)
                    {
                     if($intstr[$i+1] eq "I") { $find_ref=0; $find_var=0; last; }
                     if($intstr[$i+1] eq "D") { $find_ref=0; last; }
                    }
                } ## elseif insertion 
				#if($id=~/\/2$/ || ($flag & 0x80)) { $rid=$id; $rid=~s/\/2$//g; $rid.="\/2";  }
				#if($id=~/\/1$/ || ($flag & 0x40)) { $rid=$id; $rid=~s/\/1$//g; $rid.="\/1"; } 
				if($find_ref==1) { $count_read_ref{$id}=$t; }
				if($find_var==1) { $count_read_var{$id}=$t; }		 
					}	
			
		
			my $n_key_ref=keys %count_read_ref;
			print OUT $chr,"\t",$pos,"\t",$ref,"\t",$var,"\t",$n_key_ref,"\t","Ref-Support","\n"; 
			foreach my $rr (sort keys %count_read_ref)
				{
					print OUT $count_read_ref{$rr},"\n";
				}
			my $n_key_var=keys %count_read_var;
                        print OUT $chr,"\t",$pos,"\t",$ref,"\t",$var,"\t",$n_key_var,"\t","Var-Support","\n";
                        foreach my $rr (sort keys %count_read_var)
                                {
                                        print OUT $count_read_var{$rr},"\n";
                                }
		} ## if f_bam ##
		} ## for maf ##

close OUT;
