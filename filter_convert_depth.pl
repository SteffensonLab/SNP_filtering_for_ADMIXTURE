=pod

The perl is to convert vcf file for the format : GT:DP:AD

in addition:
(1) only bi-allel variants were kept;

(2) The non-polymorphic site (such as ancient DNA data) will be kept because the data in the following analysis may be used to be emerged with known variants

tolerance: threshold for reads depths of heterozygous calls. 0.2 means allele_1/allele_2>=0.2 or <=0.8

sample_info: including 4 columns
column1: sample ID
column2: minimum reads depth for a homozygous call
column3: minimum reads depth for both alleles of a heterozygous call
column4: maximum reads depth for a genotype call

=cut

#!/usr/bin/perl -w
use strict;
use warnings;

die '@ARGV is required' if @ARGV !=4;

my $tolerance=shift; #  for heterozygous calls. 0.2 means allele_1/allele_2>=0.2 or <=0.8
my $sample_info=shift; # s1 2 2 20 
my $vcf=shift; 
my $output=shift;

my %x;
open IN,$sample_info or die $!;
while (<IN>)
{
    chomp;
    my @a=split;
    push @{$x{$a[0]}},$a[1],$a[2],$a[3];
}
close IN;

open IN,"gzip -dc $vcf | grep -v \"^##\"|" or die $!;
open OUT,"| gzip > $output";

my $t=<IN>;
print OUT "$t";

chomp $t;
my @t=split /\s+/,$t;
my %y;
for (my $i=9;$i<@t;$i++)
{
    die "$t[$i] has no sample information" if not exists $x{$t[$i]};
    @{$y{$i}}=@{$x{$t[$i]}}; # 2 2 20
}

while (<IN>)
{
    chomp;
    my @a=split;
    next if $a[4]=~/,/;
    $a[7]=".";
    $a[8]="GT:DP:DV";
    for (my $i=9;$i<@a;$i++)
    {
        if ($a[$i]=~/^\.\/\./){$a[$i]="./.:.:."}   
        else
        {
            my $dp1=$y{$i}[0];
            my $dp2=$y{$i}[1];
            my $dp3=$y{$i}[2];

            my @b=split /:/,$a[$i];
            if ($b[2]<1 or $b[2]>$dp3){$a[$i]="./.:.:."}
            else
            {
                my @c=split /,/,$b[3];
                my $ref=$c[0];
                my $alt=$c[1];
                my $rate_ref=$ref/$b[2];
                if ($rate_ref>$tolerance)
                {
                    if ($ref>=$dp1){$a[$i]="0/0:$b[2]:$alt"}
                    else{$a[$i]="./.:.:."}
                }
                elsif($rate_ref<(1-$tolerance))
                {
                    if ($alt>=$dp1){$a[$i]="1/1:$b[2]:$alt"}
                    else{$a[$i]="./.:.:."}
                }
                else
                {
                    if ($ref>=$dp2 && $alt>=$dp2){$a[$i]="0/1:$b[2]:$alt"}
                    else{$a[$i]="./.:.:."}
                }
            }
        }
    }
    my $ln=join "\t",@a;
    print OUT "$ln\n";
}            
close IN;
close OUT;
