#!/usr/bin/perl -w

use strict;


# AT5G53240.1   Scaffold00032   81.91   282 42  2   403 681 3630975 3631250 6e-35    153

open("BR","<".$ARGV[0]);

my $cur_q ='';
my $best_eval=10;
my $best_score=0;

while(<BR>) {
    chomp;
    my ($q,$h,$pid,$a,$b,$c,$q_s,$q_e,$s_s,$s_e,$eval,$bit) =(split("\t",$_));

    if($cur_q eq $q) {
#        if($best_eval > $eval/100 || $eval == 0) {
#        if ($eval < 1e-20) {
        if($eval < 1e-20 && $bit/$best_score > 0.6){
            print $_."\n";
        }
    }
    else {
        $cur_q = $q;
        $best_eval = $eval;
        $best_score = $bit;
        print $_."\n";
    }
}