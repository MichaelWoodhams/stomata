#!/usr/bin/env perl
use strict;
use warnings;

# This is a post-hoc fix to a problem which
# will get a proper fix elsewhere, so is not needed long term.

# I did a bunch of runs of runSim.R, output many text files.
# These text files do not record in their body which
# measurement, climate or bdratio they ran on.
# This code will create a combined table which includes this information.

open(OUT,">combined.txt") or die;
print OUT "measure bdratio climate rep depth rand.corr1 rand.corr2 rand.corr3 rand.corr4 rand.corr5 rand.sd1 rand.sd2 rand.sd3 rand.sd4 rand.sd5 rand.sigma2 rand.ace.sigma2 rand.ace.se\n";

my @files = glob("results/results*txt");
for my $file (@files) {
    my @parse = split("_",$file);
    my $measure = $parse[1];
    my $bdratio = $parse[2];
    my $climate = $parse[5];
    open(IN,$file) or die;
    my $header = <IN>; # throw away header
    my @lines = <IN>;
    close(IN);
    for my $line (@lines) {
	print OUT "$measure $bdratio $climate $line";
    }
}
close(OUT);
