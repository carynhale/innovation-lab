#!/usr/bin/env perl

use warnings;
use strict;
use Getopt::Std;

my %opt;
getopts('h:s:e:l:', \%opt);

my $i = 0;
while (<STDIN>) {
    chomp;
    if ($i % 2 == 0) {
        print;
    } else {
        my $ss = $_;
        if ($opt{s}) {
            $ss = substr($ss, $opt{s});
        }
        if ($opt{e}) {
            $ss = substr($ss, 0, length($ss) - $opt{e});
        }
        if ($opt{l}) {
            $ss = substr($ss, 0, $opt{l});
        }
        print $ss;
    }
    print "\n";
    $i++;
}
