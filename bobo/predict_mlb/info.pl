#!/usr/bin/perl -w
use strict;
use YAML;

$_=`cat ../games/2015-04-05/Cardinals-Cubs`;

my @name=/[playerpage]\/\d+">(.*?)<\/a>.*?[right].*?(\d+)<\/td>/g;
my @grade=/[right][^Score].*?(\.?\d+)<\/td>/g;
print Dump join '|',@grade;
#print Dump @b;

#open FH,'<index.html';
#while(<FH>){
#    @_=/[^"normal"]>(.*)<\/h3>/g;
#    print Dump @_;
#}
#close FH;


# vi:sw=4:ts=4
