#!/usr/bin/perl -w
use strict; use YAML;use JSON;
#regex: check pattern, extract things from string, replace
# how about extract things in a ()?



#�Υ~�����O�L�X
$_=`cat ~/public_html/games/2015-04-05/Cardinals-Cubs`;
#�����Xtable�̪����
my @table = /<table class="data" width="100%" >(.*)<\/table>/g;
#die scalar @table;
#��ŦX�����e�s���y���m�W�}�C
my @name = /playerpage\/\d+">(.*?)<\/a>/g;
my @score = /<td align="right" width="\d+%">(.*)<\/td>/g;

#die scalar @score;
#print Dump @_; 
#print Dump join '|',@name;

for(my $i=0; $i<=$#name; $i++) {
    print "$i $name[$i]\n";
}
# vi:sw=4:ts=4
