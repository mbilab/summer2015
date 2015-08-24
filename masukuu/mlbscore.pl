#!/usr/bin/perl -w
use strict; use YAML;use JSON;
#regex: check pattern, extract things from string, replace
# how about extract things in a ()?



#用外部指令印出
$_=`cat ~/public_html/games/2015-04-05/Cardinals-Cubs`;
#先取出table裡的資料
my @table = /<table class="data" width="100%" >(.*)<\/table>/g;
#die scalar @table;
#把符合的內容存成球員姓名陣列
my @name = /playerpage\/\d+">(.*?)<\/a>/g;
my @score = /<td align="right" width="\d+%">(.*)<\/td>/g;

#die scalar @score;
#print Dump @_; 
#print Dump join '|',@name;

for(my $i=0; $i<=$#name; $i++) {
    print "$i $name[$i]\n";
}
# vi:sw=4:ts=4
