#!/usr/bin/perl -w
use strict; use YAML;

open FH, '<index.html'; # open file
my $text = '';
while(<FH>){$text = $text . "$_";}
close FH; # close file

my @t = $text =~ /<h3.*?h3>/g;
foreach(@t)
{
	/<h3.*?normal.*>/ and next;
	s/(<.*?>)//g; 
	print $_ . "\n";
}	

# vi:sw=4:ts=4
