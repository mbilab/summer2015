#!/usr/bin/perl -w
use strict; use YAML;
# regex: check pattern, extract things from string, replace

# check pattern: often used in input validation 

## if account need to be the combination of lower-case alphabet and numbers
$_ = "andy23512";
/^[a-z0-9]+$/ and print "account!\n"; 

## check if the string is an interger
$_ = "-1";
/^-?\d+$/ and print "interger!\n";

## pratice: check if the string is a number (interger and floating point) 
### try it yourself !!!
$_ = "-10.0";
/-?\d+(\.\d+)?/ and print "number!\n";

# extraction

$_ = "1bc465dd";
# extract single lower-case letter
@_ = /[a-z]/g;
print Dump join '|', @_;

# how about extract things in a ()?
$_ = "Tien-Hao Chang (dirty), Che-Wei Chang (tangent)";
@_ = /\((.*)\)/g;
print Dump join '|', @_; # what the ...?
# * is greedy!!!!!!
@_ = /\((.*?)\)/g; # *? => not greedy
print Dump join '|', @_; 
# another possible solution
@_ = /\(([^)]+)/g;
print Dump join '|', @_; 

# when to use /.../s?
$_ = "
(\n
dirty is handsome.\n
)\n
(\n
I am dirty but I am not dirty.\n
)\n
";
@_ = /\((.*?)\)/g;
print Dump join '|', @_; # nothing is extracted?
@_ = /\((.*?)\)/gs;
#print $_[0]; # nothing is extracted?
#print $_[1]; # nothing is extracted?


# vi:sw=4:ts=4
