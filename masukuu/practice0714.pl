#!/usr/bin/perl -w
use strict; use YAML;
#regex: check pattern, extract things from string, replace

# how about extract things in a ()?

=head
open FH,'index.html';
while(<FH>){
  @_ = /[^"normal"]>(.*)<\/h3>/g;
  print Dump @_; 
}
=cut

#用外部指令印出index.html
$_=`cat index.html`;
@_ = /[^"normal"]>(.*)<\/h3>/g;
print Dump @_; 

# vi:sw=4:ts=4
