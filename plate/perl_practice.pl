#!/usr/bin/perl -w
use strict;
use YAML;

=head
foreach(split /\n/, `/bin/cat index.html`){
	my @cols = /<h3 class=".*">(.*)<\/h3>/g;
	print Dump @cols;
}
=cut

$_ = `/bin/cat index.html`;
my @cols = /<h3.*?>(.*?)<.*?h3>/g;
print Dump @cols;









# vi:sw=4:ts=4:nowrap
