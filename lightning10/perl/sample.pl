#!/usr/bin/perl -w
use strict; # force variable delaration
use YAML; # we only use the Dump ... 
# Yet Another Markup Language, but it's not important
# other module, you can find them in cpan http://www.cpan.org/

if($ARGV[0] eq 'hello'){
	print "Hello world\n";
}
# variable
## scalar
my $scalar; # for number, string, reference

# single-line comment
=head
	multi-line comment
	http://perldoc.perl.org/perlpod.html
=cut

## array
if($ARGV[0] eq 'array'){
	my @array = (1, 2, 3); # array assignment: normal
	my @a = qw/a b c d e/; # array assignment: string array ('a', 'b', 'c', 'd', 'e')
	print '@a: (a, b, c, d, e)';
	print "\n";
	print 'length of @a: ';
	print scalar @a; # length of an array
	print "\n";
	print '$a[1]: ';
	print $a[1]; # get content of certain index of an array
	print "\n";
	print '$#a: ';
	print $#a; # get the index of the last entry of an array
	print "\n";
	# @a[1..3]: ('b', 'c', 'd') # a part of array  
	# push @a, 4;
	# pop @a;
	# shift @a;
	# unshift @a, 4;
}
if($ARGV[0] eq 'hash'){
	my %h = (a=>15, b=>16, c=>17); # hash assignment: normal
	%h = ('a', 15, 'b', 16, 'c', 17); # hash assignment: another way 
	%h = qw/a 15 b 16 c 17/; # hash assignment: yet another way 
	print '%h: (a=>15, b=>16, c=>17)';
	print "\n";
	print '$h{c}: '; 
	print $h{c}; # get value in hash by a key
	print "\n";
	print 'keys %h: ';
	print join ', ', keys %h; # get all keys of a hash 
}
if($ARGV[0] eq 'for'){
	my @a = qw/a b c d e/;
	# numeric for loop
	# for(0 .. $#a){print "$_: $a[$_]\t";}
	# for my $i(0 .. $#a){print "$i: $a[$i]\t";}
	print "\n";
	# array traverse foreach loop
	# foreach(@a){print "$_\t";}
	# actually... for==foreach in perl
	for(0 .. 10){
		if($_ < 2) { next; } # next: continue in c
		print "$_\n";
		if(5 == $_) { last; } # last: break in c
	}
}
if($ARGV[0] eq 'cond'){
	# ($ARGV[1])? print 'true': print 'false';
	print $ARGV[1]?'true':'false';
	print "\n";
	$ARGV[1] and print 'true and';
	$ARGV[1] or print 'false or';
}
if($ARGV[0] eq 'sub'){
	print &example(1, 3, 5, 7);
}
sub example{
	print join ', ', @_; # your parameters are in @_
	print "\n";
	print shift; # shift here would shift one item in @_
	print "\n";
	print join ', ', @_;
	return 9;
	# 9; # actually, the last line would be returned automatically
}
if($ARGV[0] eq 'file'){
	open FH, '<10.gtf'; # read a file
	while(<FH>){ # use <FH> to get a line from file
		print Dump "$_"; # $_ is the line we get, note that \n is in $_
		chomp; # remove "\n" in $_
		print Dump "$_";
	}
	close FH; # remember to close it after use
	open FH, '>output'; # write a file (overwrite!)
	for(0 .. 12){
		print FH $_**2; # no comma!!
		print FH "\n";
	}
	close FH;
	# open FH, '>>output' # append to a file
	# for(13 .. 15){
	#	print FH $_**2;
	#	print FH "\n";
	# }
}
if($ARGV[0] eq 'map'){
	# format: map {function} list
	print Dump map {$_*2} (0 .. 18);
}
if($ARGV[0] eq 'sort'){
	# format: sort {sort function} list
	print join ',', sort (0 .. 18); # string
	print "\n";
	print join ',', sort {$a cmp $b} (0 .. 18); # string
	print "\n";
	print join ',', sort {$a <=> $b} (0 .. 18); # numeric
	print "\n";
	# reverse: use sort function {$b <=> $a}
}
if($ARGV[0] eq 'shell'){
	print `/bin/cat output`; # alias may cause problem!!!
}
if($ARGV[0] eq 'argv'){
	print join ',', @ARGV;
	print "\n";
	my ($a, $b, $c) = @ARGV; # destructure;
	# my ($a, $b, $c) = (shift, shift, shift); # destructure;
	print "\$a=$a, \$b=$b, \$c=$c\n";
	shift; # shift here would shift one item in @ARGV;
	print join ',', @ARGV;
}
if($ARGV[0] eq 'split'){
	foreach(split /\n/, `/bin/cat 10.gtf`){ # another way to read a file
		my @cols = split /\t/; # split here is split $_ by \t
		die Dump @cols; # die here end the program, Dump help us print a structure
	}
}
# most difficult part so far: reference
if($ARGV[0] eq 'reference'){
	my $array_ref = ['a', 'b', 'c']; # notice the symbol
	my @array =     ('a', 'b', 'c'); 
	print '$array_ref: ';
	print $array_ref; # you can see that it is an address 
	print "\n";
	print Dump @{$array_ref}; # convert it to an array
	print '$array_ref->[0]: ';
	print $array_ref->[0]; # how to get value
	my $array_ref2 = \@array; # or you can get the reference of array by this

	my $hash_ref = {a=>15, b=>16, c=>17}; # notice the symbol
	my %hash =     (a=>15, b=>16, c=>17);
	print '$hash_ref: ';
	print $hash_ref; # you can see that it is an address 
	print "\n";
	print Dump %{$hash_ref}; # convert it to a hash
	print '$hash_ref->{b}: ';
	print $hash_ref->{b}; # how to get value
	my $hash_ref2 = \%hash; # or you can get the reference of hash by this
}
if($ARGV[0] eq 'other'){
	open FH, '<not_exist' or die $!; # $! save the error msg
}
# complex structure
if($ARGV[0] eq 'nested'){
	use JSON;
	## array of hash
	my @editor = (
		{name=>'vi', desc=>'god of editor'}, # each entry is a hash ref
		{name=>'emacs', desc=>'god\'s editor'}
	);
	# print JSON->new->pretty(1)->encode(\@editor);
	## use push to add an item
	## push array, reference
	push @editor, {name=>'notepad', desc=>'editor of non-human'};
	# print JSON->new->pretty(1)->encode(\@editor);
	## traverse it
	foreach(@editor){
		# print "$_->{name}: $_->{desc}\n"; # note that now $_ is a hash ref here
	}
	
	## hash of array
	my %int = (
		even=> [0, 2, 4, 6], # each value is an array ref
		odd=> [1, 3, 5, 7]
	);
	# print JSON->new->pretty(1)->encode(\%int);
	## use assignment to add key
	$int{prime} = [2, 3, 5, 7];
	# print JSON->new->pretty(1)->encode(\%int);
	## traverse it
	foreach(keys %int){
		print "$_: ";
		print join ',', @{$int{$_}};
		print "\n";
	}
	# hash of hash, array of array can be constructed in same concept
	# even more complicated structure can be constructed
	# proper data structure lead to good performace!
}


# vi:sw=4:ts=4:nowrap
