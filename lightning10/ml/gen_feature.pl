#!/usr/bin/perl -w

#use warnings;
use strict;
use YAML;
use constant PI    => 4 * atan2(1, 1);
use Time::HiRes;
our $time;

=head
sub start{
	$time = Time::HiRes::time();
}
sub end{
	print $_[0].' ';
	my $tmp = Time::HiRes::time();
	print( ($tmp - $time)* 1000000);
	print " us\n";
}
=cut
# Smith-Waterman  Algorithm

# usage statement
#die "usage: $0 <sequence 1> <sequence 2> <result_file> <prediction>\n" unless @ARGV == 4;
my ($inputFile, $resultFile, $prediction) = @ARGV; 
#my ($inputFile, $resultFile, $prediction) = ('./debug/output_final','train_data','1');
# my ($inputFile, $resultFile, $prediction) = ('train/result/1','test_new','0');
my %opt = (
	"window-size" => 55 
);

my %di2num = ( "AA" => 0, "AT" => 1, "AC" => 2, "AG" => 3, "TA" => 4, "TT" => 5, "TC" => 6, "TG" => 7, "CA" => 8, "CT" => 9, "CC" => 10, "CG" => 11, "GA" => 12, "GT" => 13, "GC" => 14, "GG" => 15);

my %tri2num = ( "TCA" => 0, "TCC" => 1, "TCG" => 2, "TCT" => 3, "TTC" => 4, "TTT" => 5, "TTA" => 6, "TTG" => 7, "TAC" => 8, "TAT" => 9, "TAA" => 10, "TAG" => 11, "TGC" => 12, "TGT" => 13, "TGA" => 14, "TGG" => 15, "CTA" => 16, "CTC" => 17, "CTG" => 18, "CTT" => 19, "CCA" => 20, "CCC" => 21, "CCG" => 22, "CCT" => 23, "CAC" => 24, "CAT" => 25, "CAA" => 26, "CAG" => 27, "CGA" => 28, "CGC" => 29, "CGG" => 30, "CGT" => 31, "ATA" => 32, "ATC" => 33, "ATT" => 34, "ATG" => 35, "ACA" => 36, "ACC" => 37, "ACG" => 38, "ACT" => 39, "AAC" => 40, "AAT" => 41, "AAA" => 42, "AAG" => 43, "AGC" => 44, "AGT" => 45, "AGA" => 46, "AGG" => 47, "GTA" => 48, "GTC" => 49, "GTG" => 50, "GTT" => 51, "GCA" => 52, "GCC" => 53, "GCG" => 54, "GCT" => 55, "GAC" => 56, "GAT" => 57, "GAA" => 58, "GAG" => 59, "GGA" => 60, "GGC" => 61, "GGG" => 62, "GGT" => 63);

my $featureNUMBER = 264;
my $RESULT_OUTPUT_PATHWAY = $resultFile;

# my ($abc, $def) = (0, 0);
# process for miRNA file

open FH, $inputFile or die "input file $inputFile can't be read";
my @col;
my (@seq1Title, @miRNAseq, @seq2Title, @UTRseq);
# &start;
while (<FH>) {
	@col = split /\t/; 
	push @seq1Title, $col[0];
	push @miRNAseq, $col[1];
	push @seq2Title, $col[2];
	push @UTRseq, $col[3];
}
# &end('read');
my $sample_num = scalar @seq1Title;

my $MATCH    =  6; # +1 for letters that match
my $MISMATCH = -3; # -1 for letters that mismatch
my $GAP_OPENNING  = -9; # -1 for any gap
my $GAP_EXTENSION = -4;
my $WOBBLE_MATCH = 3;

open(FHDW, ">$RESULT_OUTPUT_PATHWAY.txt") or die "$!\n";
open (FHDW3, ">$RESULT_OUTPUT_PATHWAY.feature") or die "$!\n";
open (FHDW5, ">$RESULT_OUTPUT_PATHWAY.out") or die "$!\n";
open (FHHIT, ">$RESULT_OUTPUT_PATHWAY.hit") or die "$!\n";

open (FHMISS, ">$RESULT_OUTPUT_PATHWAY.miss") or die "$!\n";

my $pair_number = 0;
print FHDW5 "miRNA_name\tmRNA_name\tAlignment_num\tAlignment_starts\n";

for (my $i = 0; $i <= $#seq1Title && $i <= $#seq2Title; $i++) {
	# debug info
	print "Executing (".($i+1)."/$sample_num) ".sprintf("%.2f",($i+1)/$sample_num*100)." %\n";
	print "*Sequence1:$seq1Title[$i]\t*Sequence2:$seq2Title[$i]\n";
	print "start-";

	print FHDW "*Sequence1:$seq1Title[$i]\n*Sequence2:$seq2Title[$i]\n\n";
	my $miRNA_seed = substr($miRNAseq[$i], 1, 6);
	$miRNA_seed =~ tr/ACUG/UGAC/;
	print FHDW "miRNA seed: ".$miRNA_seed."\n";
	# &start;
	my @seed_pos = &kmer($miRNA_seed, $UTRseq[$i], 6, 1); 
	# &end("kmer");
	my (@feature, @predict_starts, @instance_feature, @prev_mfe, @filt_feature);
	my $j = 1;
	#seed_pos = kmer @hitall
	foreach (@seed_pos) {
		my $align_start = $_;
		next if ($align_start - length($miRNAseq[$i]) <= 15);
		my ($nearest_dis, $nearest_neighbour) = &nearest($align_start, @seed_pos);
		($nearest_dis eq 'inf') and next;
		my ($mfe, $pred_start, @align_feature) = &smithWaterman($miRNAseq[$i], $UTRseq[$i], $align_start, $nearest_dis);
		defined($mfe) or next; #! defined mfe or next
		my $again = 0;
		if($j > 1){ #! if predict_start exists
			for (my $k = 0; $k <= $#predict_starts; $k++) {
				if (abs ($pred_start - $predict_starts[$k]) <= 15) {
					$again = 1 if ($mfe >= $prev_mfe[$k]);
				}
			}
		}
		next if ($again == 1);
		push (@predict_starts, $pred_start);
		push (@prev_mfe, $align_feature[$featureNUMBER]);
		@instance_feature = (@instance_feature, @align_feature);
		$j++;
	}
	if(1 == $j){
		print FHMISS "$seq1Title[$i]\t$miRNAseq[$i]\t$seq2Title[$i]\t$UTRseq[$i]";
		next;
	}
	print FHHIT "$seq1Title[$i]\t$miRNAseq[$i]\t$seq2Title[$i]\t$UTRseq[$i]";
	my $predict_num = $j - 1;
	my $miRNA_name = $seq1Title[$i]; #! change name
	my $mRNA_name = $seq2Title[$i]; #! change name
	print FHDW5 "$miRNA_name\t$mRNA_name\t$predict_num\t";
	for (my $i = 0; $i <= $#predict_starts; $i++) { #! change end point
		print FHDW5 $predict_starts[$i];
		print FHDW5 "(".$prev_mfe[$i].")" if ($i <= $predict_num - 1);
		print FHDW5 "," if ($i < $predict_num - 1);
	}
	print FHDW5 "\n";
	print FHDW3 "$prediction";
	# my @feature;
	my $tradeoff = $instance_feature[$featureNUMBER];
	for (my $a = $featureNUMBER; $a >= 0; $a--) {
		my ($sum, $max_b) = (0, 0);
		
		$feature[$a] = 0;
		for (my $b = 0; $b < $predict_num; $b++) {
			if ($a == $featureNUMBER) {
				if ($prediction == 1) {
					if ($instance_feature[$a + ($featureNUMBER+1) * $b] < $tradeoff) {
						$tradeoff = $instance_feature[$a + ($featureNUMBER+1) * $b];
						$max_b = $b;
					}
				} elsif ($prediction == 0) {
					if ($instance_feature[$a + ($featureNUMBER+1) * $b] > $tradeoff) {
						$tradeoff = $instance_feature[$a + ($featureNUMBER+1) * $b];
						$max_b = $b;
					}
				}
			}
			$sum += $instance_feature[$a + ($featureNUMBER+1) * $b];
		}
		$feature[$a] = $sum if ($a < 12 || ($a > 134 && $a < 137));
		if (($a >= 12 && $a < 90)||($a >= 137 && $a < ($featureNUMBER-1))){
		    $feature[$a] = $sum / $predict_num if ($predict_num != 0);
		} elsif (($a >= 90 && $a <= 134)  ||$a == ($featureNUMBER+1)||$a == $featureNUMBER) {
			$feature[$a] = $instance_feature[$a + ($featureNUMBER+1) * $max_b] if ($instance_feature[$a + ($featureNUMBER+1) * $max_b] =~ /\d/);
		}
	}
	print FHDW3 map {" ".($_+1).":".$feature[$_]} 0..$featureNUMBER;	#print Feature vector
	print FHDW3 "\n";
	print "end\n";
}

close (FHDW);
close (FHDW3);
close (FHDW5);

sub smithWaterman {
	# &start;
	my ($miRNA, $UTRseq, $align_start, $distance) = @_;
	$miRNA = uc $miRNA;
	$UTRseq = uc $UTRseq;
	my $UTRseq_part = substr($UTRseq, $align_start - length($miRNA), length($miRNA) + 7);
	#my $seq2_start = $align_start - length($miRNA) + 1;
	#print $UTRseq_part."\t".$align_start."\n";
	$miRNA = reverse $miRNA;
	
	# initialization
	my (@matrix, @D_matrix, @I_matrix);
	$matrix[0][0]{score}   = 0;
	$matrix[0][0]{pointer} = "none";
	$I_matrix[0][0]{score} = 0;
	$D_matrix[0][0]{score} = 0;
	
	for (my $i = 0; $i <= length($UTRseq_part); $i++) {
		for(my $j = 0; $j <= length($miRNA); $j++) {
			$matrix[$i][$j]{score}   = 0;
			$matrix[$i][$j]{pointer} = "none";
			$D_matrix[$i][$j]{score} = 0;
			$I_matrix[$i][$j]{score} = 0;
		}
	}
	# fill
	my $max_i     = 0;
	my $max_j     = 0;
	my $max_score = 0;
	
	for(my $i = 1; $i <= length($UTRseq_part); $i++) {
		for(my $j = 1; $j <= length($miRNA); $j++) {
			my ($diagonal_score, $left_score, $up_score, $D_left_score, $I_up_score);
			
			# calculate match score
			my $letter1 = substr($miRNA, $j-1, 1);
			my $letter2 = substr($UTRseq_part, $i-1, 1);  
			
			# for similarity
			#if ($letter1 eq $letter2) {
			#	$diagonal_score = $matrix[$i-1][$j-1]{score} + $MATCH;
			#}
			#else {
			#	$diagonal_score = $matrix[$i-1][$j-1]{score} + $MISMATCH;
			#}
			
			# for complementation
			if (($letter1=~/A/ && $letter2=~/T|U/) || ($letter1=~/T|U/ && $letter2=~/A/) || ($letter1=~/C/ && $letter2=~/G/) || ($letter1=~/G/ && $letter2=~/C/)) {
				if ($j-1 < 8) {
					$diagonal_score = $matrix[$i-1][$j-1]{score} + 2 * $MATCH;
					$diagonal_score += 10 if ($j-1 == 0);
				} else {
					$diagonal_score = $matrix[$i-1][$j-1]{score} + $MATCH;
				}
			} elsif (($letter1=~/G/ && $letter2=~/T|U/) || ($letter1=~/T|U/ && $letter2=~/G/)) {
				$diagonal_score = $matrix[$i-1][$j-1]{score} + $WOBBLE_MATCH;
			} else {
				$diagonal_score = $matrix[$i-1][$j-1]{score} + $MISMATCH;
			}
			
			# calculate gap scores
			$up_score   = $matrix[$i-1][$j]{score} + $GAP_OPENNING;
			$left_score = $matrix[$i][$j-1]{score} + $GAP_OPENNING;
			$D_left_score = $D_matrix[$i][$j-1]{score} + $GAP_EXTENSION;
			$I_up_score = $I_matrix[$i-1][$j]{score} + $GAP_EXTENSION;
			#print "\t$up_score\t\tI:$I_up_score\n$left_score\t?\tD:$D_left_score\n\n";
			if ($diagonal_score <= 0 and $up_score <= 0 and $left_score <= 0) {
				$matrix[$i][$j]{score}   = 0;
				$matrix[$i][$j]{pointer} = "none";
				$D_matrix[$i][$j]{score} = 0;
				$I_matrix[$i][$j]{score} = 0;
				next; # terminate this iteration of the loop
			}
        
			# choose best score
			if ($D_left_score >= $left_score) {
				$D_matrix[$i][$j]{score} = $D_left_score;
				$left_score = $D_left_score;
			} else {
				$D_matrix[$i][$j]{score} = $left_score;
			}
			if ($I_up_score >= $up_score) {
				$I_matrix[$i][$j]{score} = $I_up_score;
				$up_score = $I_up_score;
			} else {
				$I_matrix[$i][$j]{score} = $up_score;
			}
			if ($diagonal_score >= $up_score) {
				if ($diagonal_score >= $left_score) {
					$matrix[$i][$j]{score}   = $diagonal_score;
					$matrix[$i][$j]{pointer} = "diagonal";
				}
				else {
					$matrix[$i][$j]{score}   = $left_score;
					$matrix[$i][$j]{pointer} = "left";
				}
			} else {
				if ($up_score >= $left_score) {
					$matrix[$i][$j]{score}   = $up_score;
					$matrix[$i][$j]{pointer} = "up";
				}
				else {
					$matrix[$i][$j]{score}   = $left_score;
					$matrix[$i][$j]{pointer} = "left";
				}
			}
			
			# set maximum score
			if ($matrix[$i][$j]{score} >= $max_score && $i >= length($miRNA) + 6) {
				$max_i     = $i;
				$max_j     = $j;
				#push (@max_i, $i);
				#push (@max_j, $j);
				$max_score = $matrix[$i][$j]{score};
				#push (@max_score, $matrix[$i][$j]{score});
			}
		}
	}
	#print "MAX SCORE: $max_score\n";
	# return if ($max_score < 40 );

	# trace-back
	my (@previous_mfe, @last_align, @feature, @conseq);
	#my $align_num = 0;
	#print $chr_number[$pair_number]."\t".$chr_start[$pair_number]."\n";
	# my @eachline = split /\n/, $conseq[$pair_number+1];
	#print $_."\n" foreach (@eachline);
	$pair_number++;
	#my $MFE = -5;
	
	#for (my $aa = $#max_i-1, $bb = $#max_j-1, $cc = 0; $aa > 0 && $bb > 0; $aa--, $bb--, $cc++) {
		my $align1 = "";
		my $align2 = "";
		my $match_letter = "";
		my ($repeat, $seed_scaling) = (0, 0);
		my ($seq1_start, $seq2_start);
		my $j = $max_j; # the seq2's position of the highest score
		my $i = $max_i; # the seq1's position of the highest score
		my $key = 0;
		my $miRNA_leng = 0;
		my $UTR_leng = 0;
		my ($GC_content, $AU_content) = (0, 0);
		my ($mismatch_num, $gap_num, $wobble_num, $WC_num) = (0, 0, 0, 0);
		my $align2_letter;	
		while (1) {
			last if $matrix[$i][$j]{pointer} eq "none";
		
			if ($matrix[$i][$j]{pointer} eq "diagonal") {
				$align1 .= substr($miRNA, $j-1, 1);
				$align2 .= substr($UTRseq_part, $i-1, 1);
				$align2_letter = substr($UTRseq_part, $i-1, 1);
				if ((substr($miRNA, $j-1, 1) =~ /$align2_letter/) || (substr($miRNA, $j-1, 1) =~ /A/ && substr($UTRseq_part, $i-1, 1) =~ /C|G/) || (substr($miRNA, $j-1, 1) =~ /C|G/ && substr($UTRseq_part, $i-1, 1) =~ /A/) || (substr($miRNA, $j-1, 1) =~ /T|U/ && substr($UTRseq_part, $i-1, 1) =~ /C/) || (substr($miRNA, $j-1, 1) =~ /C/ && substr($UTRseq_part, $i-1, 1) =~ /T|U/) || (substr($miRNA, $j-1, 1) =~ /U/ && substr($UTRseq_part, $i-1, 1) =~ /T/)) {
					$match_letter .= " "; # pass mismatch
					$mismatch_num++;
				} elsif ((substr($miRNA, $j-1, 1) =~ /G/ && substr($UTRseq_part, $i-1, 1) =~ /T|U/) || (substr($miRNA, $j-1, 1) =~ /T|U/ && substr($UTRseq_part, $i-1, 1) =~ /G/)) {
					$match_letter .= ":"; # mark wobble match
					$wobble_num++;
				} else {
					$match_letter .= "|"; # mark match
					if (substr($UTRseq_part, $i-1, 1) =~ /G|C/) {
						$GC_content++;
					} elsif (substr($UTRseq_part, $i-1, 1) =~ /A|T|U/) {
						$AU_content++;
					}
				}
				if ($key =~ 0) {
					$seq1_start = length($miRNA) - $j + 1;
					#$seq2_start += (length($UTRseq_part) - ($i - 1));
					#print $seq2_start."\n";
				}
				if ($j-1 < 8) {
					$seed_scaling++;
				}
				$key = 1; #close the key of the last position of sequence alignment
				$miRNA_leng++;
				$UTR_leng++;
				$WC_num++;
				$i--; $j--;
			}
			elsif ($matrix[$i][$j]{pointer} eq "left") {
				$align1 .= substr($miRNA, $j-1, 1);
				$align2 .= "-";
				$match_letter .= " ";
				$gap_num++;	$miRNA_leng++;
				$j--;
			}
			elsif ($matrix[$i][$j]{pointer} eq "up") {
				$align1 .= "-";
				$align2 .= substr($UTRseq_part, $i-1, 1);
				$match_letter .= " ";
				$gap_num++;	$UTR_leng++;
				$i--;
			}
		}
		#next if ($miRNA_leng < 15||$mismatch_num>4||$wobble_num>6||$mismatch_num+$wobble_num>7);
		# return if ($miRNA_leng < 16||$mismatch_num>4||$wobble_num>6||$mismatch_num+$wobble_num>7);
		$align1 = reverse $align1;
		$align2 = reverse $align2;
		$match_letter = reverse $match_letter;
		
		#Delete the tail of gap at miRNA align
		if ($align1 =~ /(-+)$/) {
			$UTR_leng = $UTR_leng - length($&);
			#print "[$cc]\t$align1";
			$align1 = substr($align1, 0, length($align1) - length($&));
			$align2 = substr($align2, 0, length($align2) - length($&));
			$match_letter = substr($match_letter, 0, length($match_letter) - length($&));
			$gap_num = $gap_num - length($&);
			#$max_score = $max_score - $GAP_OPENNING - $GAP_EXTENSION * (length($&) - 1);
			#print "\t=>$align1\t\tScore:$max_score\n";
		}
		#next if($gap_num >4);
=non-used
		for (my $count = 0; $count <= $#last_align; $count++) {
			if ($align1 =~ /^$last_align[$count]{s1}$/ && $align2 =~ /^$last_align[$count]{s2}$/) {
				#print "REPEAT!!\t$cc<->$count\t$align1\t$align2\tScore:$max_score[$cc]\n";
				$repeat = 1;
				last;
			}
		}
		#$last_align[$cc]{s1} = "$align1";
		#$last_align[$cc]{s2} = "$align2";
		next if ($repeat == 1); # filter the repeated align after deleting the tail of gap at miRNA align
=cut	
		#Scaling seed region
		#if ($seq1_start =~ 1) {
		#	for (my $j = 1; $j <= 8 - $seq1_start; $j++) {
		#		if (substr($align1, $j, 1) =~ /[A-Za-z]/ && substr($match_letter, $j, 1) =~ /\S/) {
		#			$seed_scaling++;
		#		}
		#	}
		#} else {
		#	for (my $j = 0; $j <= 8 - $seq1_start; $j++) {
		#		if (substr($align1, $j, 1) =~ /[A-Za-z]/ && substr($match_letter, $j, 1) =~ /\S/) {
		#			$seed_scaling++;
		#		}
		#	}
		#}
#		system "RNAup -b -d2 --noLP -c 'S' < $RESULT_OUTPUT_PATHWAY.seq > RESULT_OUTPUT_PATHWAY.rnaup";
		#Minimum engery scores
		my ($structure, $energy, $filter) = (0, 0, 0);	
		my $align2_tmp = $align2;
		my $cc;
		$align2_tmp =~ s/-//g;
		my $align1_rev = reverse $align1;
		$align1_rev =~ s/-//g;
		open (FHDW2, ">$RESULT_OUTPUT_PATHWAY.seq") || die "$!\n";
		print FHDW2 ">align_3'UTR\t$align2\n$align2_tmp\n>align_miRNA\t$align1\n$align1_rev\n\n";
#		system "C:\\AppServ\\www\\ViennaRNA-1.8.5_win\\Progs\\RNAduplex.exe < $RESULT_OUTPUT_PATHWAY.seq > $RESULT_OUTPUT_PATHWAY.duplex";
#		system "/media/data/nfs/home/kang/public_html/miRNAServer/HumanMain/Tool/ViennaRNA-1.8.5_win/Progs/RNAduplex < $RESULT_OUTPUT_PATHWAY.seq > $RESULT_OUTPUT_PATHWAY.duplex";
#		system "C:\\AppServ\\www\\PlantServer\\Tool\\ViennaRNA-1.8.5_win\\Progs\\RNAduplex.exe < $RESULT_OUTPUT_PATHWAY.seq > $RESULT_OUTPUT_PATHWAY.duplex";
		# &end("sw");
		# &start;
		system "RNAduplex < $RESULT_OUTPUT_PATHWAY.seq > $RESULT_OUTPUT_PATHWAY.duplex";
		# &end("RNAduplex");
#		system "/nfs/master/kang/Plant/Tool/ViennaRNA-1.8.5/Progs/RNAduplex < $RESULT_OUTPUT_PATHWAY.seq > $RESULT_OUTPUT_PATHWAY.duplex";		
#		system "/media/data/nfs/home/kang/PlantMain/Tool/ViennaRNA-1.8.5/Progs/RNAduplex < $RESULT_OUTPUT_PATHWAY.seq > $RESULT_OUTPUT_PATHWAY.duplex";


		
		# &start;
		open (FHDR3, "$RESULT_OUTPUT_PATHWAY.duplex") || die "$!\n";
		while (my $line=<FHDR3>) {
			chomp $line;
			if ($line =~ /^>|\n/) {
				#print $line."\n";
				;
			} else {
				if ($line =~ /^(.+|\(+|&|\)+)/) {
					my $tmd = $1;
					$tmd =~ s/://g;
					$tmd =~ s/\s+/ /g;
					my @t = split / /, $tmd;
					#print "$t[0]\t$t[1]\t$t[2]\t$t[3]\n";
					$structure = $t[0];
					my @ener_pos_1 = split /,/, $t[1];
					my $ener_leng_1 = $ener_pos_1[1] - $ener_pos_1[0];
					my @ener_pos_2 = split /,/, $t[2];
					my $ener_leng_2 = $ener_pos_2[1] - $ener_pos_2[0];
					$t[3] =~ s/\(|\)//g;
					$energy = $t[3];
					$energy *= ($ener_leng_1 / length($align1_rev)) * ($ener_leng_2 / length($align2_tmp)) if ($energy =~ /\d/);
				}
			}
		}
		close (FHDR3);
		# return if ($energy eq "" or $energy > -5); #! next/return problem
		#next if ($energy > -5);
		push (@previous_mfe, $energy);
		#my ($UTR_struc, $miRNA_struc) =~ split /&/, $structure;
		#$miRNA_struc = reverse $miRNA_struc;
		#$miRNA_struc =~ s/\)/\(/g;
		
		#print "UTR length:".$UTR_leng."\n";
		my (@starts, @ends);
		#$seq1_start = $max_j - $miRNA_leng + 1;
		$seq2_start = $align_start - $UTR_leng + 8;
		my $seq1_end = $seq1_start + $miRNA_leng - 1;
		my $seq2_end = $seq2_start + $UTR_leng - 1;
		push (@starts, $seq2_start);
		push (@ends, $seq2_end);	
		
		#Exclude repeated alignment
		for (my $i = 0; $i <= $#starts; $i++) {
			if (abs ($seq2_start - $starts[$i]) <= 10 && abs ($seq2_end - $ends[$i]) <= 10) {
				last if ($i == $#previous_mfe);
				$repeat = 1 if ($energy >= $previous_mfe[$i]);
			}
		}
		# return if ($repeat == 1);
		
		#Distance to the nearest target site
		#for (my $count = $#starts; $count > 0; $count--) {
		#	my $diff = abs $seq2_start - $starts[$count-1];
		#	if ($diff < $distance) {
		#		$distance = $diff;
		#	}
		#}
		#my	$dis_percent = $distance/length($UTRseq);
		my $nearer_end = ($seq2_start > (length($UTRseq) - $seq2_start)) ? (length($UTRseq) - $seq2_start) : $seq2_start;
		# return if ($nearer_end <= 15);
		my $pair_num = $wobble_num + $GC_content + $AU_content;
		my $seq1_pair_proportion = $pair_num / $miRNA_leng;
		my $seq2_pair_proportion = $pair_num / $UTR_leng;

		#3 binding
		my $three_binding = 0;
		$align1_rev = reverse $align1;
		my $match_letter_rev = reverse $match_letter;
		my $down = $seq1_start + $miRNA_leng - 1;
		for (my $temp_seq1_end = length($align1_rev)-1; $down > 12 && $temp_seq1_end > 0; $temp_seq1_end--) {
			if (substr($align1_rev, $temp_seq1_end, 1) =~ /A|T|C|G|U/) {
				if ($down >= 13 && $down <= 17 && substr($match_letter_rev, $temp_seq1_end, 1) =~ /\|/) {
					$three_binding++;
				}
				$down--;
			}
		}
		
		#central 9-11
		my $central_binding = 0;
		$align1_rev = reverse $align1;
		$match_letter_rev = reverse $match_letter;
		$down = $seq1_start + $miRNA_leng - 1;
		for (my $temp_seq1_end = length($align1_rev)-1; $down > 9 && $temp_seq1_end > 0; $temp_seq1_end--) {
			if (substr($align1_rev, $temp_seq1_end, 1) =~ /A|T|C|G|U/) {
				if ($down >= 9 && $down <= 11 && substr($match_letter_rev, $temp_seq1_end, 1) =~ /\|/) {
					$central_binding++;
				}
				$down--;
			}
		}		
		
		#Determine 6mer match
		my ($six_mer, $noncanonical, $exclude) = (1, 0, 0);
		if ($seq1_start == 1) {
			for (my $j = 1; $j <= 7 - $seq1_start; $j++) {
				if (substr($align1_rev, $j, 1) =~ /[A-Za-z]/ && substr($match_letter_rev, $j, 1) =~ /\|/ && $six_mer == 1) {
					$six_mer = 1;
				} else {
					#only allow one GU or mismatch in seed region
					if (substr($align1_rev, $j, 1) =~ /[A-Z]/ && substr($match_letter_rev, $j, 1) =~ /:|\s/ && $exclude == 0) {
						$six_mer = 0;
						$exclude = 1;
					} else {
						$six_mer = 0;
						last;
					}
				}
			}
		} else {
			for (my $j = 0; $j <= 7 - $seq1_start; $j++) {
				if (substr($align1_rev, $j, 1) =~ /[A-Za-z]/ && substr($match_letter_rev, $j, 1) =~ /\|/ && $six_mer == 1) {
					$six_mer = 1;
				} else {
					if (substr($align1_rev, $j, 1) =~ /[A-Z]/ && substr($match_letter_rev, $j, 1) =~ /:|\s/ && $noncanonical == 0) {
						$six_mer = 0;
						$noncanonical = 1;
					} elsif (substr($align1_rev, $j, 1) =~ /[A-Z]/ && substr($match_letter_rev, $j, 1) =~ /:|\s/ && $noncanonical == 1) {
						$exclude = 1;
						last;
					} else {
						$six_mer = 0;
						$exclude = 1;
						last;
					}
				}
			}
		}
		#next if ($exclude == 1);
		
		#Determine 7mer-1A, 7mer-m8, 8mer
		my ($seven_mer_1A, $seven_mer_m8, $eight_mer) = (0, 0, 0);
		if ($six_mer == 1) {
			if ($seq1_start == 1 && substr($align1_rev, 0, 1) =~ /[A-Z]/ && substr($match_letter_rev, 0, 1) =~ /\|/ && substr($align2, length($align2) - 1, 1) =~ /A/) {
				$seven_mer_1A = 1;
			} else {
				$seven_mer_1A = 0;
			}
			if (substr($align1_rev, 8 - $seq1_start, 1) =~ /[A-Za-z]/ && substr($match_letter_rev, 8 - $seq1_start, 1) =~ /\|/) {
				$seven_mer_m8 = 1;
				if ($seq1_start == 1 && substr($align1_rev, 0, 1) =~ /[A-Z]/ && substr($match_letter_rev, 0, 1) =~ /\|/ && substr($align2, length($align2) - 1, 1) =~ /A/) {
					$eight_mer = 1;
				} else {
					$eight_mer = 0;
				}
			} else {
				$seven_mer_m8 = 0;
			}
		}
#		next if ($six_mer == 0 && $seven_mer_1A == 0 && $seven_mer_m8 == 0 && $eight_mer == 0 && $noncanonical == 0);
		#$filter = 1 if ($six_mer == 0 && $seven_mer_1A == 0 && $seven_mer_m8 == 0 && $eight_mer == 0 && $noncanonical == 0);
		
		#Count Watson-Crick pairing number
		#my $WC_pair = 0;
		#for (my $i = 0; $i < length($match_letter); $i++) {
		#	if (substr($match_letter, $i, 1) =~ /\|/) {
		#		$WC_pair++;
		#	}
		#}
		my $compesatory = 0;
		my ($WC_pair_six, $WC_pair_m8, $WC_pair_1A, $WC_pair_eight, $WC_pair_noncanonical) = (0, 0, 0, 0, 0);
		if ($six_mer == 1) {
			$WC_pair_six = 1 if ($three_binding >= 1);
		}
		if ($seven_mer_m8 == 1) {
			$WC_pair_m8 = 1 if ($three_binding >= 1);
		}
		if ($seven_mer_1A == 1) {
			$WC_pair_1A = 1 if ($three_binding >= 1);
		}
		if ($eight_mer == 1) {
			$WC_pair_eight = 1 if ($three_binding >= 1);
		}
		if ($compesatory == 1) {
			$WC_pair_noncanonical = 1 if ($three_binding >= 1);
		}
		
		my $local_seed_end = 7 - $seq1_start;
		my $global_seed_end = 7;
		my $seed_lnth = 6;
		if ($eight_mer == 1 || $seven_mer_m8 == 1) {
			$local_seed_end += 1;
			$global_seed_end += 1;
			$seed_lnth += 1;
			$align_start -= 1;
		}
		# print "$local_seed_end $global_seed_end $seed_lnth $align_start"; die;
		
		#Flanking AU composition scores
		#my $UTRseq_rev = reverse $UTRseq;
		my ($AU_composition_num, $AU_composition_score) = (0, 0);
		my ($Up_AU_composition_num, $Up_AU_score, $Down_AU_composition_num, $Down_AU_score) = (0, 0, 0, 0);
		for (my $i = 1; $i <= 30; $i++) {
			if (substr($UTRseq, $align_start -32 + $i, 1) =~ /A|T|U/) {
				$AU_composition_num++;
				$Up_AU_composition_num++;
				$AU_composition_score += 1/$i;
				$Up_AU_score += 1/$i;
			}
			if (($align_start + $seed_lnth + $i - 2) < length($UTRseq) and substr($UTRseq, $align_start + $seed_lnth + $i - 2, 1) =~ /A|T|U/) {
				$AU_composition_num++;
				$Down_AU_composition_num++;
				$AU_composition_score += 1/$i;
				$Down_AU_score += 1/$i;
			}
		}
		for (my $i = 1; $i < $seed_lnth; $i++) {
			if (substr($UTRseq, $align_start + $i - 1, 1) =~ /A|T|U/) {
				$AU_composition_num++;
				$AU_composition_score += 1/$i;
			}
		}
		my $flanking_AU_percent = $AU_composition_num/70;
		my $Up_flanking_AU_percent = $Up_AU_composition_num/30;
		my $Down_flanking_AU_percent = $Down_AU_composition_num/30;
		my ($flanking_AU_six, $flanking_AU_m8, $flanking_AU_1A, $flanking_AU_eight, $flanking_AU_noncanonical) = (0, 0, 0, 0, 0);
		my $Up_flanking_AU;
		my $Down_flanking_AU;
		if ($flanking_AU_percent >= 0.6) {
			$flanking_AU_six = 1 if ($six_mer == 1);
			$flanking_AU_m8 = 1 if ($seven_mer_m8 == 1);
			$flanking_AU_1A = 1 if ($seven_mer_1A == 1);
			$flanking_AU_eight = 1 if ($eight_mer == 1);
			$flanking_AU_noncanonical = 1 if ($noncanonical == 1);
		} else {
			$flanking_AU_six = 0 if ($six_mer == 1);
			$flanking_AU_m8 = 0 if ($seven_mer_m8 == 1);
			$flanking_AU_1A = 0 if ($seven_mer_1A == 1);
			$flanking_AU_eight = 0 if ($eight_mer == 1);
			$flanking_AU_noncanonical = 0 if ($noncanonical == 1);
		}
		if ($Up_flanking_AU_percent >= 0.6) {
			$Up_flanking_AU = 1;
		} else {
			$Up_flanking_AU = 0;
		}
		if ($Down_flanking_AU_percent >= 0.6) {
			$Down_flanking_AU = 1;
		} else {
			$Down_flanking_AU = 0;
		}
		# &end("RNA");
		
		#Frequency of single nuclei at seed matching site
		# &start;
		my @single_nuclei = (0, 0, 0, 0);
		my $align2_rev = reverse $align2;
		my $single_seed_end = $local_seed_end;
		my $nt2;
		if ($seq1_start == 1) {
			for (my $k = 1; $k <= $single_seed_end; $k++) {	 
				if (substr($align1_rev, $k, 1) =~ /[A-Za-z]/ && substr($match_letter_rev, $k, 1) =~ /\S/) {
					$nt2 = substr($align2_rev, $k, 1);
					$single_nuclei[0]++ if ($nt2 =~ /A/);
					$single_nuclei[1]++ if ($nt2 =~ /T|U/);
					$single_nuclei[2]++ if ($nt2 =~ /C/);
					$single_nuclei[3]++ if ($nt2 =~ /G/);
				} elsif (substr($align1_rev, $k, 1) =~ /-/) {
					$single_seed_end++;
				}
			}
		} else {
			for (my $k = 0; $k <= $single_seed_end; $k++) {	
				if (substr($align1_rev, $k, 1) =~ /[A-Za-z]/ && substr($match_letter_rev, $k, 1) =~ /\S/) {
					$nt2 = substr($align2_rev, $k, 1);
					$single_nuclei[0]++ if ($nt2 =~ /A/);
					$single_nuclei[1]++ if ($nt2 =~ /T|U/);
					$single_nuclei[2]++ if ($nt2 =~ /C/);
					$single_nuclei[3]++ if ($nt2 =~ /G/);
				} elsif (substr($align1_rev, $k, 1) =~ /-/) {
					$single_seed_end++;
				}
			} 
		}
		my $single_nuclei_sum = $single_nuclei[0] + $single_nuclei[1] + $single_nuclei[2] + $single_nuclei[3];
		$single_nuclei_sum == 0 and $single_nuclei_sum = 1;
		my $single_A = $single_nuclei[0] / $single_nuclei_sum;
		my $single_T = $single_nuclei[1] / $single_nuclei_sum;
		my $single_C = $single_nuclei[2] / $single_nuclei_sum;
		my $single_G = $single_nuclei[3] / $single_nuclei_sum;
		
###		#Frequency of Di-nuclei at seed matching site
		my @di_nuclei = (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
		my $di_nuclei_sum = 0;
		my $di_seed_end = $local_seed_end;
		my $di_nt2;
		if ($seq1_start == 1) {
			for (my $k = 1; $k < $di_seed_end; $k++) {	 
				if (substr($align1_rev, $k, 2) =~ /[A-Za-z]{2}/  && substr($match_letter_rev, $k, 2) =~ /\S\S/) {
					$di_nt2 = substr($align2_rev, $k, 2);
					$di_nuclei[0]++ if ( $di_nt2 =~ /AA/);
					$di_nuclei[1]++ if ( $di_nt2 =~ /AT/);
					$di_nuclei[2]++ if ( $di_nt2 =~ /AC/);
					$di_nuclei[3]++ if ( $di_nt2 =~ /AG/);
					$di_nuclei[4]++ if ( $di_nt2 =~ /TA/);
					$di_nuclei[5]++ if ( $di_nt2 =~ /TT/);
					$di_nuclei[6]++ if ( $di_nt2 =~ /TC/);
					$di_nuclei[7]++ if ( $di_nt2 =~ /TG/);
					$di_nuclei[8]++ if ( $di_nt2 =~ /CA/);
					$di_nuclei[9]++ if ( $di_nt2 =~ /CT/);
					$di_nuclei[10]++ if( $di_nt2 =~ /CC/);
					$di_nuclei[11]++ if( $di_nt2 =~ /CG/);
					$di_nuclei[12]++ if( $di_nt2 =~ /GA/);
					$di_nuclei[13]++ if( $di_nt2 =~ /GT/);
					$di_nuclei[14]++ if( $di_nt2 =~ /GC/);
					$di_nuclei[15]++ if( $di_nt2 =~ /GG/);
				} elsif (substr($align1_rev, $k, 2) =~ /(-+)/) {
					$di_seed_end += length($1);
				}
			}
		} else {
			for (my $k = 0; $k < $di_seed_end; $k++) {
				if (substr($align1_rev, $k, 2) =~ /[A-Za-z]{2}/  && substr($match_letter_rev, $k, 2) =~ /\S\S/) {
					$di_nt2 = substr($align2_rev, $k, 2);
					$di_nuclei[$di2num{$di_nt2}]++;
				} elsif (substr($align1_rev, $k, 2) =~ /(-+)/) {
					$di_seed_end += length($1);
				}
			} 
		}
		for (my $i = 0; $i < $#di_nuclei; $i++) {
			$di_nuclei_sum += $di_nuclei[$i];
		}
		# return if ($di_nuclei_sum == 0);
		$di_nuclei_sum == 0 and $di_nuclei_sum = 1;
		my $di_AA = $di_nuclei[0]/$di_nuclei_sum;
		my $di_AT = $di_nuclei[1]/$di_nuclei_sum;
		my $di_AC = $di_nuclei[2]/$di_nuclei_sum;
		my $di_AG = $di_nuclei[3]/$di_nuclei_sum;
		my $di_TA = $di_nuclei[4]/$di_nuclei_sum;
		my $di_TT = $di_nuclei[5]/$di_nuclei_sum;
		my $di_TC = $di_nuclei[6]/$di_nuclei_sum;
		my $di_TG = $di_nuclei[7]/$di_nuclei_sum;
		my $di_CA = $di_nuclei[8]/$di_nuclei_sum;
		my $di_CT = $di_nuclei[9]/$di_nuclei_sum;
		my $di_CC = $di_nuclei[10]/$di_nuclei_sum;
		my $di_CG = $di_nuclei[11]/$di_nuclei_sum;
		my $di_GA = $di_nuclei[12]/$di_nuclei_sum;
		my $di_GT = $di_nuclei[13]/$di_nuclei_sum;
		my $di_GC = $di_nuclei[14]/$di_nuclei_sum;
		my $di_GG = $di_nuclei[15]/$di_nuclei_sum;
		
		#Frequency of single nuclei at out-seed flanking site
		my ($single_A_out, $single_T_out, $single_C_out, $single_G_out) = (0, 0, 0, 0);
		my $single_nt_out;
		for (my $i = 1; $i <= 70 and $seq2_start - 32 + $i < length($UTRseq); $i++) {
			$single_nt_out = substr($UTRseq, $seq2_start -32 + $i, 1);
			$single_A_out++ if ( $single_nt_out =~ /A/);
			$single_T_out++ if ( $single_nt_out =~ /U|T/);
			$single_C_out++ if ( $single_nt_out =~ /C/);
			$single_G_out++ if ( $single_nt_out =~ /G/);
		}
		
###		#Frequency of di-single nucleis at out-seed flanking site
		my @di_nuclei_out = (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
		my $di_nuclei_out_sum = 0;
		my $di_nt_out;
		for (my $i = 1; $i < 70 and $seq2_start - 32 + $i <= length($UTRseq) - 2; $i++) {
			$di_nt_out = substr($UTRseq, $seq2_start -32 + $i, 2);
			exists $di2num{$di_nt_out} and $di_nuclei_out[$di2num{$di_nt_out}]++;
		}
		for (my $i = 0; $i <= $#di_nuclei_out; $i++) {
			$di_nuclei_out_sum += $di_nuclei_out[$i];
		}
		next if ($di_nuclei_out_sum == 0);
		my $di_AA_out = $di_nuclei_out[0]/$di_nuclei_out_sum;
		my $di_AT_out = $di_nuclei_out[1]/$di_nuclei_out_sum;
		my $di_AC_out = $di_nuclei_out[2]/$di_nuclei_out_sum;
		my $di_AG_out = $di_nuclei_out[3]/$di_nuclei_out_sum;
		my $di_TA_out = $di_nuclei_out[4]/$di_nuclei_out_sum;
		my $di_TT_out = $di_nuclei_out[5]/$di_nuclei_out_sum;
		my $di_TC_out = $di_nuclei_out[6]/$di_nuclei_out_sum;
		my $di_TG_out = $di_nuclei_out[7]/$di_nuclei_out_sum;
		my $di_CA_out = $di_nuclei_out[8]/$di_nuclei_out_sum;
		my $di_CT_out = $di_nuclei_out[9]/$di_nuclei_out_sum;
		my $di_CC_out = $di_nuclei_out[10]/$di_nuclei_out_sum;
		my $di_CG_out = $di_nuclei_out[11]/$di_nuclei_out_sum;
		my $di_GA_out = $di_nuclei_out[12]/$di_nuclei_out_sum;
		my $di_GT_out = $di_nuclei_out[13]/$di_nuclei_out_sum;
		my $di_GC_out = $di_nuclei_out[14]/$di_nuclei_out_sum;
		my $di_GG_out = $di_nuclei_out[15]/$di_nuclei_out_sum;
###################################################################################################################################################
		#Frequency of Tri-nuclei at seed matching site
		#Frequency of Di-nuclei at seed matching site
		my @tri_nuclei = (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
		my $tri_nuclei_sum = 0;
		my $dash = 0;
		#my $tri_seed_end = $local_seed_end;
		my $tri_nt;
		if ($seq1_start==1){
			for (my $k = 1; $k < length($align1_rev)+1 + $dash && $k <= length($align1_rev) - 3; $k++) { #! prevent out bound	 
				if (substr($align1_rev, $k, 3) =~ /[A-Za-z]{3}/  && substr($match_letter_rev, $k, 3) =~ /\S\S\S/) {
					$tri_nt = substr($align2_rev, $k, 3);
					$tri_nuclei[$tri2num{$tri_nt}]++;
				} elsif (substr($align1_rev, $k, 3) =~ /(-+)/) {
					$dash += length($1);
				}
			}
		} else {
			for (my $k = 0; $k < length($align1_rev)+$dash && $k <= length($align1_rev) - 3 ; $k++) { #! prevent out bound	 
				if (substr($align1_rev, $k, 3) =~ /[A-Za-z]{3}/  && substr($match_letter_rev, $k, 3) =~ /\S\S\S/) {
					$tri_nt = substr($align2_rev, $k, 3);
					$tri_nuclei[$tri2num{$tri_nt}]++;
				} elsif (substr($align1_rev, $k, 3) =~ /(-+)/) {
					$dash += length($1);
				}		
			}
		}
		for (my $i = 0; $i < 64; $i++) {
			$tri_nuclei_sum += $tri_nuclei[$i];
		}
		($tri_nuclei_sum == 0) and $tri_nuclei_sum = 1;
		my @tri_nuclei_frequency = (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
		for (my $tri=0; $tri < 64; $tri++){
			$tri_nuclei_frequency[$tri] = $tri_nuclei[$tri]/$tri_nuclei_sum;
		}
		#print "\nalign2_rev: $align2_rev\n";
		#print "\nmatchl_rev: $match_letter_rev\n";
		#print "\nalign1_rev: $align1_rev\n";
		#print "\nTrinucleiSUM: $tri_nuclei_sum\n";
		#print "tri_nuclei_frequency: @tri_nuclei_frequency\n\n";
		my $tricountlength = length($align1_rev)+1;
		#print "-----------$tricountlength-------------";
		#@tri_frequency;
		#print "OUT1: @tri_nuclei_frequency\n\n";
		#Frequency of tri-single nucleis at out-seed flanking site		
		my @tri_nuclei_out = (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
		my $tri_nt_out;
		for (my $i = 1; $i < 70 and $seq2_start - 32 + $i <= length($UTRseq) - 3; $i++) {	#1 to 70
			$tri_nt_out = substr($UTRseq, $seq2_start -32 + $i, 3);
            exists $tri2num{$tri_nt_out} and $tri_nuclei_out[$tri2num{$tri_nt_out}]++;
		}		
		my $tri_nuclei_out_sum;
		for (my $i = 0; $i < 64; $i++) {
			$tri_nuclei_out_sum += $tri_nuclei_out[$i];
		}
		($tri_nuclei_out_sum == 0) and $tri_nuclei_out_sum = 1;
		my @tri_nuclei_out_frequency = (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
		for (my $tri=0; $tri < 64; $tri++){
			$tri_nuclei_out_frequency[$tri] = $tri_nuclei_out[$tri]/$tri_nuclei_out_sum;
		}
###################################################################################################################################################		
		#Frequency of 6 types of base pairing at seed matching site
		my @base_pair_type = (0, 0, 0, 0, 0, 0);
		my $pair_seed_end = $local_seed_end;
		if ($seq1_start == 1) {
			for (my $k = 1; $k <= $pair_seed_end; $k++) {	 
				if (substr($align1_rev, $k, 1) =~ /[A-Za-z]/ && substr($match_letter_rev, $k, 1) =~ /\S/) {
					$base_pair_type[0]++ if (substr($align1_rev, $k, 1) =~ /A/ && substr($align2_rev, $k, 1) =~ /U|T/);
					$base_pair_type[1]++ if (substr($align1_rev, $k, 1) =~ /T|U/ && substr($align2_rev, $k, 1) =~ /A/);
					$base_pair_type[2]++ if (substr($align1_rev, $k, 1) =~ /G/ && substr($align2_rev, $k, 1) =~ /C/);
					$base_pair_type[3]++ if (substr($align1_rev, $k, 1) =~ /C/ && substr($align2_rev, $k, 1) =~ /G/);
					$base_pair_type[4]++ if (substr($align1_rev, $k, 1) =~ /G/ && substr($align2_rev, $k, 1) =~ /U|T/);
					$base_pair_type[5]++ if (substr($align1_rev, $k, 1) =~ /U|T/ && substr($align2_rev, $k, 1) =~ /G/);
				} elsif (substr($align1_rev, $k, 1) =~ /-/) {
					$pair_seed_end++;
				}
			}
		} else {
			for (my $k = 0; $k <= $pair_seed_end; $k++) {	
				if (substr($align1_rev, $k, 1) =~ /[A-Za-z]/ && substr($match_letter_rev, $k, 1) =~ /\S/) {
					$base_pair_type[0]++ if (substr($align1_rev, $k, 1) =~ /A/ && substr($align2_rev, $k, 1) =~ /U|T/);
					$base_pair_type[1]++ if (substr($align1_rev, $k, 1) =~ /T|U/ && substr($align2_rev, $k, 1) =~ /A/);
					$base_pair_type[2]++ if (substr($align1_rev, $k, 1) =~ /G/ && substr($align2_rev, $k, 1) =~ /C/);
					$base_pair_type[3]++ if (substr($align1_rev, $k, 1) =~ /C/ && substr($align2_rev, $k, 1) =~ /G/);
					$base_pair_type[4]++ if (substr($align1_rev, $k, 1) =~ /G/ && substr($align2_rev, $k, 1) =~ /U|T/);
					$base_pair_type[5]++ if (substr($align1_rev, $k, 1) =~ /U|T/ && substr($align2_rev, $k, 1) =~ /G/);
				} elsif (substr($align1_rev, $k, 1) =~ /-/) {
					$pair_seed_end++;
				}
			} 
		}
		my $base_pair_type_sum = $base_pair_type[0] + $base_pair_type[1] + $base_pair_type[2] + $base_pair_type[3] + $base_pair_type[4] + $base_pair_type[5];
		$base_pair_type_sum == 0 and $base_pair_type_sum = 1;
		my $pair_AU = $base_pair_type[0] / $base_pair_type_sum;
		my $pair_UA = $base_pair_type[1] / $base_pair_type_sum;
		my $pair_GC = $base_pair_type[2] / $base_pair_type_sum;
		my $pair_CG = $base_pair_type[3] / $base_pair_type_sum;
		my $pair_GU = $base_pair_type[4] / $base_pair_type_sum;
		my $pair_UG = $base_pair_type[5] / $base_pair_type_sum;
		
		#Frequency of 32 types of bi-nucleis base pairing at seed matching site
		my @bi_base_pair = (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
		my ($bi_seed_end, $bi_base_pair_sum) = ($local_seed_end, 0);
		if ($seq1_start == 1) {
			for (my $k = 1; $k < $bi_seed_end; $k++) {	 
				if (substr($align1_rev, $k, 2) =~ /[A-Za-z]{2}/  && substr($match_letter_rev, $k, 2) =~ /\S\S/) {
					$bi_base_pair[0]++ if (substr($align1_rev, $k, 2) =~ /AA/ && substr($align2_rev, $k, 2) =~ /UU|TT/);
					$bi_base_pair[1]++ if (substr($align1_rev, $k, 2) =~ /AT/ && substr($align2_rev, $k, 2) =~ /UA|TA/);
					$bi_base_pair[2]++ if (substr($align1_rev, $k, 2) =~ /AT/ && substr($align2_rev, $k, 2) =~ /UG|TG/);
					$bi_base_pair[3]++ if (substr($align1_rev, $k, 2) =~ /AC/ && substr($align2_rev, $k, 2) =~ /UG|TG/);
					$bi_base_pair[4]++ if (substr($align1_rev, $k, 2) =~ /AG/ && substr($align2_rev, $k, 2) =~ /UC|TC/);
					$bi_base_pair[5]++ if (substr($align1_rev, $k, 2) =~ /AG/ && substr($align2_rev, $k, 2) =~ /UU|TT/);
					$bi_base_pair[6]++ if (substr($align1_rev, $k, 2) =~ /TA/ && substr($align2_rev, $k, 2) =~ /AU|AT/);
					$bi_base_pair[7]++ if (substr($align1_rev, $k, 2) =~ /TA/ && substr($align2_rev, $k, 2) =~ /GU|GT/);
					$bi_base_pair[8]++ if (substr($align1_rev, $k, 2) =~ /TT/ && substr($align2_rev, $k, 2) =~ /AA/);
					$bi_base_pair[9]++ if (substr($align1_rev, $k, 2) =~ /TT/ && substr($align2_rev, $k, 2) =~ /GA/);
					$bi_base_pair[10]++ if (substr($align1_rev, $k, 2) =~ /TT/ && substr($align2_rev, $k, 2) =~ /AG/);
					$bi_base_pair[11]++ if (substr($align1_rev, $k, 2) =~ /TC/ && substr($align2_rev, $k, 2) =~ /AG/);
					$bi_base_pair[12]++ if (substr($align1_rev, $k, 2) =~ /TC/ && substr($align2_rev, $k, 2) =~ /GG/);
					$bi_base_pair[13]++ if (substr($align1_rev, $k, 2) =~ /TG/ && substr($align2_rev, $k, 2) =~ /AC/);
					$bi_base_pair[14]++ if (substr($align1_rev, $k, 2) =~ /TG/ && substr($align2_rev, $k, 2) =~ /GC/);
					$bi_base_pair[15]++ if (substr($align1_rev, $k, 2) =~ /TG/ && substr($align2_rev, $k, 2) =~ /AU|AT/);
					$bi_base_pair[16]++ if (substr($align1_rev, $k, 2) =~ /CA/ && substr($align2_rev, $k, 2) =~ /GU|GT/);
					$bi_base_pair[17]++ if (substr($align1_rev, $k, 2) =~ /CT/ && substr($align2_rev, $k, 2) =~ /GA/);
					$bi_base_pair[18]++ if (substr($align1_rev, $k, 2) =~ /CT/ && substr($align2_rev, $k, 2) =~ /GG/);
					$bi_base_pair[19]++ if (substr($align1_rev, $k, 2) =~ /CC/ && substr($align2_rev, $k, 2) =~ /GG/);
					$bi_base_pair[20]++ if (substr($align1_rev, $k, 2) =~ /CG/ && substr($align2_rev, $k, 2) =~ /GC/);
					$bi_base_pair[21]++ if (substr($align1_rev, $k, 2) =~ /CG/ && substr($align2_rev, $k, 2) =~ /GU|GT/);
					$bi_base_pair[22]++ if (substr($align1_rev, $k, 2) =~ /GA/ && substr($align2_rev, $k, 2) =~ /CG|CG/);
					$bi_base_pair[23]++ if (substr($align1_rev, $k, 2) =~ /GA/ && substr($align2_rev, $k, 2) =~ /UG|TG/);
					$bi_base_pair[24]++ if (substr($align1_rev, $k, 2) =~ /GT/ && substr($align2_rev, $k, 2) =~ /CA/);
					$bi_base_pair[25]++ if (substr($align1_rev, $k, 2) =~ /GT/ && substr($align2_rev, $k, 2) =~ /UA|TA/);
					$bi_base_pair[26]++ if (substr($align1_rev, $k, 2) =~ /GT/ && substr($align2_rev, $k, 2) =~ /CG/);
					$bi_base_pair[27]++ if (substr($align1_rev, $k, 2) =~ /GC/ && substr($align2_rev, $k, 2) =~ /UG|TG/);
					$bi_base_pair[28]++ if (substr($align1_rev, $k, 2) =~ /GC/ && substr($align2_rev, $k, 2) =~ /CG/);
					$bi_base_pair[29]++ if (substr($align1_rev, $k, 2) =~ /GG/ && substr($align2_rev, $k, 2) =~ /CC/);
					$bi_base_pair[30]++ if (substr($align1_rev, $k, 2) =~ /GG/ && substr($align2_rev, $k, 2) =~ /UC|TC/);
					$bi_base_pair[31]++ if (substr($align1_rev, $k, 2) =~ /GG/ && substr($align2_rev, $k, 2) =~ /CU|CT/);
				} elsif (substr($align1_rev, $k, 2) =~ /(-+)/) {
					$bi_seed_end += length($1);
				}
			}
		} else {
			for (my $k = 0; $k < $bi_seed_end; $k++) {
				if (substr($align1_rev, $k, 2) =~ /[A-Za-z]{2}/  && substr($match_letter_rev, $k, 2) =~ /\S\S/) {
					$bi_base_pair[0]++ if (substr($align1_rev, $k, 2) =~ /AA/ && substr($align2_rev, $k, 2) =~ /UU|TT/);
					$bi_base_pair[1]++ if (substr($align1_rev, $k, 2) =~ /AT/ && substr($align2_rev, $k, 2) =~ /UA|TA/);
					$bi_base_pair[2]++ if (substr($align1_rev, $k, 2) =~ /AT/ && substr($align2_rev, $k, 2) =~ /UG|TG/);
					$bi_base_pair[3]++ if (substr($align1_rev, $k, 2) =~ /AC/ && substr($align2_rev, $k, 2) =~ /UG|TG/);
					$bi_base_pair[4]++ if (substr($align1_rev, $k, 2) =~ /AG/ && substr($align2_rev, $k, 2) =~ /UC|TC/);
					$bi_base_pair[5]++ if (substr($align1_rev, $k, 2) =~ /AG/ && substr($align2_rev, $k, 2) =~ /UU|TT/);
					$bi_base_pair[6]++ if (substr($align1_rev, $k, 2) =~ /TA/ && substr($align2_rev, $k, 2) =~ /AU|AT/);
					$bi_base_pair[7]++ if (substr($align1_rev, $k, 2) =~ /TA/ && substr($align2_rev, $k, 2) =~ /GU|GT/);
					$bi_base_pair[8]++ if (substr($align1_rev, $k, 2) =~ /TT/ && substr($align2_rev, $k, 2) =~ /AA/);
					$bi_base_pair[9]++ if (substr($align1_rev, $k, 2) =~ /TT/ && substr($align2_rev, $k, 2) =~ /GA/);
					$bi_base_pair[10]++ if (substr($align1_rev, $k, 2) =~ /TT/ && substr($align2_rev, $k, 2) =~ /AG/);
					$bi_base_pair[11]++ if (substr($align1_rev, $k, 2) =~ /TC/ && substr($align2_rev, $k, 2) =~ /AG/);
					$bi_base_pair[12]++ if (substr($align1_rev, $k, 2) =~ /TC/ && substr($align2_rev, $k, 2) =~ /GG/);
					$bi_base_pair[13]++ if (substr($align1_rev, $k, 2) =~ /TG/ && substr($align2_rev, $k, 2) =~ /AC/);
					$bi_base_pair[14]++ if (substr($align1_rev, $k, 2) =~ /TG/ && substr($align2_rev, $k, 2) =~ /GC/);
					$bi_base_pair[15]++ if (substr($align1_rev, $k, 2) =~ /TG/ && substr($align2_rev, $k, 2) =~ /AU|AT/);
					$bi_base_pair[16]++ if (substr($align1_rev, $k, 2) =~ /CA/ && substr($align2_rev, $k, 2) =~ /GU|GT/);
					$bi_base_pair[17]++ if (substr($align1_rev, $k, 2) =~ /CT/ && substr($align2_rev, $k, 2) =~ /GA/);
					$bi_base_pair[18]++ if (substr($align1_rev, $k, 2) =~ /CT/ && substr($align2_rev, $k, 2) =~ /GG/);
					$bi_base_pair[19]++ if (substr($align1_rev, $k, 2) =~ /CC/ && substr($align2_rev, $k, 2) =~ /GG/);
					$bi_base_pair[20]++ if (substr($align1_rev, $k, 2) =~ /CG/ && substr($align2_rev, $k, 2) =~ /GC/);
					$bi_base_pair[21]++ if (substr($align1_rev, $k, 2) =~ /CG/ && substr($align2_rev, $k, 2) =~ /GU|GT/);
					$bi_base_pair[22]++ if (substr($align1_rev, $k, 2) =~ /GA/ && substr($align2_rev, $k, 2) =~ /CG|CG/);
					$bi_base_pair[23]++ if (substr($align1_rev, $k, 2) =~ /GA/ && substr($align2_rev, $k, 2) =~ /UG|TG/);
					$bi_base_pair[24]++ if (substr($align1_rev, $k, 2) =~ /GT/ && substr($align2_rev, $k, 2) =~ /CA/);
					$bi_base_pair[25]++ if (substr($align1_rev, $k, 2) =~ /GT/ && substr($align2_rev, $k, 2) =~ /UA|TA/);
					$bi_base_pair[26]++ if (substr($align1_rev, $k, 2) =~ /GT/ && substr($align2_rev, $k, 2) =~ /CG/);
					$bi_base_pair[27]++ if (substr($align1_rev, $k, 2) =~ /GC/ && substr($align2_rev, $k, 2) =~ /UG|TG/);
					$bi_base_pair[28]++ if (substr($align1_rev, $k, 2) =~ /GC/ && substr($align2_rev, $k, 2) =~ /CG/);
					$bi_base_pair[29]++ if (substr($align1_rev, $k, 2) =~ /GG/ && substr($align2_rev, $k, 2) =~ /CC/);
					$bi_base_pair[30]++ if (substr($align1_rev, $k, 2) =~ /GG/ && substr($align2_rev, $k, 2) =~ /UC|TC/);
					$bi_base_pair[31]++ if (substr($align1_rev, $k, 2) =~ /GG/ && substr($align2_rev, $k, 2) =~ /CU|CT/);
				} elsif (substr($align1_rev, $k, 2) =~ m/(-+)/) {
					$bi_seed_end += length($1);
				}
			} 
		}
		for (my $i = 0; $i < $#bi_base_pair; $i++) {
			$bi_base_pair_sum += $bi_base_pair[$i];
		}
		# &end('seq');
		#next if ($bi_base_pair_sum == 0);
		
		#my $align_filename = 't/data/$resultFile.seq';
		#generate the species tree automatically using a Bio::DB::Taxonomy database
		#$tdb = Bio::DB::Taxonomy->new(-source => 'entrez');
		#@features = $factory->run("$resultFile.seq", $tdb);
		
		
		
		
		
		#Secondary Structure Accessibility
		my $window_size = 160;
		#for (my $i = 0; $i < 40; $i++) {
			#my $seq_access = substr($UTRseq, ($seq2_start - 20 + $i) - $window_size * 0.5, $window_size);
		my $judge = ($seq2_start - $window_size * 0.5 >= 0) ? ($seq2_start - $window_size * 0.5) : 0;
		my $seq_access = substr($UTRseq, $judge, $window_size);
		open (FHDW4, ">$RESULT_OUTPUT_PATHWAY.tmp") || die "$!\n";
		print FHDW4 $seq_access."\n";
			
#		system "C:\\AppServ\\www\\ViennaRNA-1.8.5_win\\Progs\\RNAplfold.exe -W 80 -L 40 -u 8 < $RESULT_OUTPUT_PATHWAY.tmp";
#		system "/media/data/nfs/home/kang/public_html/miRNAServer/HumanMain/Tool/ViennaRNA-1.8.5_win/Progs/RNAplfold -W 80 -L 40 -u 8 < $RESULT_OUTPUT_PATHWAY.tmp";
#		system "C:\\AppServ\\www\\PlantServer\\Tool\\ViennaRNA-1.8.5_win\\Progs\\RNAplfold.exe -W 80 -L 40 -u 8 < $RESULT_OUTPUT_PATHWAY.tmp";
		# &start;
		system "RNAplfold -W ".$opt{"window-size"}." -L 40 -u 8 < $RESULT_OUTPUT_PATHWAY.tmp";
		# &end("RNAplfold");
#		system "/nfs/master/kang/Plant/Tool/ViennaRNA-1.8.5/Progs/RNAplfold -W 80 -L 40 -u 8 < $RESULT_OUTPUT_PATHWAY.tmp";
#		system "/media/data/nfs/home/kang/PlantMain/Tool/ViennaRNA-1.8.5/Progs/RNAplfold -W 80 -L 40 -u 8 < $RESULT_OUTPUT_PATHWAY.tmp";
		
		open (FHDR4, "plfold_lunp") || die "$!\n";
		my $access_count = 0;
		my @accessibility = (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
		#print "----------------Alignment $align_num:--------------------\n";
		while (<FHDR4>) {
			chomp;
			my $label = 0;
			$label = $1 if (/^(\d+)\s/);
			if ($label >= 60 && $label <= 100) {
				#for (my $i = 1; $i <= 40; $i++) {
				$accessibility[$access_count] = $1 if (/^.*\s(.*)\s.*\s.*\s.*\s.*\s.*\s.*\s.*\s$/);
				#print "$label: $accessibility[$access_count]\n";
				$access_count++;
			}
		}
		#my $access_sum = 0;
		#$access_sum += $_ foreach (@accessibility);
		#print "Averaged accessibility score:".($access_sum/$#accessibility)."\t($access_sum/$#accessibility)\n";
		close (FHDW4);
		close (FHDR4);
		
		
		
		
#------------------------------------------------------------------------------------------------------------------------------------------		
#		#Retrieve Conservation score phastCons
#		my ($chr_location, $conservation_score_seed, $conservation_score_seedout, $seed_leng) = (0, 0, 0, 0);
#		
#		for (my $i = 0; $i < $UTR_leng; $i++) {
#			#if ($eachline[$seq2_start + $i] =~ m/^$chr_location/) {
#			my @eachtab = split /\t/, $eachline[$seq2_start + $i + 1];
#			$chr_location = $chr_start[$pair_number-1] + $seq2_start + $i;
#			#print "chr:$eachtab[0]\tget its scores => $eachtab[1]\n";
#			#print $chr_location."\n";
#			if ($seq1_start =~ 1) {
#				if ($i >= 1 && $i <= 7) {
#					$conservation_score_seed += $eachtab[1];
#					$seed_leng++;
#					#print $eachtab[1]."\tV\n";
#				} elsif ($i > 6) {
#					$conservation_score_seedout += $eachtab[1];
#					#print $eachtab[1]."\n";
#				}
#			} else {
#				if ($i >= 0 && $i <= 6) {
#					$conservation_score_seed += $eachtab[1];
#					$seed_leng++;
#					#print $eachtab[1]."\tV\n";
#				} elsif ($i > 6) {
#					$conservation_score_seedout += $eachtab[1];
#					#print $eachtab[1]."\n";
#				}
#			}
#		}
#		#print "Seed/Seedout Length: $seed_leng\t".($UTR_leng - $seed_leng)."\n";
#		$conservation_score_seed /= $seed_leng;
#		$conservation_score_seedout /= ($UTR_leng - $seed_leng);
#		#print "Conservation Score: $conservation_score_seed\t$conservation_score_seedout\n";
#		#if ($again =~ 0) {
#		#next if ($conservation_score_seed == 0 && $conservation_score_seedout == 0);
#		#}
#---------------------------------------------------------------------------------------------------------------------------------------		
		# &start;
		my $compactness;
		if ($UTR_leng/$miRNA_leng < 1) {
			$compactness = ($seq1_pair_proportion + $seq2_pair_proportion * ( cos((PI * (1 - $UTR_leng/$miRNA_leng) / 2)) ** 0.5)) / 2;
		} else {
			$compactness = ($seq1_pair_proportion + $seq2_pair_proportion * 1) / 2;
		}
		
		#position specific feature
		my $bin_leng = int(length($UTRseq)/10) + 1;
		my $position = int($seq2_start / $bin_leng) + 1;
		#print "length of bin is $bin_leng\n";
		#print "alignment is located at #$position bin\n";
		
		#bulge-related features
		my ($bulge_leng_miRNA_max, $bulge_leng_miRNA_tmp, $bulge_leng_UTR_max, $bulge_leng_UTR_tmp, $bulge_num, $bulge_miRNA_num, $bulge_UTR_num) = (0, 0, 0, 0, 0, 0, 0);
		for (my $i = 0; $i < length($align1); $i++) {
			if (substr($align1, $i, 1) =~ /-/ && substr($align1, $i - 1, 1) =~ /[^\-]/) {
				$bulge_leng_miRNA_tmp = 1;
				$bulge_miRNA_num++;
				$bulge_num++;
				for (my $j = 1; $j < length($align1) - $i; $j++) {
					$bulge_leng_miRNA_tmp++ if (substr($align1, $i + $j, 1) =~ /-/);
					last if (substr($align1, $i + $j, 1) =~ /[^\-]/);
				}
				$bulge_leng_miRNA_max = $bulge_leng_miRNA_tmp if ($bulge_leng_miRNA_tmp > $bulge_leng_miRNA_max);
			}
		} 
		for (my $i = 0; $i < length($align2); $i++) {
			if (substr($align2, $i, 1) =~ /-/ && substr($align2, $i - 1, 1) =~ /[^\-]/) {
				$bulge_num++;
				$bulge_UTR_num++;
				$bulge_leng_UTR_tmp = 1;
				for (my $j = 1; $j < length($align2) - $i; $j++) {
					$bulge_leng_UTR_tmp++ if (substr($align2, $i + $j, 1) =~ /-/);
					last if (substr($align2, $i + $j, 1) =~ /[^\-]/);
				}
				$bulge_leng_UTR_max = $bulge_leng_UTR_tmp if ($bulge_leng_UTR_tmp > $bulge_leng_UTR_max);
			}
		} 
		#print "max bulge length miRNA=$bulge_leng_miRNA_max\tmax bulge length UTR=$bulge_leng_UTR_max\tnumber of bulge=$bulge_num\n";
		
		#if ($again =~ 0) {
			#print FHDW "Seq1 start-end at $seq1_start - $seq1_end;\t\tPairs Number/miRNA Alignment Length: $pair_num/$miRNA_leng = $seq1_pair_proportion\n";			
			print FHDW "--------------------------------------------------------------------------------------------------------------------------------\n";
			print FHDW "Seq1 start-end at $seq1_start - $seq1_end;\t\tPairs Number/miRNA Alignment Length: $pair_num/$miRNA_leng = $seq1_pair_proportion\n";
			print FHDW "Seq2 start-end at $seq2_start - $seq2_end;\tPairs Number/3'UTR Alignment Length: $pair_num/$UTR_leng = $seq2_pair_proportion\n";
			
			print FHDW "Seq1: 3' $align1 5'\n";
			print FHDW "         $match_letter\n";
			print FHDW "Seq2: 5' $align2 3'\nScore: $max_score\n";
			
			
			
			#print FHDW "Seq1:\t3'\t$align1\t5'\n";
			#print FHDW "\t\t\t$match_letter\n";
			#print FHDW "Seq2:\t5'\t$align2\t3'\nScore: $max_score\n";
			print FHDW "Mismatches: $mismatch_num;\tGaps: $gap_num;\tWobble pairs: $wobble_num\n";
			print FHDW "G+C content of target site: $GC_content; A+T content of target site: $AU_content\n";
			#print FHDW "Distance: $distance/".length($UTRseq)." = $dis_percent\n";
			print FHDW "Distance to the nearer end of 3'UTR: $nearer_end\tRelative Position in UTR: ".$seq2_start/length($UTRseq)."\t($seq2_start/".length($UTRseq).")\n";
			print FHDW "3' Binding/3' Supplementary Site: $three_binding/5\n";
			print FHDW "Structure: $structure\t\tEnergy: $energy\n";
			
			print FHDW "Local Flanking AU Composition: $flanking_AU_six(6_mer)\t$flanking_AU_m8(7mer_m8)\t$flanking_AU_1A(7mer_1A)\t$flanking_AU_eight(8mer)\n";
			print FHDW "Frequency of Single Nuclei at Seed Matching Site:\tA: $single_nuclei[0]\tT: $single_nuclei[1]\tC: $single_nuclei[2]\tG: $single_nuclei[3]\n";
			print FHDW "Frequency of Di-Nucleis at Seed Matching Site:\tAA:$di_nuclei[0]\tAT:$di_nuclei[1]\tAC:$di_nuclei[2]\tAG:$di_nuclei[3]\tTA:$di_nuclei[4]\tTT:$di_nuclei[5]\tTC:$di_nuclei[6]\tTG:$di_nuclei[7]";
			print FHDW "\tCA:$di_nuclei[8]\tCT:$di_nuclei[9]\tCC:$di_nuclei[10]\tCG:$di_nuclei[11]\tGA:$di_nuclei[12]\tGT:$di_nuclei[13]\tGC:$di_nuclei[14]\tGG:$di_nuclei[15]\n";
			
			#print FHDW "Frequency of Single Nuclei in Seed-out Matching:\tA: $single_nuclei_out[0]\tT: $single_nuclei_out[1]\tC: $single_nuclei_out[2]\tG: $single_nuclei_out[3]\n";
			#print FHDW "Frequency of Di-Nucleis in Seed-out Matching:\tAA:$di_nuclei_out[0]\tAT:$di_nuclei_out[1]\tAC:$di_nuclei_out[2]\tAG:$di_nuclei_out[3]\tTA:$di_nuclei_out[4]\tTT:$di_nuclei_out[5]\tTC:$di_nuclei_out[6]\tTG:$di_nuclei_out[7]";
			#print FHDW "\tCA:$di_nuclei_out[8]\tCT:$di_nuclei_out[9]\tCC:$di_nuclei_out[10]\tCG:$di_nuclei_out[11]\tGA:$di_nuclei_out[12]\tGT:$di_nuclei_out[13]\tGC:$di_nuclei_out[14]\tGG:$di_nuclei_out[15]\n";
			
			print FHDW "Frequency of Base Pairing at Seed Matching Site:\tA-t: $pair_AU\tU-a: $pair_UA\tC-g: $pair_CG\tG-c: $pair_GC\tG-t: $pair_GU\tU-g: $pair_UG\n";
			print FHDW "6mer match:$six_mer\t7mer_m8 match:$seven_mer_m8\t7mer_1A match:$seven_mer_1A\t8mer match:$eight_mer\tNoncanonical:$noncanonical\n";
#			print FHDW "Conservation Score of Seed: $conservation_score_seed\tConservation Score of Seed-out: $conservation_score_seedout\n";
			print FHDW "Watson-Crick pairing: $WC_pair_six(6mer)\t$WC_pair_m8(7mer_m8)\t$WC_pair_1A(7mer_1A)\t$WC_pair_eight(8mer)\n";
			print FHDW "alignment is located at #$position bin\n";
			print FHDW "max bulge length miRNA=$bulge_leng_miRNA_max\tmax bulge length UTR=$bulge_leng_UTR_max\tnumber of bulge=$bulge_num\tnumber of bulge on seed site=$bulge_miRNA_num\tnumber of bulge on target site=$bulge_UTR_num\n\n";
			
		#}
		#create feature
		$feature[0] = $six_mer;
		$feature[1] = $flanking_AU_six;
		$feature[2] = $WC_pair_six;
		$feature[3] = $seven_mer_m8;
		$feature[4] = $flanking_AU_m8;
		$feature[5] = $WC_pair_m8;
		$feature[6] = $seven_mer_1A;
		$feature[7] = $flanking_AU_1A;
		$feature[8] = $WC_pair_1A;
		$feature[9] = $eight_mer;
		$feature[10] = $flanking_AU_eight;
		$feature[11] = $WC_pair_eight;
		
		$feature[12] = $single_A;
		$feature[13] = $single_T;
		$feature[14] = $single_C;
		$feature[15] = $single_G;

		$feature[16] = $single_A_out;
		$feature[17] = $single_T_out;
		$feature[18] = $single_C_out;
		$feature[19] = $single_G_out;
		
		$feature[20] = $di_AA;
		$feature[21] = $di_AT;
		$feature[22] = $di_AC;
		$feature[23] = $di_AG;
		$feature[24] = $di_TA;
		$feature[25] = $di_TT;
		$feature[26] = $di_TC;
		$feature[27] = $di_TG;
		$feature[28] = $di_CA;
		$feature[29] = $di_CT;
		$feature[30] = $di_CC;
		$feature[31] = $di_CG;
		$feature[32] = $di_GA;
		$feature[33] = $di_GT;
		$feature[34] = $di_GC;
		$feature[35] = $di_GG;
		
		$feature[36] = $di_AA_out;
		$feature[37] = $di_AT_out;
		$feature[38] = $di_AC_out;
		$feature[39] = $di_AG_out;
		$feature[40] = $di_TA_out;
		$feature[41] = $di_TT_out;
		$feature[42] = $di_TC_out;
		$feature[43] = $di_TG_out;
		$feature[44] = $di_CA_out;
		$feature[45] = $di_CT_out;
		$feature[46] = $di_CC_out;
		$feature[47] = $di_CG_out;
		$feature[48] = $di_GA_out;
		$feature[49] = $di_GT_out;
		$feature[50] = $di_GC_out;
		$feature[51] = $di_GG_out;

		
		$feature[52] = $pair_AU;
		$feature[53] = $pair_UA;
		$feature[54] = $pair_CG;
		$feature[55] = $pair_GC;
		$feature[56] = $pair_GU;
		$feature[57] = $pair_UG;
		$feature[58] = $bi_base_pair[0]/$di_nuclei_sum;
		$feature[59] = $bi_base_pair[1]/$di_nuclei_sum;
		$feature[60] = $bi_base_pair[2]/$di_nuclei_sum;
		$feature[61] = $bi_base_pair[3]/$di_nuclei_sum;
		$feature[62] = $bi_base_pair[4]/$di_nuclei_sum;
		$feature[63] = $bi_base_pair[5]/$di_nuclei_sum;
		$feature[64] = $bi_base_pair[6]/$di_nuclei_sum;
		$feature[65] = $bi_base_pair[7]/$di_nuclei_sum;
		$feature[66] = $bi_base_pair[8]/$di_nuclei_sum;
		$feature[67] = $bi_base_pair[9]/$di_nuclei_sum;
		$feature[68] = $bi_base_pair[10]/$di_nuclei_sum;
		$feature[69] = $bi_base_pair[11]/$di_nuclei_sum;
		$feature[70] = $bi_base_pair[12]/$di_nuclei_sum;
		$feature[71] = $bi_base_pair[13]/$di_nuclei_sum;
		$feature[72] = $bi_base_pair[14]/$di_nuclei_sum;
		$feature[73] = $bi_base_pair[15]/$di_nuclei_sum;
		$feature[74] = $bi_base_pair[16]/$di_nuclei_sum;
		$feature[75] = $bi_base_pair[17]/$di_nuclei_sum;
		$feature[76] = $bi_base_pair[18]/$di_nuclei_sum;
		$feature[77] = $bi_base_pair[19]/$di_nuclei_sum;
		$feature[78] = $bi_base_pair[20]/$di_nuclei_sum;
		$feature[79] = $bi_base_pair[21]/$di_nuclei_sum;
		$feature[80] = $bi_base_pair[22]/$di_nuclei_sum;
		$feature[81] = $bi_base_pair[23]/$di_nuclei_sum;
		$feature[82] = $bi_base_pair[24]/$di_nuclei_sum;
		$feature[83] = $bi_base_pair[25]/$di_nuclei_sum;
		$feature[84] = $bi_base_pair[26]/$di_nuclei_sum;
		$feature[85] = $bi_base_pair[27]/$di_nuclei_sum;
		$feature[86] = $bi_base_pair[28]/$di_nuclei_sum;
		$feature[87] = $bi_base_pair[29]/$di_nuclei_sum;
		$feature[88] = $bi_base_pair[30]/$di_nuclei_sum;
		$feature[89] = $bi_base_pair[31]/$di_nuclei_sum;

		$feature[90] = $compactness;
		$feature[91] = $mismatch_num/$pair_num;
		$feature[92] = $gap_num/$pair_num;
		$feature[93] = $wobble_num/$pair_num;
		$feature[94] = $GC_content/$pair_num;
		$feature[95] = $AU_content/$pair_num;
		$feature[96] = $max_score;
		$feature[97] = $Up_flanking_AU;
		$feature[98] = $Down_flanking_AU;
		$feature[99] = $nearer_end;
		$feature[100] = $position;
		$feature[101] = $three_binding/5;
		$feature[102] = $AU_composition_score;
		$feature[103] = $Up_AU_score;
		$feature[104] = $Down_AU_score;
		$feature[105] = ($accessibility[0] + $accessibility[1]) / 2;
		$feature[106] = ($accessibility[2] + $accessibility[3]) / 2;
		$feature[107] = ($accessibility[4] + $accessibility[5]) / 2;
		$feature[108] = ($accessibility[6] + $accessibility[7]) / 2;
		$feature[109] = ($accessibility[8] + $accessibility[9]) / 2;
		$feature[110] = ($accessibility[10] + $accessibility[11]) / 2;
		$feature[111] = ($accessibility[12] + $accessibility[13]) / 2;
		$feature[112] = ($accessibility[14] + $accessibility[15]) / 2;
		$feature[113] = ($accessibility[16] + $accessibility[17]) / 2;
		$feature[114] = ($accessibility[18] + $accessibility[19]) / 2;
		$feature[115] = ($accessibility[20] + $accessibility[21]) / 2;
		$feature[116] = ($accessibility[22] + $accessibility[23]) / 2;
		$feature[117] = ($accessibility[24] + $accessibility[25]) / 2;
		$feature[118] = ($accessibility[26] + $accessibility[27]) / 2;
		$feature[119] = ($accessibility[28] + $accessibility[29]) / 2;
		$feature[120] = ($accessibility[30] + $accessibility[31]) / 2;
		$feature[121] = ($accessibility[32] + $accessibility[33]) / 2;
		$feature[122] = ($accessibility[34] + $accessibility[35]) / 2;
		$feature[123] = ($accessibility[36] + $accessibility[37]) / 2;
		$feature[124] = ($accessibility[38] + $accessibility[39]) / 2;
#		$feature[125] = 0;
#		$feature[126] = 0;
#		$feature[125] = $conservation_score_seed;
#		$feature[126] = $conservation_score_seedout;
		$feature[125] = length($UTRseq); # 127
		$feature[126] = $bulge_leng_miRNA_max; # 128
		$feature[127] = $bulge_leng_UTR_max; # 129
		$feature[128] = $bulge_num; # 130
		$feature[129] = $bulge_miRNA_num; # 131
		$feature[130] = $bulge_UTR_num; # 132
		$feature[131] = $distance; # 133
		$feature[132] = $noncanonical; # 134
		$feature[133] = $flanking_AU_noncanonical; # 135
		$feature[134] = $WC_pair_noncanonical; # 136

		
		for($i=135; $i< 199; $i++){ # 137 201
			$feature[$i] = $tri_nuclei_frequency[$i-135];
			#$feature[$i] = 0;
			#print FHDW "$feature[$i]";
		}
		
		for($i=199; $i< 263; $i++){ # 201 265
			$feature[$i] = $tri_nuclei_out_frequency[$i-201];
			#$$feature[$i] = 0;
		}
		#$feature[265] = 0;
		$feature[264] = $central_binding/3;
		$feature[265] = $energy;		
		#$feature[266] = 0;		
=head
		for($i=0; $i< 12; $i++){
			#$feature[$i] = $tri_nuclei_out_frequency[$i-201];
			$feature[$i] = 0;
		}
		for($i=58; $i< 137; $i++){
			#$feature[$i] = $tri_nuclei_out_frequency[$i-201];
			$feature[$i] = 0;
		}		
=cut	
		#}
		#$align_num++;
		#push (@align_starts, $seq2_start);
	#}
	#print FHDW3 "$prediction";
	#for (my $f = 1; $f <= 131; $f++) {
	#	$feature[$f] /= $align_num if ($f >= 13 && $f <= 90 && $align_num != 0);
	#	$feature[$f] = 0 if ($align_num =~ 0);
	#	print FHDW3 "\t$f:$feature[$f]";
	#}
	#print FHDW3 "\n";
	close (FHDW2);	
	#$align_num, @align_starts;
	#print $seq2_start."\n";
	# &end("aa");
	return ($energy, $seq2_start, @feature);
}

sub kmer {
	my ($seed, $seq, $k, $allow) = @_;
	$seed = reverse $seed;
	$seed =~ s/U/T/g;
	my (@hit_all, @hit_compensatory);
	if ($allow) {
		my @t1 = split '', $seed;
		my @t2 = split '', $seq;
		my %err_count;
		for my $i(0..@t2-length($seed)) { # mRNA_length - seed_length times
			my $ne = 0;
			for my $j(0..@t1-1) { # seed_length times
				if ($t1[$j] ne $t2[$i + $j]) {
					$ne++; # mismatch count
				}
				last if ($ne > 1); # if mismatch count > 2 => drop
			}
			if($ne<=1){
				$err_count{$i} = $ne;
			}
		}
		# die Dump %err_count;
		foreach (sort keys %err_count) {
			my $value; # save the perferct match in mRNA
			if ($err_count{$_} == 0) {		#perfect match
				$value = substr($seq, ${_}, length($seed)); # value to save the match mRNA sequence
				push (@hit_all, ${_});
				if(scalar @hit_compensatory>0){
					for (my $j = 0; $j <= $#hit_compensatory; ++$j) { # <hit_compensatory length
						if (abs(${_} - $hit_compensatory[$j]) < 20) { # position within 20
							for my $k($j..$#hit_compensatory-1) {
								$hit_compensatory[$k] = $hit_compensatory[$k+1];
							}
							pop @hit_compensatory;
							--$j;
						}
					}
				}
			}
			if ($err_count{$_} == 1) {		#one mismatch
				$value = substr($seq, ${_}, length($seed));
				my $loc = ${_};
				my $there = 0;
				for (my $j = 0; $j <= $#hit_all; $j++) {
					if (abs($loc - $hit_all[$j]) < length($seed)) {
						$there = 1;
						last;
					}
				}
				for (my $j = 0; $j < $#hit_compensatory; $j++) {
					if (abs($loc - $hit_compensatory[$j]) < length($seed)) {
						$there = 1;
						last;
					}
				}
				if ($there == 0){
					push (@hit_compensatory, $loc);
				} 
			}
		}
	}
	(@hit_all, @hit_compensatory);
}

sub nearest {
	my ($subject, @neighbours)= @_;
	my $distance = "inf";
	my $nearest = -1;
	foreach (@neighbours) {
		if(!$_){
			next;
		}
		if (abs ($subject - $_) < $distance && abs ($subject - $_) != 0) {
			$distance = abs ($subject - $_);
			$nearest = $_;
		}
	}
	$distance, $nearest;
}
