#!/usr/bin/perl

chomp $ARGV[0];
$file=$ARGV[0];

chomp $ARGV[1];
$pattern=$ARGV[1];

$file=~/.+\/([^\.]+)/;
$pre=$1;
# $new=$pre.".fas";

open (IN, $file) || die "can not open $file\n";
while ($line=<IN>){
	chomp $line;
	if ($line=~/^>(.+)/){
		$acc=$1 if (!$exist{$1});
		$acc="" if ($exist{$1});
		#print "$acc\n";
		$exist{$acc}=1;
	}
	else {
		$seq{$acc}.=$line;
	}
}
close IN;

foreach $acc (keys %seq){
	$rev=reverse $seq{$acc};
	$count=0;
	if (($pattern) && ($acc=~/$pattern/)){
		while ($rev){
			$count++;
			$char=chop $rev;
			$char.=chop $rev;
			$char.=chop $rev;
			$pos[$count].="$char ";
		}
	}
	elsif (!$pattern){
		while ($rev){
			$count++;
			$char=chop $rev;
			$char.=chop $rev;
			$char.=chop $rev;
			$pos[$count].="$char ";
		}
	}

}
$count=0;
$pos=0;
$new=$file.".consensus.unaligned";
$new_aln=$file.".consensus";
print "new files are $new and $new_aln\n";
open (OUT, ">$new") || die "can not open $new\n";
open (ALN, ">$new_aln") || die "can not open $new_aln\n";

print OUT ">$pre Consensus Sequence\n";
print ALN ">$pre Consensus Sequence\n";
foreach $i (@pos){
	$pos++;
	@chars=split(/ /, $i);
	%count=();
	foreach $char (@chars){
		$count{$char}++;
	}
	#print "$i %count"; <STDIN>;
	@sorted=sort by_value keys (%count);
	#if ($pos==1518){
	#	foreach $s (@sorted){
	#		print "$s $count{$s}\n";
	#	}
	#}
	print OUT "$sorted[0]" unless ($sorted[0] eq "---");
	print ALN "$sorted[0]"; #unless ($sorted[0] eq "-");
}
print OUT "\n";
print ALN "\n";

close OUT;
close ALN;

sub by_value {
   ($count{$b} <=> $count{$a}) || ($b cmp $a);
}

