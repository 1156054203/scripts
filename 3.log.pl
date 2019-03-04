open IN,"lncRNA_diff_expr.txt";
open OUT,">log_lncRNA_diff_expr.txt";
while(<IN>){
	chomp;
	if($. == 1){
		print OUT "probeid\t$_\n";
		next;
	}
	@line = split "\t";
	print OUT "$line[0]";
	foreach $count (1..@line-1){
		$temp = log2($line[$count]);
		print OUT "\t$temp";
	}
	print OUT "\n";
}
close IN;
close OUT;

sub log2 {
my $n = shift;
return log($n)/log(2);
}