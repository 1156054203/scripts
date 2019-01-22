#! /usr/bin/perl

my ($file,$opt_t)=@ARGV;

my $scriptdir="/online/home/cheny/Apollo_v2/bigdata";
my $usage= <<END;

perl $0 CHGID_count.txt hg19/hg38

END
die  $usage unless @ARGV==2;

if($file=~/(CH\w\d+)_count/){
	$id=$1;
}

my $dir="/offline/data-sharing/vcf/$id";
my $cloud="/home/data/Reports/$id";

# circos in $var/healthcare
if($opt_t eq "hg19"){
system("python $scriptdir/circlize.py  -n $id -s $dir/anno/$id\.snp.tsv -e $dir/anno/$id\.indel.tsv -c $dir/anno/$id\.cnv.tsv -v $dir/anno/$id\.sv.tsv -o $dir/healthcare/IndividualPics") == 0
or die "generate hg19 circos figure failed!";}
if($opt_t eq "hg38"){
system("python $scriptdir/circlize.py -3 -n $id -s $dir/anno/$id\.snp.tsv -e $dir/anno/$id\.indel.tsv -c $dir/anno/$id\.cnv.tsv -v $dir/anno/$id\.sv.tsv -o $dir/healthcare/IndividualPics") == 0
or die "generate hg38 circos figure failed!";
}

# big data pies
system("/online/software/R-3.4.0/lib64/R/bin/Rscript $scriptdir/bigdata_Pie.r $file $dir/healthcare/IndividualPics $id");
#system("Rscript $scriptdir/bigdata_Pie.r $file $dir/healthcare/IndividualPics $id");

open IN,$file;
while(<IN>){
	chomp;
	my @array=split/\t/,$_;
	if($array[0] eq "SNVs"){
		$snv{"count"}=format_counts($array[1]);
		$snv{"gene"}=format_counts($array[2]);
		$snv{"inter"}=format_counts($array[3]);
		$snv{"het"}=format_counts($array[4]);
		$snv{"hom"}=format_counts($array[5]);
	}elsif($array[0] eq "InDels"){
		$indel{"count"}=format_counts($array[1]);
		$indel{"gene"}=format_counts($array[2]);
		$indel{"inter"}=format_counts($array[3]);
		$indel{"het"}=format_counts($array[4]);
		$indel{"hom"}=format_counts($array[5]);
	}elsif($array[0] eq "CNVs"){
		$cnv{"count"}=format_counts($array[1]);
		$cnv{"gene"}=format_counts($array[2]);
		$cnv{"inter"}=format_counts($array[3]);
		$cnv{"dup"}=format_counts($array[4]);
		$cnv{"del"}=format_counts($array[5]);
	}elsif($array[0] eq "SVs"){
		$sv{"count"}=format_counts($array[1]);
		$sv{"gene"}=format_counts($array[2]);
		$sv{"inter"}=format_counts($array[3]);
		$sv{"dup"}=format_counts($array[4]);
		$sv{"del"}=format_counts($array[5]);
		$sv{"inv"}=format_counts($array[6]);
		$sv{"tran"}=format_counts($array[7]);
	}
}
close IN;

my $bigdata= <<END;
%-----------------bigdata
\\CVSection{您的基因全图谱}
尊贵的VIP客户, 您所提交到云健康的生物样本,经全球超高通量的人类全基因组检测平台测序,以及高质量的数据处理和分析,最后由国际化的生物学专家团队进行健康风险评估得到最终的报告。(此报告仅列出与各系统疾病相关的基因变异)\\\\\\\\
云健康已为您检测了人类全部约\\orangetext{25,000}个编码基因、近\\orangetext{100,000}个非编码基因以及基因间区，包括所有的大约\\orangetext{32亿}个基因位点，并为您绘制了全世界独一无二的全基因组“\\orangetext{基因身份证}”。\\\\
\\begin{center}
\\bluetext{您的“基因身份证”\\\\}
\\includegraphics[width=0.7\\textwidth]{$cloud/bundle/circos-$id\.eps}
\\end{center}
\\LargeSep
\\begin{elaboration}
\\bluetext{从最外层开始，每层各代表:}\\\\\\\\
\\ding{192} 染色体整长上的位置标尺。\\\\
\\ding{193} 每个染色体上的条带信息。\\\\
\\ding{194} 全基因组范围上非同义突变的展示，\\redtext{Tv}: 颠换；\\rbluetext{Ti}: 转换。\\\\
\\ding{195} 全基因组范围上的片段\\redtext{增加。}\\\\
\\ding{196} 全基因组范围上的片段\\greentext{缺失。}\\\\
\\ding{197} 全基因组范围上的跨染色体的结构变异。
\\end{elaboration}
\\newpage

%---Gene Big Data

\\CVSection{您的全基因组大数据}
{\\scriptsize \\uline{注：变异数目与疾病风险没有显著相关性。}}\\\\
\\begin{minipage}{0.6\\linewidth}
  \\centerline{\\includegraphics[width=\\textwidth]{$cloud/bundle/$id\_counts.eps}}
\\end{minipage}
\\hfill
\\begin{minipage}{0.5\\linewidth}
  \\scriptsize
  \\pbox{6cm}{
   \\bluetext{总变异数统计}\\\\\\\\
  {\\Large \\pictexta{$cnv{"count"}}}  拷贝数变异(CNVs) \\\\
  {\\Large \\pictextb{$indel{"count"}}}  插入缺失变异(InDels)\\\\
  {\\Large \\pictextc{$snv{"count"}}}  单位点变异(SNVs)\\\\
  {\\Large \\pictextd{$sv{"count"}}}  结构变异(SVs)\\\\
  }
\\end{minipage}


\\begin{minipage}{0.5\\linewidth}
  \\scriptsize
  \\pbox{6cm}{
  \\bluetext{基因区变异统计}\\\\\\\\
  {\\Large \\pictexta{$cnv{"gene"}}}  拷贝数变异(CNVs) \\\\
  {\\Large \\pictextb{$indel{"gene"}}}  插入缺失变异(InDels)\\\\
  {\\Large \\pictextc{$snv{"gene"}}}  单位点变异(SNVs)\\\\
  {\\Large \\pictextd{$sv{"gene"}}}  结构变异(SVs)\\\\
  }
\\end{minipage}
\\hfill
\\begin{minipage}{0.6\\linewidth}
  \\centerline{\\includegraphics[width=\\textwidth]{$cloud/bundle/$id\_genic.eps}}
\\end{minipage}



\\begin{minipage}{0.6\\linewidth}
  \\centerline{\\includegraphics[width=\\textwidth]{$cloud/bundle/$id\_intergenic.eps}}
\\end{minipage}
\\hfill
\\begin{minipage}{0.5\\linewidth}
  {\\scriptsize
  \\pbox{6cm}{
  \\bluetext{基因间区变异统计}\\\\\\\\
  {\\Large \\pictexta{$cnv{"inter"}}}  拷贝数变异(CNVs) \\\\
  {\\Large \\pictextb{$indel{"inter"}}}  插入缺失变异(InDels)\\\\
  {\\Large \\pictextc{$snv{"inter"}}}  单位点变异(SNVs)\\\\
  {\\Large \\pictextd{$sv{"inter"}}}  结构变异(SVs)\\\\
  }}
\\end{minipage}


%----- statistical vars
\\newpage

\\begin{center}  %  单位点变异
\\bluetext{\\LARGE $snv{"count"}} 个单位点变异~(SNVs)\\\\
\\end{center}
\\begin{minipage}{0.55\\linewidth}
  \\centerline{\\includegraphics[width=\\textwidth]{$cloud/bundle/$id\_snvg.eps}}
  \\centerline{\\pbox{6cm}{\\scriptsize
 \\bluebullet  {\\normalsize \\pictexta{$snv{"gene"}}} 基因区(Genic)\\\\
 \\bluebullet {\\normalsize \\pictextc{$snv{"inter"}}} 基因间区(Intergenic)\\\\
  }}
\\end{minipage}
\\hfill
\\begin{minipage}{0.55\\linewidth}
  \\centerline{\\includegraphics[width=\\textwidth]{$cloud/bundle/$id\_snvh.eps}}
  \\centerline{\\pbox{6cm}{\\scriptsize
  \\bluebullet {\\normalsize \\pictexta{$snv{"het"}}} 杂合变异(Het)\\\\
 \\bluebullet  {\\normalsize \\pictextc{$snv{"hom"}}} 纯和变异(Hom)\\\\
  }}
\\end{minipage}
\\LargeSep
\\begin{center}  % 插入缺失变异
\\bluetext{\\LARGE $indel{"count"}} 个插入缺失变异~(Indels)\\\\
\\end{center}
\\begin{minipage}{0.55\\linewidth}
  \\centerline{\\includegraphics[width=\\textwidth]{$cloud/bundle/$id\_indg.eps}}
  \\centerline{\\pbox{6cm}{\\scriptsize
  \\bluebullet {\\normalsize \\pictexta{$indel{"gene"}}} 基因区(Genic)\\\\
  \\bluebullet {\\normalsize \\pictextc{$indel{"inter"}}} 基因间区(Intergenic)\\\\
  }}
\\end{minipage}
\\hfill
\\begin{minipage}{0.55\\linewidth}
  \\centerline{\\includegraphics[width=\\textwidth]{$cloud/bundle/$id\_indh.eps}}
  \\centerline{\\pbox{6cm}{\\scriptsize
  \\bluebullet {\\normalsize \\pictexta{$indel{"het"}}} 杂合变异(Het)\\\\
  \\bluebullet {\\normalsize \\pictextc{$indel{"hom"}}} 纯和变异(Hom)\\\\
  }}
\\end{minipage}
\\MidSep
\\MidSep
% --- 知识小天地
\\begin{elaboration}
\\scriptsize
\\CVItem{单位点变异(SNVs) }{个体间基因组DNA序列同一位置单个核苷酸变化(替代、插入或缺失)所导致的变异。不同物种、个体基因组DNA序列同一位置上的单个核苷酸都可能存在差别。}\\\\
\\CVItem{插入缺失变异(InDels)}{插入和缺失都是相对于人类参考序列来说的。缺失，指在参考序列中的某一段DNA丢失；插入，指相对于参考序列在某个位置多了一段DNA序列，即在这个位置插入了一段序列。}
\\end{elaboration}
\\newpage


\\begin{center}  % CNV变异
\\bluetext{\\LARGE $cnv{"count"}} 个拷贝数变异~(CNVs)\\\\
\\end{center}
\\begin{minipage}{0.55\\linewidth}
  \\centerline{\\includegraphics[width=\\textwidth]{$cloud/bundle/$id\_cnvg.eps}}
  \\centerline{\\pbox{6cm}{\\scriptsize
  \\bluebullet {\\normalsize \\pictexta{$cnv{"gene"}}}基因区(Genic)\\\\
  \\bluebullet {\\normalsize \\pictextc{$cnv{"inter"}}}基因间区(Intergenic)\\\\
  }}
\\end{minipage}
\\hfill
\\begin{minipage}{0.55\\linewidth}
  \\centerline{\\includegraphics[width=\\textwidth]{$cloud/bundle/$id\_cnvi.eps}}
  \\centerline{\\pbox{6cm}{\\scriptsize
  \\bluebullet {\\normalsize \\pictexta{$cnv{"del"}}}缺失变异(Deletion)\\\\
  \\bluebullet {\\normalsize \\pictextc{$cnv{"dup"}}}插入变异(Insertion)\\\\
  }}
\\end{minipage}
\\LargeSep
\\begin{center}  % SV变异
\\bluetext{\\LARGE $sv{"count"}} 个结构变异~(SVs)\\\\
\\end{center}
\\begin{minipage}{0.55\\linewidth}
  \\centerline{\\includegraphics[width=\\textwidth]{$cloud/bundle/$id\_svg.eps}}
  \\centerline{\\pbox{6cm}{\\scriptsize
  \\bluebullet {\\normalsize \\pictexta{$sv{"gene"}}}基因区(Genic)\\\\
  \\bluebullet {\\normalsize \\pictextc{$sv{"inter"}}}基因间区(Intergenic)\\\\
  }}
\\end{minipage}
\\hfill
\\begin{minipage}{0.55\\linewidth}
  \\centerline{\\includegraphics[width=\\textwidth]{$cloud/bundle/$id\_svi.eps}}
  \\centerline{\\pbox{6cm}{\\scriptsize
  \\bluebullet {\\normalsize \\pictexta{$sv{"del"}}}缺失变异(Deletion)\\\\
  \\bluebullet {\\normalsize \\pictextb{$sv{"dup"}}}插入变异(Insertion)\\\\
  \\bluebullet {\\normalsize \\pictextc{$sv{"inv"}}}倒位变异(Invertion)\\\\
  \\bluebullet {\\normalsize \\pictextd{$sv{"tran"}}}易位变异(Translocation)\\\\
  }}
\\end{minipage}
\\MidSep
%---- 知识小天地
\\begin{elaboration}
\\scriptsize
\\CVItem{拷贝数变异(CNVs)}{拷贝数目变异也称拷贝数目多态，是一种大小介于1kb至3Mb的DNA片段的变异。例如人类正常染色体拷贝数是2，有些染色体区域拷贝数变成1或3，这样，该区域就发生拷贝数缺失或增加，位于该区域内的基因表达量也会受到影响。}\\\\
\\CVItem{结构变异(SVs)}{染色体结构变异是指在染色体上发生了大片段的变异。主要包括染色体大片段的插入和缺失（CNV），染色体内部的某块区域发生翻转颠换，两条染色体之间发生重组等。}
\\end{elaboration}
\\newpage

END


$out=$dir."/healthcare/$id\_bigdata.txt";

open OUT,">$out";
print OUT $bigdata;
close OUT;



sub format_counts{
	my $tmp = $_[0];
	$tmp =~ s/(\d+?)(?=(?:\d{3})+$)/$1,/g;
	return $tmp;
}
