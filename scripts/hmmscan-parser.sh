#!/usr/bin/env sh

# Yanbin Yin
# 07/21/2015
# hmmscan output parser
# Usage: sh hmmscan-parser.sh hmmscan-output-file

# 1. take hmmer3 --domtblout output and extract necessary columns
# 2. sort on the protein position columns
# 3. remove overlapped/redundant hmm matches and keep the one with the lower e-values
# 4. calculate the covered fraction of hmm &
#    apply the E-value cutoff and the covered faction cutoff

cat $1 | grep -v '^#' | awk '{print $1,$3,$4,$6,$13,$16,$17,$18,$19}' | sed 's/ /\t/g' | \
	sort -k 3,3 -k 8n -k 9n | \
	perl -e 'while(<>){chomp;@a=split;next if $a[-1]==$a[-2];push(@{$b{$a[2]}},$_);}foreach(sort keys %b){@a=@{$b{$_}};for($i=0;$i<$#a;$i++){@b=split(/\t/,$a[$i]);@c=split(/\t/,$a[$i+1]);$len1=$b[-1]-$b[-2];$len2=$c[-1]-$c[-2];$len3=$b[-1]-$c[-2];if($len3>0 and ($len3/$len1>0.5 or $len3/$len2>0.5)){if($b[4]<$c[4]){splice(@a,$i+1,1);}else{splice(@a,$i,1);}$i=$i-1;}}foreach(@a){print $_."\n";}}' | \
	perl -e 'while(<>){chomp;@a=split(/\t/,$_);if(($a[-1]-$a[-2])>80){print $_,"\t",($a[-3]-$a[-4])/$a[1],"\n" if $a[4]<1e-5;}else{print $_,"\t",($a[-3]-$a[-4])/$a[1],"\n" if $a[4]<1e-3;}}' | awk '$NF>0.3' | sort -k 3 -k 8,9g
