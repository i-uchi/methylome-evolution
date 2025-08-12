#!/usr/bin/perl -s

use FileHandle;
use File::Basename;

$seq_file = $ARGV[0];	## converted sequence file
$align_dir = $ARGV[1];	## directory containig reference alignments
$outdir = $ARGV[2];	## output directory
$cluster_file = $ARGV[3];	## cluster file by DomClust necessary for domain information

$sepChar = ':' if (! $sepChar); ## a separator between species name and gene name
#$sepCharAli = $sepChar;
$sepCharAli = '__' if (! $sepCharAli);

$ali_suff = ".aln";

$clust_type = 'amino'; ## ortholog clustering is based on amino acid sequence
$align_type = 'nucleotide'; ## alignment is based on nucleotide sequence

if (-f $spid_list) {
	&read_spid_list($spid_list);
}

die "Usage: seq_file align_dir output_dir\n" if (! $seq_file || !  $align_dir || ! $outdir);
mkdir $outdir if (! -d $outdir);


print STDERR "read converted sequences\n";
&read_seq($seq_file);

if (-f $cluster_file) {
	open(F, $cluster_file) || die;

	while(<F>){
		chomp;
		if (/^Cluster\s(\S+)/) {
			if ($clid) {
				&conv_align($clid, \%cluster);
				undef %cluster;
			}
			$clid = $1;
			next;
		}
		next if (/^#/);

		($domname, $from, $to) = split;
		($sp,$name,$dom) = &parse_name($domname, $sepChar);
		$domname =~ s/\(\d+\)//;

		$cluster{$domname} = {name=>$name, dom=>$dom, from=>$from, to=>$to};
	}
	close(F);
	if ($clid) {
		&conv_align($clid, \%cluster);
	}
} else {
	foreach $alif (<${align_dir}/*${ali_suff}>) {
		$clid = basename($alif, $ali_suff);
		&conv_align($clid);
	}
}


sub conv_align {
	my($clid, $cluster) = @_;
	open(ALI, "$align_dir/${clid}${ali_suff}") || die "Can't open $align_dir/${clid}${ali_suff}\n";
	my($seq, $name,$clst, $sqidx);
	my($aliseq);
	my(@align, @seqname);
	my($skip_flag);
#	print "## Clsuter $clid\n";
	$OutH = FileHandle->new(">$outdir/$clid.aln");
	my($out_domname, $domname);
	while(<ALI>){
		if (/^>(\S+)/) {
			if ($aliseq) {
				#	&output_seq($out_domname, $aliseq, $OutH);
				push(@align, $aliseq);
			}
			$domname = $1;
			$clst = $cluster->{$domname};
			my($sp, $name, $dom) = &parse_name($domname, $sepCharAli);
			$out_domname = &make_name($sp,$name,$dom,$sepCharAli);
			push(@seqname, $out_domname);
			if (%ConvSp) {
				$sp = $ConvSp{$sp};
				if (! $out_origname) {
					$out_domname = &make_name($sp,$name,$dom,$sepCharAli);
				}
			}
			$name = join($sepChar, $sp, $name);
			$seq = $Seq{$name};
			#print STDERR $clst->{from}," ", $clst->{to}, " ", $clst->{dom}, "; $clst\n" if ($clst->{dom});
			if (! $seq) {
				$skip_flag = 1;
			} else {
				$skip_flag = 0;
			}
			if ($clst->{dom}) {
				# extract domain
				my($seglen) = $clst->{to} - $clst->{from} + 1;
				if ($clust_type eq 'amino') {
					my($seqlen) = length($seq) / 3;
					#print STDERR "DOM:$seqlen;", $clst->{from}+$seglen-1, "; $seglen\n";
					if ($clst->{from} + $seglen == $seqlen) {
						## add stop codon
						#print STDERR "ADD STOP\n";
						$seglen ++;
					}
					$newseq = '-' x (($clst->{from} - 1) * 3);
					$newseq .= substr($seq, ($clst->{from} - 1) * 3, $seglen * 3);
					$seq = $newseq;
				} else {
					$seq = substr($seq, ($clst->{from} - 1), $seglen);
				}
			}
			$sqidx = 0;
			$aliseq = '';
		} elsif (! $skip_flag) {
			chomp;
			my($rseq) = $_;
			foreach $c (split(//, $rseq)) {
				if ($c eq '-') {
					if ($align_type eq 'amino') {
						$aliseq .= '---';
					} else {
						$aliseq .= '-';
					}
				} else {
					if ($align_type eq 'amino') {
						$aliseq .= substr($seq, $sqidx*3, 3);
					} else {
						$aliseq .= substr($seq, $sqidx, 1);
					}
					$sqidx++;
				}
			}
		}
	}
	close(ALI);
	if ($aliseq) {
		push(@align, $aliseq);
	}
	if (@align) {
		$align_new = &delete_allgap_col(\@align);
		#	$align_new = \@align;
		&output_align(\@seqname, $align_new, $OutH);
		##	&output_seq($out_domname, $aliout, $OutH);
	}
}

sub read_seq {
	my($filename) = @_;
	open(SEQ, $filename) || die "Can't open file: $filename\n";
	while(<SEQ>) {
		if (/^>\s*(\S+)/) {
			$name = $1;
			if ($Seq{$name}) {
				warn "seqence $name alreadly exists\n";
			}
		} else {
			chomp;
			s/\s//g;
			$Seq{$name} .= $_;
		}
	}
}

sub delete_allgap_col {
	my($align) = @_;
	my($maxlen) = 0;
	my(@gapnum, @totnum);
	my($i);

	foreach $aliseq (@{$align}) {
		for ($i = 0; $i < length($aliseq); $i++) {
			if (substr($aliseq, $i, 1) eq '-') {
				$gapnum[$i]++;
			}
			$totnum[$i]++;
		}
	}

	my(@new_align);
	for ($i = 0; $i < @totnum; $i++) {
		if ($gapnum[$i] == $totnum[$i]) {
			## all gap column
		} else {
			for (my $j = 0; $j < @{$align}; $j++) {
				$new_align[$j] .= substr($align->[$j], $i, 1);
			}
		}
	}
	\@new_align;
}
sub output_align {
	my($name, $align, $OutH) = @_;
	for ($i = 0; $i < @{$align}; $i++) {
		$aliout = $align->[$i];
		$out_domname = $name->[$i];
		&output_seq($out_domname, $aliout, $OutH);
	}
}
sub output_seq {
	my($name, $seq, $outH) = @_;
	my($LineLen) = 60;
	$outH->print(">$name\n");
	my($seqlen) = length($seq);
	for (my $i = 0; $i < $seqlen; $i+=$LineLen) {
		$outH->print(substr($seq, $i, $LineLen) . "\n");
	}
}

sub read_spid_list {
	my($file)= @_;
	open(F, $file) || die "Can't open file: $file\n";
	while(<F>){
		($sp,$spid,$name) = split;
		$ConvSp{$sp} = $name;
	}
	close(F);
}
sub parse_name {
	my($domname, $sep) = @_;
	$sep = $sepChar if (! $sep);
	my($sp,$name) = split($sep, $domname);
	my($dom);
	if ($name =~ /\((\d+)\)/) {
		$dom = $1;
		$name =~ s/\($dom\)//;
	}
	$sp, $name, $dom;

}
sub make_name {
	my($sp, $name, $dom, $sep) = @_;
	$sep = $sepChar if (! $sep);
	$domname = $sp . $sep . $name;
	$domname = "$domname($dom)" if ($dom);
	$domname;
}
