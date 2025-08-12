#!/usr/bin/perl -s

$coretab = $ARGV[0];
$seqfile = $ARGV[1];
$outdir = $ARGV[2];

$sepCharAli = '__' if (! $sepCharAli);

$ClustalCmd = "clustalo";

$MINLEN = 20;

if (! $outdir) {
	$outdir = "align_core";
}

mkdir $outdir if (! -d $outdir);

$ali_suff = ".out" if (! $ali_suff);

open(SEQ, $seqfile) || die;
while(<SEQ>){
	if (/^>(\S+)/) {
		$name = $1;
	} else {
		chomp;
		$Seq{$name} .= $_;
	}
}
close(SEQ);

open(F, $coretab) || die "Can't open $coretab\n";
while(<F>){
	chomp;
	if (/^Cluster (\S+)/) {
		if ($clustid) {
			&ext_seq($clustid, \%FoundList);
		}
		$clustid = $1;
		undef %FoundList;
	} elsif (/^#/) {
		next;
	} else {
		my($domname, $from, $to) = split;
		my($sp,$name,$dom) = &parse_name($domname);
		$spname = "${sp}:${name}";
		if ($FoundList{$spname}) {
			&addSegment($FoundList{$spname}, $from, $to, $dom);
		} else {
			$FoundList{$spname} = [ {from=>$from, to=>$to, dom=>$dom} ];
		}
		close(O);
	}
}
close(F);
if ($clustid) {
	&ext_seq($clustid, \%FoundList);
}

sub ext_seq {
	my($clustid, $foundList) = @_;
	my($orig_clustid) = $clustid;
	$orig_clustid =~ s/_\d+$//;

	my($seqfile) = "$outdir/$clustid.seq";
	my($alifile) = "$outdir/$clustid.aln";
	open(O, ">$seqfile") || die;
	foreach $spname (keys %{$foundList}) {
		foreach $seg (@{$foundList->{$spname}}) {
			if ($dom) {
				$seqseg = substr($Seq{$spname}, $seg->{from} - 1, $seg->{to} - $seg->{from} + 1);
			} else {
				$seqseg = $Seq{$spname};
			}
			next if (! $seqseg);
			if ($sepCharAli) {
				my($sp,$name) = split(/:/, $spname);
				$out_spname = join($sepCharAli, $sp, $name);
			} else {
				$out_spname = $spname;
			}
			print O ">$out_spname";
			print O " ($dom) [$from $to]" if ($dom);
			print O "\n";
			print O "$seqseg\n";
		}
	}
	close(O);
#	if ($exec_align) {
#		system("$ClustalCmd -i $seqfile -o $alifile");
#	}
}

sub parse_name {
	my($domname) = @_;
	my($sp,$name) = split(/:/, $domname);
	my($dom) = 0;
	if ($name =~ /\((\d+)\)/) {
		$dom = $1;
		$name =~ s/\($dom\)//;
	}
	$sp, $name, $dom;

}


sub addSegment {
	my($foundList, $from, $to, $dom) = @_;
	foreach $seg (@{$foundList}) {
		if ($from - $seg->{to} < $MINLEN) {
			$seg->{to} = $to;
		}
	}
	push(@{$foundList}, {from=>$from, to=>$to, dom=>$dom});
}
