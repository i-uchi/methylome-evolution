#!/usr/bin/env -S perl -s

use Bio::Seq;
use Bio::SeqIO;
use Bio::Tools::IUPAC;
use Bio::Tools::GFF;

$motfile = $ARGV[0];
$seqfile = $ARGV[1];

$minr = 0.7 if (! $minr);
##$min_score = '' if (! $min_score); 	## default cutoff score seems to be 20

$seqfmt = 'genbank' if (! $seqfmt);

$circular = 1;
$addlen = 15;

if (! $seqfile) {
	die "Usage: $0 [-minr=#] [-modfile -min_score=# -skip_no_motif] [-limit_motif_file=#] motif seqfile\n";
}

if ($outfile) {
	open(O, ">$outfile") || die;
} else {
	open(O, ">&STDOUT") || die;
}

if ($DEBUG) {
	$VERBOSE = 1;
}

if ($modfile) {
	$FoundMotAll = &read_modfile($motfile);
	my(@k) = keys %$FoundMotAll;
	 if (@k== 1) {
		$FoundMotOne = $FoundMotAll->{ $k[0] };
	}
} else {
	my($limit_mot);
	if ($limit_motif_file) {
		$limit_mot = &read_limit_mot($limit_motif_file, $limit_cut);
	}
	$Motifs = &read_motifs($motfile, $limit_mot);

	if (@{$Motifs} == 0) {
		print STDERR "No significant motif included: $motfile\n";
	}
}

if ($seqfile =~ /\.(fasta|fas|fna|fa)$/) {
	$seqfmt = 'fasta';
}
$seqio = Bio::SeqIO->new(-file=>$seqfile, -format=>$seqfmt);
while ($seqobj = $seqio->next_seq()) {
	if ($Motifs) {
		$FoundMot = &search_motif($seqobj, $Motifs);
	} else {
		$seqid = $seqobj->display_id;
		$FoundMot = $FoundMotAll->{$seqid};
		if (! $FoundMot && $FoundMotOne) {
			## only one sequence
			$FoundMot = $FoundMotOne;
		}
	}
	&conv_seq($seqobj, $FoundMot);
}

sub read_motifs {
	my($motfile, $limit_hash) = @_;
	my(@Motifs);
	open(F, $motfile) || die "Can't open $motfile\n";
	while(<F>) {
		chomp;
		next if (/^#/ || /^motif/i);
		my($motif, $pos, $mod, $r_detect) = split(/\t/);
		if ($limit_hash && ! $limit_hash->{$motif}->{$pos}) {
			next;
		}
		my($motif_regex) = &to_regex($motif);
		my($modbase) = substr($motif, $pos-1, 1);
		if ($minr && $r_detect < $minr) {
			if ($VERBOSE) {
				print STDERR "Low detection rate: SKIP: $_\n";
			}
			next;
		}
		if ($mod eq 'unknown') {
			## skip unknown site
			if ($VERBOSE) {
				print STDERR "Unknown modification: SKIP: $_\n";
			}
			next;
		} elsif (! $mod) {
			if ($modbase =~ /A/i) {
				$mod = '6mA';
			} elsif ($modbase =~ /C/i) {
				$mod = '4mC';
			} else {
				print STDERR "Unrecognized modified base:  $mod\n";
			}
		}
		push(@Motifs, {motif_regex=>$motif_regex, motif=>$motif, len=>length($motif), pos=>$pos, mod=>$mod});
#		print "$motif_regex\n";
	}
	close(F);
	\@Motifs;
}
sub read_limit_mot {
	my($motlist, $limit_cut) = @_;
	my(%Limit);
	open(F, $motlist) || die;
	while(<F>) {
		next if (/^#/);
		my($motpos, $cnt) = split(/\t/);
		my($motif, $pos) = split(/,/, $motpos);
		if ($cnt >= $limit_cut) {
			$Limit{$motif}->{$pos}++;
		}
	}
	close(F);
	\%Limit;
}

sub read_modfile {
	my($modfile) = @_;
	my($gffi) = Bio::Tools::GFF->new(-file => $modfile, -gff_version=>3);
	my(%Found);
	while($feat = $gffi->next_feature()) {
		my($modtype) = $feat->primary_tag;
		## eliminate unknown modifications
		next if ($modtype !~ /m6A|m4C|m5C/);
#print STDERR "." if (++$cnt % 1000 == 0);

		my($seqid) = $feat->seq_id;
		my($start) = $feat->start;
		my($end) = $feat->end;
		my($strand) = $feat->strand;
		my($revcomp);
		$revcomp = ($feat->strand > 0) ? 0 : 1;
		if ($start != $end) {
			next;
		}
		my($methylpos) = $start;
		my($mot, $info) = ({}, {});
		if ($min_score && $feat->score < $min_score) {
			if ($VERBOSE) {
#				print STDERR "Low score: SKIP: ", $feat->gff_string($gffi), "\n";
			}
			next;
		}
		if ($skip_no_motif) {
			my($has_motif) = $feat->has_tag('motif');
			if (! $has_motif) {
				if ($VERBOSE) {
					print STDERR "No motif: SKIP: ", $feat->gff_string($gffi), "\n";
				}
				next
			}
		}
		if ($feat->has_tag('context')) {
			my($context) = $feat->get_tag_values('context');
			$info->{context} = $context;
		}
#print "methylpos=>$methylpos, mod=>$modtype, revcomp=>$revcomp\n";
		push(@{$Found{$seqid}}, {motdata=>$mot, seqid=>$seqid, methylpos=>$methylpos, mod=>$modtype, revcomp=>$revcomp, info=>$info});
	}
	 $gffi->close();
	\%Found;
}

sub to_regex {
	my($motif) = @_;
	my($motif_seq) = Bio::PrimarySeq->new(-seq => $motif, -alphabet => 'dna');
	my($motif_reg) = Bio::Tools::IUPAC->new(-seq=>$motif_seq)->regexp();
	$motif_reg;
}

sub search_motif {
	my($seqobj, $Motifs) = @_;

	my($seq) = $seqobj->seq();
	my($seqname) = $seqobj->primary_id();

	my($revseq) = $seqobj->revcom()->seq();

	my($seqlen) = length($seq);
	my(@Found);

	if ($circular) {
		my($addlen0) = ($addlen > $seqlen) ? $seqlen : $addlen;
		$seq .= substr($seq, 0, $addlen0);
		$revseq .= substr($revseq, 0, $addlen0);
	}
	foreach my $mot (@{$Motifs}) {
#print ">$seq; $seqlen\n";
		while ($seq =~ /$mot->{motif_regex}/g) {
			$hitpos = (pos $seq) - $mot->{len} + 1;
			$methylpos = $hitpos + $mot->{pos} - 1;

			last if ($hitpos >= $seqlen);
			$methylpos = ($methylpos - 1) % $seqlen + 1;

			push(@Found, {motdata=>$mot, methylpos=>$methylpos, hitpos=>$hitpos, mod=>$mot->{mod}, revcomp=>0});

			if ($DEBUG) {
				$hitseq = substr($seq, $hitpos-1, $mot->{len});
				$methyl_str = substr($seq, $methylpos-1, 1);
				print STDERR "HitFwd: $seqname: $methylpos $hitpos $mot->{motif} $mot->{mod}: $hitseq; $methyl_str\n";
		    	}
		}
#print ">$revseq\n";
		while ($revseq =~ /$mot->{motif_regex}/g) {
			$hitpos_r = (pos $revseq) - $mot->{len} + 1;
			$methylpos_r = $hitpos_r + $mot->{pos} - 1;
			$hitpos = $seqlen - $hitpos_r + 1;
			$methylpos = $seqlen - $methylpos_r + 1;

			next if ($hitpos >= $seqlen || $hitpos < 0);
			$methylpos = ($methylpos - 1) % $seqlen + 1;


			push(@Found, {motdata=>$mot, methylpos=>$methylpos, hitpos=>$hitpos, mod=>$mot->{mod}, revcomp=>1});

			if ($DEBUG) {
##				$p = pos $revseq;
				$methyl_str = substr($revseq, $methylpos_r-1, 1);
				$methyl_str2 = substr($seq, $methylpos-1, 1);
				print STDERR "HitRev: $seqname: $methylpos($methylpos_r) $hitpos($hitpos_r) $mot->{motif} $mot->{mod}; $methyl_str,$methyl_str2\n";
			}
		}
	}
	\@Found;
}

sub conv_seq {
	my($seqobj, $Found) = @_;
	my($seq) = $seqobj->seq();
	$seq = uc($seq);
	my($seqname) = $seqobj->display_id();
	my($newseq) = $seq;
	my(@Found) = sort {$a->{methylpos} <=> $b->{methylpos}} @{$Found};
	my($errflag);
	foreach $motfnd (@Found) {
		if ($check_seqid && $seqname ne $motfnd->{seqid}) {
			## data for a different sequence
#			print $seqname, " ", $motfnd->{seqid},"\n";
			next;
		}
		my($mot) = $motfnd->{motdata};
		$methyl_ch = substr($seq, $motfnd->{methylpos} - 1, 1);
		## replace
		my($origstr) = substr($seq, $motfnd->{methylpos} - 1, 1);
		my($rev) = $motfnd->{revcomp} ? '(rev)' : '';
		if ($DEBUG) {
			my($tmpch) = substr($seq, $motfnd->{methylpos} - 1, 1);
			my($tmpstr) =  substr($seq, $motfnd->{methylpos} - 21, 41);
			print STDERR ">", join("\t", $motfnd->{methylpos}, $motfnd->{revcomp}, $tmpch, $tmpstr, $motfnd->{info}->{context}),"<\n";
		}
		if ($motfnd->{mod} eq 'm6A') {
			if ($DEBUG) {
				print STDERR "m6A>$origstr<$motfnd->{revcomp}\n";
			}
			if ( ($motfnd->{revcomp} == 0 && $origstr ne 'A') ||
				($motfnd->{revcomp} == 1 && $origstr ne 'T') ) {
				print STDERR "Illegal base: $motfnd->{mod}$rev: $motfnd->{methylpos}: $origstr; >$motfnd->{hitpos};$mot->{motif_regex}<\n";
				$errflag++;
			}
			if (! $motfnd->{revcomp}) {
				substr($newseq, $motfnd->{methylpos} - 1, 1) = 'E';
			} else {
				substr($newseq, $motfnd->{methylpos} - 1, 1) = 'F';
			}
		} elsif ($motfnd->{mod} eq 'm4C') {
			if ($DEBUG) {
				print STDERR "m4C>$origstr<$motfnd->{revcomp}\n";
			}
			if ( ($motfnd->{revcomp} == 0 && $origstr ne 'C') ||
				($motfnd->{revcomp} == 1 && $origstr ne 'G') ) {
				print STDERR "Illegal base: $motfnd->{mod}$rev: $motfnd->{methylpos}: $origstr\n";
				$errflag++;
			}
			if (! $motfnd->{revcomp}) {
				substr($newseq, $motfnd->{methylpos} - 1, 1) = 'I';
			} else {
				substr($newseq, $motfnd->{methylpos} - 1, 1) = 'J';
			}
		} elsif ($motfnd->{mod} eq 'm5C') {
			if ($DEBUG) {
				print STDERR "m4C>$origstr<$motfnd->{revcomp}\n";
			}
			if ( ($motfnd->{revcomp} == 0 && $origstr ne 'C') ||
				($motfnd->{revcomp} == 1 && $origstr ne 'G') ) {
				print STDERR "Illegal base: $motfnd->{mod}$rev: $motfnd->{methylpos}: $origstr\n";
				$errflag++;
			}
			if (! $motfnd->{revcomp}) {
				substr($newseq, $motfnd->{methylpos} - 1, 1) = 'P';
			} else {
				substr($newseq, $motfnd->{methylpos} - 1, 1) = 'Q';
			}
		} elsif ($motfnd->{mod} eq 'modified_base') {
			if ($VERBOSE) {
				print STDERR "Unknown modification $motfnd->{mod}, $methyl_ch\n";
			}
		} else {
			print STDERR "Unrecognized modification $motfnd->{mod}, $methyl_ch\n";
		}
#		print join("\t", $motfnd->{methylpos}, $mot->{motif}, $methyl_ch),"\n";
		last if ($errflag > 10);
	}
	if ($errflag) {
		unlink($outfile) if (-f $outfile);
		die "Error found: $errflag: motfile=$motfile, seqfile=$seqfile\n";
	}
	$newseq =~ s/(.{60})/$1\n/g;
	print O ">$seqname\n";
	print O $newseq,"\n";
}
