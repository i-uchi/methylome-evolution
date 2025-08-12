#!/usr/bin/perl -s

use Bio::Seq;
use Bio::SeqIO;
use Bio::Tools::IUPAC;

$motfile = $ARGV[0];
$seqfile = $ARGV[1];

$minr = 0.7 if (! $minr);

$circular = 1;
$addlen = 15;

die "Usage: $0 motif seqfile\n" if (! $seqfile);

if ($DEBUG) {
	$VERBOSE = 1;
}

open(F, $motfile) || die "Can't open $motfile\n";
while(<F>) {
	chomp;
	next if (/^#/ || /^motif/i);
	($motif, $pos, $mod, $r_detect) = split(/\t/);
	$motif_regex = &to_regex($motif);
	$modbase = substr($motif, $pos-1, 1);
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
#	print "$motif_regex\n";
}
close(F);
if (@Motifs == 0) {
	print STDERR "No significant motif included: $motfile\n";
}

$seqio = Bio::SeqIO->new(-file=>$seqfile, -format=>fasta);
while ($seqobj = $seqio->next_seq()) {
	&conv_seq($seqobj);
}


sub to_regex {
	my($motif) = @_;
	my($motif_seq) = Bio::PrimarySeq->new(-seq => $motif, -alphabet => 'dna');
	my($motif_reg) = Bio::Tools::IUPAC->new(-seq=>$motif_seq)->regexp();
	$motif_reg;
}

sub conv_seq {
	my($seqobj) = @_;

	my($seq) = $seqobj->seq();
	my($seqname) = $seqobj->primary_id();

	my($revseq) = $seqobj->revcom()->seq();

	my($seqlen) = length($seq);
	my($newseq) = $seq;
	my(@Found);

	if ($circular) {
		my($addlen0) = ($addlen > $seqlen) ? $seqlen : $addlen;
		$seq .= substr($seq, 0, $addlen0);
		$revseq .= substr($revseq, 0, $addlen0);
	}
	foreach my $mot (@Motifs) {
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

	@Found = sort {$a->{methylpos} <=> $b->{methylpos}} @Found;
	foreach $motfnd (@Found) {
		$mot = $motfnd->{motdata};
		$methyl_ch = substr($seq, $motfnd->{methylpos} - 1, 1);
		## replace
		$origstr = substr($newseq, $motfnd->{methylpos} - 1, 1);
		if ($motfnd->{mod} eq 'm6A') {
			if ($DEBUG) {
				print STDERR "m6A>$origstr<$motfnd->{revcomp}\n";
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
			if (! $motfnd->{revcomp}) {
				substr($newseq, $motfnd->{methylpos} - 1, 1) = 'I';
			} else {
				substr($newseq, $motfnd->{methylpos} - 1, 1) = 'J';
			}
		} elsif ($motfnd->{mod} eq 'm5C') {
			if (! $motfnd->{revcomp}) {
				substr($newseq, $motfnd->{methylpos} - 1, 1) = 'P';
			} else {
				substr($newseq, $motfnd->{methylpos} - 1, 1) = 'Q';
			}
		} elsif ($motfnd->{mod} eq 'modified_base') {
		} else {
			print STDERR "Unrecognized modification $motfnd->{mod}, $methyl_ch\n";
		}
#		print join("\t", $motfnd->{methylpos}, $mot->{motif}, $methyl_ch),"\n";
	}
	$newseq =~ s/(.{60})/$1\n/g;
	print ">$seqname\n";
	print $newseq,"\n";
}
