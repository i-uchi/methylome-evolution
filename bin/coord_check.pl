#!/usr/bin/perl -s

use Bio::SeqIO;
use Bio::Tools::GFF;

$motfile = $ARGV[0];
$seqfile = $ARGV[1];

$seqi = Bio::SeqIO->new(-file=>$seqfile);
while ($seqobj = $seqi->next_seq) {
	$id = $seqobj->primary_id();
	if ($match_name) {
		$Seq{$id} = $seqobj;
	} else {
		push(@Seq, $seqobj);
	}
}

$ConvInfo = &get_shift($motfile);
&conv_data($motfile, $ConvInfo);

sub get_shift {
	my($motfile) = @_;
	my($moti) = Bio::Tools::GFF->new(-file=>$motfile, -gff_version=>3);
	my($shift, $reldir);
	my(%prev_shift);
	my(%ConvInfo);
	my($seqnum) = scalar(keys %Seq);
	my($info_num);

	while (my $feat = $moti->next_feature) {
		my($context) = $feat->get_tag_values('context');
		my($seqid) = $feat->seq_id;
		next if ($ConvInfo{$seqid});
		if ($match_name) {
			@Seq = ($Seq{$seqid});
		}
		my($strand) = $feat->strand;
		my($match_seq, @match);
		my($seqlen);
		foreach $seqobj (@Seq) {
			my($seq) = $seqobj->seq;
			$seqlen = $seqobj->length;
			@match = ($seq =~ m/$context/gi);
			my($match_dir) = 1;
			if (@match==0) {
				## no hit in forward seq check revcomp seq
				my($seqobj) = Bio::Seq->new(-seq=>$context);
				$rev_context = $seqobj->revcom->seq;
				@match = ($seq =~ m/$rev_context/gi);
				$match_dir = -1;
			}
#print "cont: $context, $strand\n";
#print "cont_rev: $rev_context\n";
			if (@match == 1) {
				$match_seq++;
				if ($match_dir < 0) {
					$seq =~ m/$rev_context/gi;
				} else {
					$seq =~ m/$context/gi;
				}
				my($pos) = pos $seq;
#print ">>>$pos\n";
				$pos -= ( (length($context) - 1) / 2 );	## central position in the context string

				my($origpos) = $feat->start;
				$reldir = $strand * $match_dir;

#if ($match_dir > 0) {
#print join("\t", $pos, substr($seq, $pos-1-20, 41), $feat->primary_tag, $strand * $reldir, length($context), $context),"<\n";
#} else {
##print join("\t", $pos, substr($seq, $pos-1, 1), $feat->primary_tag, $strand, length($rev_context), $rev_context),"<\n";
#print join("\t", $pos, substr($seq, $pos-1-20, 41), $feat->primary_tag, $strand * $reldir, length($rev_context), $context),"<\n";
#}

				if ($reldir < 0) {
					$origpos = &revpos($origpos, $seqlen);
				}
				$shift = ($pos - $origpos) % $seqlen;
				if ($DEBUG) {
					print "pos>>$pos, $origpos: $shift, $reldir; $strand, $match_dir\n";
				}
			}
		}
		if ($match_seq==1) {
			if ($prev_shift{$seqid}) {
				if  ($shift == $prev_shift{$seqid}) {
					$ConvInfo{$seqid} = {shift=>$shift, reldir=>$reldir, seqlen=>$seqlen};
					last if (++$info_num >= $seqnum);
				} else {
					print STDERR "Different shift: $shift, $prev_shift{$seqid}\n";
				}
			}
			$prev_shift{$seqid} = $shift;
		}
	}
	$moti->close();
	\%ConvInfo;
}
sub conv_data {
	my($motfile, $ConvInfo) = @_;
	my($moti) = Bio::Tools::GFF->new(-file=>$motfile, -gff_version=>3);
	while (my $feat = $moti->next_feature) {
		my($seqid) = $feat->seq_id;
		my($convInfo) = $ConvInfo->{$seqid};
		if (! $convInfo) {
			if (! $notFound{$seqid}) {
				print STDERR "sequence not found: $seqid\n";
				$notFound{$seqid} = 1;
			}
			next;
		}
		my($new_start) = $feat->start;
		my($new_end) = $feat->end;
		my($new_strand) = $feat->strand * $convInfo->{reldir};
		if ($convInfo->{reldir} < 0) {
			$new_start = &revpos($new_start, $convInfo->{seqlen});
			$new_end = &revpos($new_end, $convInfo->{seqlen});
		}
		$new_start = ($new_start + $convInfo->{shift}) % $convInfo->{seqlen};
		$new_end = ($new_end + $convInfo->{shift}) % $convInfo->{seqlen};
		$feat->start( $new_start );
		$feat->end( $new_end );
		$feat->strand($new_strand);
#print substr($Seq[0]->seq, $new_start-1, 5),"   $new_start, $new_end, $new_strand   ";
		print $feat->gff_string($moti), "\n";
	}
	$moti->close();
}
sub revpos {
	my($pos, $len) = @_;
	$len - $pos + 1;
}
