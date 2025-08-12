#!/usr/local/bin/perl

package ExtSequence;
#use Sequence::CodonTable;
use FileHandle;

sub new {
	my($class, $seq) = @_;
	my($this) = {};
	$this->{dir} = 1;
	bless $this, $class;
	if ($seq) {
		$this->setseq($seq);
	}
	$this;
}
sub setseq {
	my($this, $seq) = @_;
	$this->{seq} = $seq;
}
sub settype {
	my($this, $type) = @_;
	if ($type) {
		$this->{type} = $type;
	} else {
		$this->{type} = $this->infer_type;
	}
}
sub gettype {
	my($this) = @_;
	if (! $this->{type}) {
		$this->settype;
	}
	return $this->{type};
}
sub infer_type {
	my($this) = @_;
	return if (! $this->{seq});
	my %comp = $this->composition;
	my($nt);
	foreach $c (keys %comp) {
		if ($c =~ /[atgcu]/i) {
			$nt += $comp{$c}
		}
	}
	if ($nt / $this->length >= 0.8) {
		'nt';
	}  else {
		'aa';
	}
}
sub catseq {
	my($this, $seq) = @_;
	if (! ref($seq)) {
		$this->{seq} .= $seq;
	} elsif (ref($seq) eq 'Sequence') {
		$this->{seq} .= $seq->getseq;
	}
}
sub getseq {
	my($this) = @_;
	if ($this->{dir} < 0) {
		return &calc_complement($this->{seq});
	} else {
		return $this->{seq};
	}
}
sub length {
	my($this) = @_;
	$this->{length} = length($this->{seq});
	return $this->{length};
}
sub complement {
	my($this) = @_;
	$this->{dir} *= -1;
}
# 1-based coordinate
sub get_subseq {
	my($this, $from, $to, $dir) = @_;
	$this->get_subseq0($from-1, $to-$from+1, $dir);
}
sub get_subseq0 {
	my($this, $from, $len, $dir) = @_;
	my $subseq = substr($this->{seq}, $from, $len);
	if ($dir < 0) {
		$subseq = &calc_complement($subseq);
	}
	$subseq;
}

sub get_complement {
	my($this) = @_;
	&calc_complement($this->{$seq});
}
sub calc_complement {
	my( $org_seq ) = @_;
	my( $new_seq ) = $org_seq;

	$new_seq =~ tr/ACGTURYMKWSBDHVNXacgturymkwsbdhvnx/TGCAAYRKMWSVHDBNXtgcaayrkmwsvhdbnx/;
	$new_seq =~ tr/EFIJPQefijpq/FEJIQPfejiqp/;
	$new_seq =  reverse $new_seq;
	return $new_seq;
}

#sub translate {
#	my ($this, $transl_id) = @_;
#	if ($this->gettype ne 'nt') {
#		return;
#	}
#	return &translation($this->getseq, $transl_id);
#
#}
#sub translation {
#	my($seq, $transl_id) = @_;
#	my($trans_seq);
#	$transl_id = 1 if (! $transl_id);
#	my($codtab) = Sequence::CodonTable->new($transl_id);
#	my($length) = length($seq);
#	for (my $i= 0; $i < $length; $i+=3) {
#		my $codon = substr($seq, $i, 3);
#		$trans_seq .= $codtab->trans($codon);
#	} 
#	$trans_seq;
#}
sub composition {
	my($this) = @_;
	my(%Count);
	foreach $s (split(//, $this->{seq})) {
		$Count{$s}++;
	}
	return %Count;
}

sub print_fasta {
	my($this, $comment, $opt) = @_;
	my $LINELEN = 60;
	my $fh;
	my $seq = $this->getseq;
	if ($opt->{fh}) {
		$fh = $opt->{fh};
	} else {
		$fh = FileHandle->new(">&STDOUT");
	}
	$comment = $this->{comment} if (! $comment);
	if ($comment) {
		print $fh ">", $comment, "\n";
	}
	for (my $i = 0; $i < $this->length; $i+=$LINELEN) {
		print $fh substr($seq, $i, $LINELEN), "\n";
	}
}
1;
