#!/usr/bin/env -S perl -s
use Bio::Tools::GFF;
use Bio::SeqIO;
use ExtSequence;

$gbkfile = $ARGV[0];

$seqfile = $ARGV[1];

if ($seqfile) {
	open(F, $seqfile) || die "Can't open $seqfile\n";
	while(<F>){
		if (/^>\s*(\S+)/) {
			if ($seq) {
				$Seq{$name} = ExtSequence->new($seq);
			}
			$name = $1;
			$seq = '';
		} else {
			chomp;
			$seq .= $_;
		}
	}
	close(F);
	if ($seq) {
		$Seq{$name} = ExtSequence->new($seq);
	}
}

$seqio = Bio::SeqIO->new(-file => $gbkfile, -format=>'genbank') || die;
while ($seq = $seqio->next_seq) {
	$name = $seq->display_id;
	foreach $feat ($seq->get_SeqFeatures()) {
		next if ($feat->primary_tag ne 'CDS');
		if (%Seq) {
#			$name = $feat->seq_id;
#print STDERR ">$name,",$feat->start,"<\n";
			$subseq = $Seq{$name}->get_subseq($feat->start, $feat->end, $feat->strand);
		} else {
			$subseq_obj = $seq->trunc($feat->start, $feat->end);
			if ($feat->strand < 0) {
				$subseq_obj = $subseq_obj->revcom;
			}
			$subseq = $subseq_obj->seq;
		}
		($product) = $feat->get_tag_values('product');
		($locus_tag) = $feat->get_tag_values('locus_tag');
		if ($spname) {
			$locus_tag = "$spname:$locus_tag";
		}
		print ">$locus_tag $product";
		print " [", join(" ", $feat->start, $feat->end, $feat->strand), "]\n";
		print "$subseq\n";
	}
	
}

$seqio->close;
