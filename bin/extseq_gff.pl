#!/usr/bin/env -S perl -s
use FindBin;
use lib $FindBin::Bin;
use Bio::Tools::GFF;
use Bio::SeqIO;
use ExtSequence;

$gfffile = $ARGV[0];
$seqfile = $ARGV[1];

#$seqio = Bio::SeqIO->new(-file => $seqfile, -format=>'Fasta') || die;
#while ($seq = $seqio->next_seq) {
#	$name = $seq->display_id;
#	$Seq{$name} = $seq;
#}
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

$gffio = Bio::Tools::GFF->new(-file => $gfffile, -gff_version => 3);

while($feat = $gffio->next_feature()) {
	next if ($feat->primary_tag ne 'CDS');
	$name = $feat->seq_id;
	$subseq = $Seq{$name}->get_subseq($feat->start, $feat->end, $feat->strand);

#	$feat->attach_seq($Seq{$name});
#	$fseq = $feat->seq;
	($locus_tag) = $feat->get_tag_values('locus_tag');
	($product) = $feat->get_tag_values('product');
	if ($spname) {
		$locus_tag = "$spname:$locus_tag";
	}
	print ">$locus_tag $product";
	print " [", join(" ", $feat->start, $feat->end, $feat->strand), "]\n";
	print "$subseq\n";
#	print $fseq->seq,"\n";
}
$gffio->close;
exit;

