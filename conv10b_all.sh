#!/bin/sh

indir=orig_data
conv_seqdir=seq10b
tmp_alidir=tmp_align
orig_alidir=orig_align
conv_alidir=ali10b

origseq_file=origseq.fna	## Genomic sequence file in FASTA format
annotation_file=origseq.gff	## Annotation file in GFF format
gbk_file=origseq.gbk	## GenBank format file (alternative to annotation file in GFF format)

orthogroup_file=core.sample.o0	## Ortholog group file in DomClust "o0" format
#orthogroup_file=core.all.o0

motlist_file=motlist_0.7.out	## motif list file listing all the motifs identified in either strain.

## cutoff
use_motif=1	# use motif to determine methylated bases
minr=0.7	# min coverage for significant motif 
min_score=50	# min score for modification identification
limit_cut=24	## consider common motifs only (motifs that appear more than this number of strains)

origseq=origseq.cds.fas
convseq=conv10b.cds.fas

if [ ! -d $indir ]; then
	echo "No input directory"
	exit
fi
mkdir -p $conv_seqdir

echo "## 1. Convert genome sequence"

for dir in $indir/*; do
	strain=`basename $dir`
	odir=$conv_seqdir/$strain
	mkdir -p $odir
	echo $strain

	if [ "$use_motif" == 1 -a -f $dir/motif.data ]; then
		## use motif file
		if [ -f $motlist_file ]; then
			bin/conv10b.pl -minr=$minr -limit_motif_file=$motlist_file -limit_cut=$limit_cut $dir/motif.data $dir/$origseq_file > $odir/conv10b.fas
		else
			bin/conv10b.pl -minr=$minr $dir/motif.data $dir/${origseq_file}  > $odir/conv10b.fas
		fi
	elif [ -f $dir/modifications.gff ]; then
		## use modification file
		bin/conv10b.pl -modfile -min_score=$min_score $dir/modifications.gff $dir/${origseq_file} > $odir/conv10b.fas
	fi
done

echo "## 2. Extract CDS"

for str_dir in $indir/*; do
	str_name=`basename $str_dir`
	echo $str_name
	if [ -f $str_dir/${annotation_file} ]; then
		bin/extseq_gff.pl -spname=$str_name $str_dir/${annotation_file} $str_dir/${origseq_file} > $str_dir/origseq.cds.fna
		c_str_dir=$conv_seqdir/$str_name
		for infas in $c_str_dir/conv10b*.fas; do
			if [[ "$infas" =~ cds ]]; then
#			       echo "skip CDS $infas"
				continue
			fi
			outfas=`echo $infas |sed s/.fas/.cds.fas/`
			bin/extseq_gff.pl -spname=$str_name $str_dir/${annotation_file} $infas > $outfas
		done
	elif [ -f $str_dir/${gbk_file} ]; then
		bin/extseq_gbk.pl -spname=$str_name $str_dir/${gbk_file}  > $str_dir/origseq.cds.fna

		for infas in $c_str_dir/conv10b*.fas; do
			if [[ "$infas" =~ cds ]]; then
#			       echo "skip CDS $infas"
				continue
			fi
			outfas=`echo $infas |sed s/.fas/.cds.fas/`
			bin/extseq_gbk.pl -spname=$str_name $str_dir/${gbk_file} $infas > $outfas
		done
	fi
done

cat $indir/*/origseq.cds.fna > $origseq
cat $conv_seqdir/*/conv10b.cds.fas > $convseq

echo "## 3. Create alignments"
bin/extract_core_seq.pl $orthogroup_file $origseq $tmp_alidir
bin/exec_align.sh $tmp_alidir $orig_alidir


echo "## 4. Create convered alignments"

bin/conv_align.pl $convseq $orig_alidir $conv_alidir $orthogroup_file
