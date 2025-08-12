#!/bin/sh

indir=$1
outdir=$2
if (( $# < 2 )); then
	echo Usage $0 indir [outdir]
	exit 1
fi
if [ "$outdir" = "" ]; then
	outdir=$indir
fi

if [ ! -d $outdir ]; then
	mkdir -p $outdir
fi

for f in $indir/*.seq; do
	clustid=`basename $f .seq`
	outf=$clustid.aln
	if [ -s $outdir/$outf ]; then
		continue
	fi
	echo $clustid
	clustalo --in $f --out $outdir/$outf
done

