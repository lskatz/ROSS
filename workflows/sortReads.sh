#!/bin/bash

set -e

f=$1
r=$2

if [ "$r" == "" ]; then
  echo "$0: sorts reads by GC content in-place. This"
  echo "hypothetically reduces the gzip file size."
  echo
  echo "Usage: $0 f.fastq.gz r.fastq.gz"
  exit 1
fi;

# use ramdisk to go fast
tmpdir=/dev/shm/$USER
tmpF=$tmpdir/$(basename $f);
tmpR=$tmpdir/$(basename $r);
trap " { rm -f $tmpF $tmpR; } " EXIT

# 1. shuffle
# 2. each paired end goes on one line (eight joins in `paste`)
# 3. Calculate GC content and unshift that value on each line.
# 4. Sort by GC
# 5. Remove GC with `cut`
# 6. go back to eight lines per read pair
# 7. deshuffle to a temporary file
friends_joey.pl $f $r | \
  paste - - - - - - - - | \
  perl -F'\t' -lane '@GC=("$F[1]$F[5]"=~/[GCgc]/g); print join("\t",scalar(@GC)/length("$F[1]$F[5]"), @F);' | \
  sort -k1,1n | \
  cut -f2- | \
  tr '\t' '\n' | \
  friends_joey.pl -d -gz > $tmpF 2> $tmpR

# Move these reads back, so that this script edits in-place.
mv -v $tmpF $f 
mv -v $tmpR $r
