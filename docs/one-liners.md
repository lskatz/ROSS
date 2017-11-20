# One-liners

These are random one-liners or similar that make use of ROSS.
Some or many of them are large and so they will be displayed on multiple lines for readability.

## In-place fastq compression

This works by sorting by GC content, which plays into how the gzip algorithm works.
`friends_joey.pl` is used to shuffle the forward and reverse reads and then used to
deshuffle them.  The actual reads remain the same but are sorted differently, yielding
an estimated 10-30% compression gain.  Please see Torsten's blog for more details.  http://thegenomefactory.blogspot.com.au/2012/11/sorting-fastq-files-by-sequence-length.html?m=1

    friends_joey.pl Forward.fastq.gz Reverse.fastq.gz | \
      paste - - - - - - - - | \
      perl -F'\t' -lane '@GC=("$F[1]$F[5]"=~/[GCgc]/g); print join("\t",scalar(@GC)/length("$F[1]$F[5]"), @F);' | \
      sort -k1,1n | \
      cut -f2- | \
      tr '\t' '\n' | \
      friends_joey.pl -d -gz > 1.fastq.gz 2>2.fastq.gz

Test that the content hasn't changed by sorting and hash-summing

    zcat 1.fastq.gz       | sort | md5sum
    zcat Forward.fastq.gz | sort | md5sum
    zcat 2.fastq.gz       | sort | md5sum
    zcat Reverse.fastq.gz | sort | md5sum

