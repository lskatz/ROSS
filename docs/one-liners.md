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

## Downsample reads to the same number of nucleotides as other reads

This problem came about because I truncated a 2x250bp fastq file to a 2x150bp fastq file.  However, I needed to backtrack and downsample the 2x250bp read set to the same coverage level as the new 2x150bp fastq file.  This uses an in-place bash trick to run `friends_rachel.pl` to find the number of nucleotides in the 150bp file.  It also uses a fast intermediate storage space in memory, in `/dev/shm/$USER`.

    mkdir /dev/shm/$USER
    ls reads/shuffled/*.fastq.gz | \
      xargs -P 24 -n 1 bash -c '
        b=$(basename $0); 
        friends_ursula.pl $0 --nobin --size $(friends_rachel.pl reads.150bp/shuffled/$b --fast --numcpus 1 | \
          tail -n 1 | \
          cut -f 3) | \
          gzip -c > /dev/shm/$USER/$b && \
        mv -nv /dev/shm/$USER/$b reads/downsampled/$b
      '
      
