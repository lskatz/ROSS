#!/usr/bin/env perl
# Kmer counter with no dependencies.
# I wanted to increase capatibility with Perl and have a more
# standalone script
# Author: Lee Katz <lkatz@cdc.gov>

use strict;
use warnings;
use Getopt::Long qw/GetOptions/;
use Data::Dumper qw/Dumper/;
use File::Basename qw/basename fileparse/;
use File::Temp qw/tempdir/;
use List::Util qw/max/;

use FindBin qw/$RealBin/;
use lib "$RealBin/../lib/perl5";
use FRIENDS qw/logmsg mktempdir openFastq/;

use threads;
use Thread::Queue;

# http://perldoc.perl.org/perlop.html#Symbolic-Unary-Operators
# +Inf and -Inf will be a binary complement of all zeros
use constant MAXINT =>  ~0;
use constant MININT => -~0;

local $0=basename $0;

exit main();

sub main{
  my $settings={};
  GetOptions($settings,qw(help sortby=s kmerlength|kmer=i gt=i kmerCounter=s delta=i tempdir=s numcpus=i)) or die $!;
  $$settings{kmerlength} ||=21;
  $$settings{kmerCounter}||="";
  $$settings{delta}      ||=100;
  $$settings{gt}         ||=1;
  $$settings{tempdir}    ||=mktempdir();
  $$settings{numcpus}    ||=1;
  $$settings{sortby}     ||="";

  my($fastq)=@ARGV;
  die usage() if(!$fastq || $$settings{help});
  die "ERROR: I could not find fastq at $fastq" if(!-e $fastq);

  # Count kmers and make a histogram.  kmerToHist()
  # will optionally write the histogram file.
  my $kmercount=countKmers($fastq,$$settings{kmerlength},$settings);

  # Print the output
  if($$settings{sortby} =~ /kmer/i){
    for my $kmer(sort {$a cmp $b} keys(%$kmercount)){
      print $kmer."\t".$$kmercount{$kmer}."\n";
    }
  } elsif($$settings{sortby} =~ /count/i){
    for my $kmer(sort {$$kmercount{$a} <=> $$kmercount{$b}} keys(%$kmercount)){
      print $kmer."\t".$$kmercount{$kmer}."\n";
    }
  } else {
    while(my($kmer,$count)=each(%$kmercount)){
      print $kmer."\t".$count."\n";
    }
  }

  return 0;
}

# TODO: count kmers with faster programs in this order of
# priority: jellyfish, KAnalyze, Khmer
# and lastly, pure perl.
sub countKmers{
  my($fastq,$kmerlength,$settings)=@_;

  my $kmerHash={};

  if($$settings{kmerCounter} =~ /^\s*$/){
    logmsg "Auto-detecting which kmer counter program to use.";
    if(which("jellyfish")){
      $$settings{kmerCounter}="jellyfish";
    } else {
      $$settings{kmerCounter}="pureperl";
    }
  }

  # try/catch kmer counting and if it fails, do the
  # pure perl method.  Don't redo pure perl if
  # it was already set up that way.
  eval{
    if($$settings{kmerCounter}=~ /(pure)?.*perl/i){
      logmsg "Counting with pure perl";
      $kmerHash=countKmersPurePerl($fastq,$kmerlength,$settings);
    } elsif($$settings{kmerCounter} =~ /jellyfish/i){
      logmsg "Counting with Jellyfish";
      $kmerHash=countKmersJellyfish($fastq,$kmerlength,$settings);
    } else {
      die "ERROR: I do not understand the kmer counter $$settings{kmerCounter}";
    }
  };

  if($@){
    if($$settings{kmerCounter}=~ /(pure)?.*perl/i){
      die "ERROR counting kmers with pure perl: $!";
    }
    logmsg "Error detected.  Trying again by counting with pure perl";
    $kmerHash=countKmersPurePerl($fastq,$kmerlength,$settings);
  }

  return $kmerHash;
}

sub countKmersPurePerl{
  my($fastq,$kmerlength,$settings)=@_;

  # Multithreading
  my $seqQ=Thread::Queue->new;
  my @thr;
  for(0..$$settings{numcpus}-1){
    $thr[$_]=threads->new(\&countKmersPurePerlWorker,$kmerlength,$seqQ,$settings);
  }

  # Pure perl to make this standalone... the only reason
  # we are counting kmers in Perl instead of C.
  my $fastqFh=openFastq($fastq,$settings);
  my $i=0;
  my @buffer=();
  while(<$fastqFh>){ # burn the read ID line
    $i++;
    my $seq=<$fastqFh>;
    push(@buffer, $seq);

    if($i % 1000000 == 0){
      logmsg "Enqueuing ".scalar(@buffer)." reads for kmer counting";
      $seqQ->enqueue(@buffer);
      @buffer=();
    }
    
    # Burn the quality score lines
    <$fastqFh>;
    <$fastqFh>;
  }
  close $fastqFh;

  logmsg "Enqueuing ".scalar(@buffer)." reads for kmer counting";
  $seqQ->enqueue(@buffer);
  
  # Send the termination signal
  $seqQ->enqueue(undef) for(@thr);

  while($seqQ->pending > @thr){
    logmsg $seqQ->pending . " reads still pending...";
    for(1..60){
      last if($seqQ->pending <= @thr);
      sleep 1;
    }
  }

  # Join the threads and put everything into a large kmer hash
  logmsg "Combining kmers from multiple threads";
  my %kmer=();
  for(@thr){
    my $threadKmer=$_->join;
    for my $kmer(keys(%$threadKmer)){
      $kmer{$kmer}+=$$threadKmer{$kmer};
    }
  }

  # Filtering step
  logmsg "Filtering";
  while(my($kmer,$count)=each(%kmer)){
    delete($kmer{$kmer}) if($count < $$settings{gt});
  }

  return \%kmer;
}

sub countKmersPurePerlWorker{
  my($kmerlength,$seqQ,$settings)=@_; 

  my %kmer;
  while(defined(my $seq=$seqQ->dequeue)){

    my $numKmersInRead=length($seq)-$kmerlength;

    # Count kmers in a sliding window.
    # We must keep this loop optimized for speed.
    for(my $j=0;$j<$numKmersInRead;$j++){
      $kmer{substr($seq,$j,$kmerlength)}++;
    }

  }

  return \%kmer;
}


sub countKmersJellyfish{
  my($fastq,$kmerlength,$settings)=@_;
  my $basename=basename($fastq);
  my %kmer=();

  # Version checking
  my $jfVersion=`jellyfish --version`;
  # e.g., jellyfish 2.2.6
  if($jfVersion =~ /(jellyfish\s+)?(\d+)?/){
    my $majorVersion=$2;
    if($majorVersion < 2){
      die "ERROR: Jellyfish v2 or greater is required";
    }
  }
  
  my $outprefix="$$settings{tempdir}/$basename.mer_counts";
  my $jfDb="$$settings{tempdir}/$basename.merged.jf";
  my $kmerTsv="$$settings{tempdir}/$basename.jf.tsv";

  # Counting
  logmsg "Counting kmers in $fastq";
  my $jellyfishCountOptions="-s 10000000 -m $kmerlength -o $outprefix -t $$settings{numcpus}";
  my $uncompressedFastq="$$settings{tempdir}/$basename.fastq";
  if($fastq=~/\.gz$/i){
    logmsg "Decompressing fastq for jellyfish into $uncompressedFastq";
    system("zcat $fastq > $uncompressedFastq"); die if $?;
    system("jellyfish count $jellyfishCountOptions $uncompressedFastq");
  } else {
    system("jellyfish count $jellyfishCountOptions $fastq");
  }
  die "Error: problem with jellyfish" if $?;

  logmsg "Jellyfish dump $outprefix";
  my $lowerCount=$$settings{gt}-1;
  system("jellyfish dump --lower-count=$lowerCount --column --tab -o $kmerTsv $outprefix");
  die if $?;

  # Load kmers to memory
  logmsg "Reading jellyfish kmers to memory";
  open(TSV,$kmerTsv) or die "ERROR: Could not open $kmerTsv: $!";
  while(<TSV>){
    chomp;
    my @F=split /\t/;
    $kmer{$F[0]}=$F[1];
  }
  close TSV;

  # cleanup
  for($jfDb, $kmerTsv, $outprefix, $uncompressedFastq){
    unlink $_ if($_ && -e $_);
  }

  return \%kmer;
}

# http://www.perlmonks.org/?node_id=761662
sub which{
  my($exe,$settings)=@_;
  
  my $tool_path="";
  for my $path ( split /:/, $ENV{PATH} ) {
      if ( -f "$path/$exe" && -x "$path/$exe" ) {
          $tool_path = "$path/$exe";
          last;
      }
  }
  
  return $tool_path;
}

sub usage{
  "
  $0: 
    Print a tab-delimited count of kmers from a reads file

  Usage: $0 file.fastq[.gz]
  --kmer     21     kmer length
  --delta    100    How different the counts have to be to
                    detect a valley or peak
  --numcpus  1

  MISC
  --kmerCounter ''  The kmer counting program to use.
                    Default: (empty string) auto-choose
                    Options: perl, jellyfish
  --gt       1      Print kmers greater than this count.
  --sortby   ''     Default: unsorted. Options: kmer, count
  "
}

