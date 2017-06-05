#!/usr/bin/env perl

# TODO use ART for advanced things
# TODO allow for basic ratios of ATCGN
# TODO give better odds to certain read scores

use strict;
use warnings;
use File::Basename qw/basename/;
use Getopt::Long;
use POSIX qw/floor/;
use List::Util qw/shuffle/;

use FindBin qw/$RealBin/;
use lib "$RealBin/../lib/perl5";
use Friends qw/logmsg mktempdir openFastq/;

use threads;
use Thread::Queue;

local $0=basename $0;
my @nt=qw(A T C G);
my $numNt=@nt;
my @qual=();
my $numQual=@qual;

exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(numcpus=i readlength=i numreads=i numbases=i tempdir=s help)) or die $!;
  $$settings{numcpus}||=1;
  $$settings{readlength}||=150;
  $$settings{numreads}||=1000;
  $$settings{tempdir}||=mktempdir();
  $$settings{numbases}||=0;

  @qual=initialQualScores();

  die usage() if($$settings{help});
  
  my $reads;
  if(@ARGV){
    $reads=randomizeReads(\@ARGV,$settings);
  } 
  # If no reads are given, then just shuffle
  else {
    $reads=generateRandomReads($$settings{numreads},$$settings{readlength},$settings);
  }
  print $reads;

  return 0;
}

sub initialQualScores{
  my @qual;
  for(0..40){
    push(@qual,chr($_+33));
  }
  $numQual=@qual;
  return @qual;
}

sub generateRandomReads{
  my($numreads,$readlength,$settings)=@_;

  my $readStr="";
  for(my $i=0;$i<$numreads;$i++){
    my $seq ="";
    my $qual="";
    for(my $j=0;$j<$readlength;$j++){
      $seq .=$nt[ floor(rand(4)) ];
      $qual.=$qual[ floor(rand($numQual)) ];
    }

    $readStr.='@read'.$i."\n".$seq."\n"."+\n".$qual."\n";
  }

  return $readStr;
}

sub randomizeReads{
  my($infiles,$settings)=@_;

  my @readEntry;
  my %basesPerRead;
  for my $infile(@$infiles){
    my $fh=openFastq($infile,$settings);
    while(my $id=<$fh>){
      my $sequence=<$fh>;
      my $plus=<$fh>;
      my $qual=<$fh>;
      my $read=$id.$sequence.$plus.$qual;
      push(@readEntry,$read);

      # Record the number of bases per read
      if($$settings{numbases}){
        chomp($sequence);
        $basesPerRead{$read}=length($sequence);
      }
    }
  }

  my $reads="";
  my $numreads=0;
  my $numBases=0;
  for(shuffle(@readEntry)){
    $reads.=$_;
    last if(++$numreads == $$settings{numreads});
    
    # Wrap it up if we have enough bases
    if($$settings{numbases}){
      $numBases+=$basesPerRead{$_};
      last if($numBases > $$settings{numbases});
    }
  }
  return $reads;
}

sub usage{
  "$0: Generate random reads or shuffle reads
  Usage: 
         # Completely random reads
         $0 --numreads 100 --readlength 150

         # Get reads from a real fastq file; do not truncate.
         $0 --numreads 100 file.fastq[.gz]

         # 1x coverage for many bacterial genomes
         $0 --numbases 5000000 file.fastq[.gz]
  "
}

