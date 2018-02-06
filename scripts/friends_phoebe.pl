#!/usr/bin/env perl

# TODO use ART for advanced things
# TODO allow for basic ratios of ATCGN
# TODO give better odds to certain read scores
# TODO PE-interleaved awareness

use strict;
use warnings;
use File::Basename qw/basename/;
use Getopt::Long;
use POSIX qw/floor/;
use List::Util qw/shuffle/;
use Data::Dumper;

use FindBin qw/$RealBin/;
use lib "$RealBin/../lib/perl5";
use Friends qw/logmsg mktempdir openFastq/;

use threads;
use Thread::Queue;

local $0=basename $0;
my @nt=qw(A T C G);
my $numNt=@nt;

exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(pe|paired-end numcpus=i readlength=i numreads=i numbases=i tempdir=s help)) or die $!;
  $$settings{numcpus}||=1;
  $$settings{readlength}||=150;
  $$settings{tempdir}||=mktempdir();
  $$settings{numbases}||=0;
  $$settings{numreads}||=0;

  die usage() if($$settings{help});
  
  my $reads;
  if(@ARGV){
    $reads=randomizeReads(\@ARGV,$settings);
  } 
  # If no reads are given, then just shuffle
  else {
    $reads=generateRandomReads($$settings{readlength},$settings);
  }
  print $reads;

  return 0;
}

sub initialQualScores{
  my @qual;
  for(0..40){
    push(@qual,chr($_+33));
  }
  return @qual;
}

sub generateRandomReads{
  my($readlength,$settings)=@_;

  # If number of bases is specified but not number of reads,
  # set the number of reads to the max integer.
  if($$settings{numbases} > 1 && $$settings{numreads}==0){
    $$settings{numreads} = ~0; # max int
  } 
  # Same thing happens to the number of bases if the number
  # of reads have been set
  elsif( $$settings{numreads} > 1 && $$settings{numbases}==0){
    $$settings{numbases} = ~0;
  } 
  # If neither are set, then just set them both to 1.
  elsif( $$settings{numreads}==0 && $$settings{numbases}==0){
    $$settings{numreads}=1;
    $$settings{numbases}=1;
  }

  my @qual=initialQualScores();
  my $numQual=@qual;

  my $readStr="";
  my $basesCount=0;
  for(my $i=0;$i<$$settings{numreads};$i++){
    my $seq  ="";
    my $qual ="";
    my $seq2 ="";
    my $qual2="";
    for(my $j=0;$j<$readlength;$j++){
      $seq  .=$nt[ floor(rand(4)) ];
      $qual .=$qual[ floor(rand($numQual)) ];
      $seq2 .=$nt[ floor(rand(4)) ];
      $qual2.=$qual[ floor(rand($numQual)) ];
    }
    my $id =$i;
    my $id2=$i;
    if($$settings{pe}){
      $id ="$i/1";
      $id2="$i/2";
    }

    $readStr.='@read'.$id."\n".$seq."\n"."+\n".$qual."\n";
    $basesCount+=length($seq);

    if($$settings{pe}){
      $readStr.='@read'.$id2."\n".$seq2."\n"."+\n".$qual2."\n";
    $basesCount+=length($seq2);
    }

    last if($basesCount > $$settings{numbases});
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

      if($$settings{pe}){
        my $id2=<$fh>;
        my $sequence2=<$fh>;
        my $plus2=<$fh>;
        my $qual2=<$fh>;
        $read.=$id2.$sequence2.$plus2.$qual2;
        # keep track of the number of bases
        chomp($sequence2);
        $sequence.=$sequence2;
      }
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
         $0 --numreads 100 --numbases 9999 --readlength 150

         # A single random read
         $0 

         # Get reads from a real fastq file; do not truncate.
         $0 --numreads 100 file.fastq[.gz]

         # 1x coverage for many bacterial genomes
         $0 --numbases 5000000 file.fastq[.gz]

  OPTIONS
  --pe              paired end
  --numreads  0     Maximum number of reads.
  --numbases  0     Maximum number of nucleotides.
  "
}

