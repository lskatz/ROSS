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
  GetOptions($settings,qw(numcpus=i readlength=i numreads=i tempdir=s help)) or die $!;
  $$settings{numcpus}||=1;
  $$settings{readlength}||=150;
  $$settings{numreads}||=1000;
  $$settings{tempdir}||=mktempdir();

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
  for my $infile(@$infiles){
    my $fh=openFastq($infile,$settings);
    while(my $read=<$fh>.<$fh>.<$fh>.<$fh>){
      push(@readEntry,$read);
    }
  }

  my $reads="";
  for(shuffle(@readEntry)){
    $reads.=$_;
  }
  return $reads;
}

sub usage{
  "$0: Generate random reads or shuffle reads
  Usage: $0 --numreads 100 --readlength 150
  "
}
