#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

use FindBin qw/$RealBin/;
use lib "$RealBin/../lib/perl5";
use ROSS qw/logmsg mktempdir/;

# http://perldoc.perl.org/perlop.html#Symbolic-Unary-Operators
# +Inf and -Inf will be a binary complement of all zeros
use constant MAXINT =>  ~0;
use constant MININT => -~0;

exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help numcpus=i tempdir=s kmerlength|kmer=i delta=i)) or die $!;
  $$settings{kmerlength} ||=21;
  $$settings{delta}      ||=100;
  $$settings{tempdir}    ||=mktempdir();
  $$settings{numcpus}    ||=1;

  my($fastq)=@ARGV;
  die usage() if(!$fastq || $$settings{help});
  die "ERROR: I could not find fastq at $fastq" if(!-e $fastq);

  # Count kmers and make a histogram.  kmerToHist()
  # will optionally write the histogram file.
  my $kmercount=countKmers($fastq,$$settings{kmerlength},$settings);
  $histogram=kmerToHist($kmercount,$settings);

  # Find peaks and valleys using the simplified
  # 'delta' algorithm.
  my $peaks=findThePeaksAndValleys($histogram,$$settings{delta},$settings);
  die Dumper $peaks

  # Configure the output
  my @outputPos=();
  push(@outputPos, @{$$peaks{valleys}}) if($$settings{valleys});
  push(@outputPos, @{$$peaks{peaks}})   if($$settings{peaks});
  @outputPos=sort {$$a[0] <=> $$b[0]} @outputPos;

  if(!@outputPos){
    logmsg "WARNING: no peaks or valleys were reported";
  }

  # Header for output table
  print join("\t",qw(kmer count))."\n";
  # Finish off the table of output.
  for my $position (@outputPos){
    print join("\t",@$position)."\n";
  }

  return 0;
}

sub findThePeaksAndValleys{
  my($hist, $delta, $settings)=@_;

  # http://perldoc.perl.org/perlop.html#Symbolic-Unary-Operators
  # +Inf and -Inf will be a binary complement of all zeros
  my($min,$max)=(-~0, ~0);
  my($minPos,$maxPos)=(0,0)
  my @maxTab=();
  my @minTab=();

  my $lookForMax=1;

  my $numZeros=0; # If we see too many counts of zero, then exit.

  for(my $kmerCount=$$settings{gt}+1;$kmerCount<@$hist;$kmerCount++){
    my $countOfCounts=$$hist[$kmerCount];
    if($countOfCounts == 0){
      $numZeros++;
    }
    if($countOfCounts > $max){
      $max=$countOfCounts;
      $maxPos=$kmerCount;
    }
    if($countOfCounts < $min){
      $min=$countOfCounts;
      $minPos=$kmerCount;
    }

    if($lookForMax){
      if($countOfCounts < $max - $delta){
        push(@maxTab,[$maxPos,$max]);
        $min=$countOfCounts;
        $minPos=$kmerCount;
        $lookForMax=0;
      }
    }
    else{
      if($countOfCounts > $min + $delta){
        push(@minTab,[$minPos,$min]);
        $max=$countOfCounts;
        $maxPos=$kmerCount;
        $lookForMax=1;
      }
    }

    last if($numZeros > 9999);
  }

  return {peaks=>\@maxTab, valleys=>\@minTab};
}

