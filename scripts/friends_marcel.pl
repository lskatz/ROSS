#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Statistics::Descriptive;
use POSIX qw/floor ceil/;
use File::Basename qw/basename fileparse/;

use FindBin qw/$RealBin/;
use lib "$RealBin/../lib/perl5";
use Friends qw/logmsg mktempdir openFastq/;
use Kmer;

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
  logmsg "Counting kmers";
  my $kmer=Kmer->new($fastq,{kmerlength=>$$settings{kmerlength},numcpus=>$$settings{numcpus}});
  my $kmercount=$kmer->count;
  my $histogram=$kmer->histogram;

  #logmsg "Histogram:\n".join("\n",@$histogram);

  # Find peaks and valleys using the simplified
  # 'delta' algorithm.
  logmsg "Finding peaks and valleys in kmer coverage";
  my $peaks=findThePeaksAndValleys($histogram,$$settings{delta},$settings);
  logmsg "Peaks   at k={".join(",",map{$$_[0]} @{ $$peaks{peaks} })."}";
  logmsg "Valleys at k={".join(",",map{$$_[0]} @{ $$peaks{valleys} })."}";

  # re-score reads according to kmer
  logmsg "Determining score adjustments";
  my $kmerCountScore=scoreKmerCounts($peaks,$kmercount,$histogram,$settings);
  logmsg "Rescoring reads";
  my $reads=rescoreReads($fastq,$kmerCountScore,$kmercount,$settings);

  print $reads;

  return 0;
}

sub scoreKmerCounts{
  my($peaks,$kmercount,$histogram,$settings)=@_;

  # From valley to the next valley
  my @kmerCountScore=(0) x scalar(@$histogram);
  for(my $peakCount=0;$peakCount<scalar(@{$$peaks{peaks}});$peakCount++){
    
    # Define the peaks as being from valley to valley.
    my $firstValley=0;
    if($peakCount > 0){
      $firstValley=$$peaks{valleys}[$peakCount-1][0];
    }
    my $firstPeak=$$peaks{peaks}[$peakCount][0];
    my $secondValley=$$peaks{valleys}[$peakCount][0];
    if(!defined $secondValley){
      $secondValley=scalar(@$histogram)-1;
    }

    # The peak has been defined, so find the mean and stdev
    my $stat=Statistics::Descriptive::Full->new();
    for(my $j=$firstValley;$j<=$secondValley;$j++){
      next if($j==0);
      while(my($kmer,$count)=each(%$kmercount)){
        $stat->add_data($count) if($count==$j);
      }
    }

    # Score kmer counts:
    #   -1 if it is outside of 2 stdev
    #   +1 if it is inside  of 1 stdev
    for(my $j=$firstValley;$j<$secondValley;$j++){
      if($j < floor($stat->mean - 2*$stat->standard_deviation)){
        $kmerCountScore[$j]=-1;
      } elsif($j >= floor($stat->mean - $stat->standard_deviation) && $j <= ceil($stat->mean + 1*$stat->standard_deviation)){
        $kmerCountScore[$j]= 1;
      } elsif($j > ceil($stat->mean + 2*$stat->standard_deviation)){
        $kmerCountScore[$j]=-1;
      }
    }
  }

  return \@kmerCountScore;
}

sub findThePeaksAndValleys{
  my($hist, $delta, $settings)=@_;

  # http://perldoc.perl.org/perlop.html#Symbolic-Unary-Operators
  # +Inf and -Inf will be a binary complement of all zeros
  my($min,$max)=($$hist[0], $$hist[0]);
  my($minPos,$maxPos)=(0,0);
  my @maxTab=();
  my @minTab=();

  my $lookForMax=1;

  my $numZeros=0; # If we see too many counts of zero, then exit.

  for(my $kmerCount=0;$kmerCount<@$hist;$kmerCount++){
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

sub rescoreReads{
  my($fastq,$kmerCountScore,$kmercount,$settings)=@_;

  my $kmerlength=$$settings{kmerlength};

  my $updatedReads="";

  my $fh=openFastq($fastq);
  while(my $id=<$fh>){
    my $sequence=<$fh>;
    my $plus=<$fh>;
    my $qual=<$fh>;

    chomp($id,$sequence,$plus,$qual);
    
    my $numKmers=length($sequence)-$kmerlength;
    my @qual=map{ ord($_)-33 } split(//,$qual);
    for(my $i=0;$i<$numKmers;$i++){
      # add the relative score
      for(my $j=$i;$j<$i+$kmerlength;$j++){
        $qual[$j]+=$$kmerCountScore[ $$kmercount{ substr($sequence,$i,$kmerlength) } ];
      }
    }

    # After rescoring, enforce some min/max limits
    for(@qual){
      if($_ > 40){
        $_=40;
      }elsif($_ < 0){
        $_=0;
      }
    }

    $updatedReads.="$id\n$sequence\n$plus\n".join("",map{ chr($_+33) } @qual)."\n";
  }

  return $updatedReads;
}

sub usage{
  local $0=basename $0;
  "$0: alters quality scores of a fastq file based on kmer abundance
  Usage: $0 file.fastq

  Steps: count kmers
         make histogram
         Find areas between two 'valleys' of kmer abundance
         Calculate mean/stdev for each area
         For areas within 1 stdev of a mean, add +1 score for each base's phred score
         For areas further than 2 stdevs of a mean, subtract 1.
         Traverse through each read to uncover kmers and make the score readjustment
  "
}
