#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use File::Basename qw/basename/;
use File::Temp qw/tempdir/;
use Getopt::Long qw/GetOptions/;
use Array::IntSpan;

sub logmsg{ print STDERR basename($0).": @_\n"; }
exit main();

sub main{
  my $settings={};
  GetOptions($settings,qw(totalNs=i tempdir=s numcpus=i posFrequencies=s nsPerRead=s help)) or die $!;
  die usage() if($$settings{help});
  $$settings{totalNs}||=0;
  $$settings{tempdir}||=tempdir("charles.XXXXXX", TMPDIR=>1, CLEANUP=>1);
  $$settings{numcpus}||=1;

  # Required settings
  for (qw(posFrequencies nsPerRead)){
    $$settings{$_}||=die("ERROR: need parameter $_\n".usage());
  }
  
  # Read the frequencies files
  my $posFrequencies=readPositionFrequencies($$settings{posFrequencies});
  my $ncountPerReadFrequencies=nsPerRead($$settings{nsPerRead});

  my $numNs = introduceNs($posFrequencies,$ncountPerReadFrequencies,$settings);
  logmsg "Introduced $numNs Ns.";

  exit 0;
}

sub readPositionFrequencies{
  my($infile,$settings)=@_;

  my $intspan=Array::IntSpan->new();
  my $bp=0;
  my $currRangeMin=0;
  open(my $fh,$infile) or die "ERROR: could not read $infile: $!";
  while(my $datum=<$fh>){
    chomp($datum);
    my($bp, $count)=split(/\t/, $datum);
    $intspan->set_range($currRangeMin, $count+$currRangeMin, $bp);
    $currRangeMin+=$count+1;
    $bp++;
  }
  close $fh; 
  return $intspan;
}

sub nsPerRead{
  my($infile,$settings)=@_;

  my $intspan=Array::IntSpan->new();

  my $currRangeMin=0;
  open(my $fh, $infile) or die "ERROR: could not read $infile: $!";
  while(my $datum = <$fh>){
    chomp $datum;
    my($bin,$count)=split(/\t/,$datum);
    next if($bin < 1);
    $intspan->set_range($currRangeMin, $count+$currRangeMin, $bin);
    $currRangeMin += $count + 1;
  }
  close $fh;

  return $intspan;
}

sub introduceNs{
  my($posFrequencies,$ncountPerReadFrequencies,$settings)=@_;
  my $totalNs=$$settings{totalNs};
  my $nsThusFar=0;

  # Get maximum number in the intspan so that
  # we have a parameter for rand().
  my @posRanges = sort {$$a[1] <=> $$b[1]} $posFrequencies->get_range_list;
  my $maxPosInt = $posRanges[-1][1];
  my @nsRanges  = sort {$$a[1] <=> $$b[1]} $ncountPerReadFrequencies->get_range_list;
  my $maxNsInt  = $nsRanges[-1][1];

  while(my $id=<>){
    my $seq =<>;
    my $plus=<>;
    my $qual=<>;

    chomp($seq);

    # How many Ns will be added?
    my $numNs = $ncountPerReadFrequencies->lookup(int(rand($maxNsInt)+1));

    my $numTries=0;
    for(my $i=0;$i<$numNs;$i++){
      $numTries++;

      last if($seq=~/^N+$/);

      # Introduce an N somewhere in the read using a random
      # number between 1 and that maximum number.
      my $pos = $posFrequencies->lookup(int(rand($maxPosInt + 1)));

      # Try again if the position has already been modified.
      if(substr($seq, $pos, 1) eq "N"){
        if($numTries >= 100000){
          logmsg "WARNING: tried substituting N $numTries times for a total of $numNs Ns but could not find a free position in $seq (length: ".length($seq).")";
          last;
        }
        $i--;
        next;
      }

      substr($seq, $pos, 1) = "N";
      $nsThusFar++;

      last if($nsThusFar >= $totalNs);
    }

    print "$id$seq\n$plus$qual";

    # If we have finished with the total number of Ns,
    # then just copy the rest of the input file.
    if($nsThusFar >= $totalNs){
      while(<>){
        print $_;
      }
    }
  }

  return $nsThusFar;
}

sub usage{
  "$0: introduce Ns into a read set randomly, according to a profile.
       Ns are introduced in a read set in order, and so it is up to
       the user to randomize the read order.
  Usage: zcat file.fastq.gz | $0 [options] | gzip -c > out.fastq.gz

  Options:
  --posFrequencies  A file with counts of Ns per line. Will be
                    used in a probabalistic fashion to introduce Ns.
  --nsPerRead       A file with counts of how many Ns there are per 
                    read from a real dataset.
  --totalNs         How many Ns total to introduce. Default: 0.

  For both --posFrequencies and --nsPerRead, the file format is
  tab delimited with first column bin, second column count.
  The count is a whole number and not a frequency.
  "
}

