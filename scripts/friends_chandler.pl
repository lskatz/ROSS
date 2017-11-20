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
use Friends qw/logmsg mktempdir openFastq/;
use Bio::Kmer;

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
  GetOptions($settings,qw(help sortby|sort-by=s kmerlength|kmer=i gt=i kmerCounter=s delta=i tempdir=s numcpus=i)) or die $!;
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

  # Auto-choose if the kmer counter is not given
  if(!$$settings{kmerCounter}){
    if(which("jellyfish")){
      $$settings{kmerCounter}="jellyfish";
    } else {
      $$settings{kmerCounter}="perl";
    }
  }

  # Kmer counting object
  my $kmer=Kmer->new($fastq,{
      numcpus=>$$settings{numcpus}, 
      kmerlength=>$$settings{kmerlength}, 
      gt=>$$settings{gt},
      kmercounter=>$$settings{kmerCounter},
  });
  # Run the counting
  my $kmercount=$kmer->count();

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

