#!/usr/bin/env perl

# Fix fastq files as best I can

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

# Some globals to make this thing go faster
my @validNt = qw(A T G C N);
my $validNt = join("", @validNt);
my $validRe = qr/[$validNt]/;

my @validQScore;
for my $int(0..60){ # usually, 0 to 40 though
  push(@validQScore, chr($int + 33));
}
my $validScoreStr=join("",@validQScore);
my $validScoreRe = qr/[$validScoreStr]/;
my $worstScore = $validQScore[0];

local $0=basename $0;

exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help)) or die $!;

  die usage() if($$settings{help});
  
  for my $f(@ARGV){
    printFixedFastq($f,$settings);
  }

  return 0;
}

sub printFixedFastq{
  my($file, $settings)=@_;

  my $fh = openFastq($file);
  while(my $id = <$fh>){
    my $seq  = <$fh>;
    my $plus = <$fh>;
    my $qual = <$fh>;
    # Don't accept incomplete entries
    last if(!$id || !$seq || !$plus || !$qual);
    chomp($id, $seq, $plus, $qual);

    last if(length($seq) < 1 || length($qual) < 1);

    $id  = fixId  ($id,  $settings);
    $seq = fixSeq ($seq, $settings);
    $plus= fixPlus($plus, $settings);
    $qual= fixQual($qual, length($seq), $settings);

    print "$id\n$seq\n$plus\n$qual\n";
  }
  close $fh;
}

sub fixId{
  my($id, $settings)=@_;
  $id=~s/^@?/@/; # ensure there is a single @ at the beginning
  return $id;
}
sub fixSeq{
  my($seq, $settings)=@_;

  my $fixed="";
  for my $nt(split(//, $seq)){
    if($nt !~ $validRe){
      $nt="N";
    }
    $fixed.=$nt;
  }
  return $fixed;
}
sub fixPlus{
  my($plus, $settings)=@_;
  $plus=~s/^\+*/+/; # ensure there is a single + at the beginning
  return $plus;
}
sub fixQual{
  my($qual, $length, $settings)=@_;
  
  my $fixed="";
  my $fixedLength=0;
  for my $q(split(//, $qual)){
    if($q !~ $validScoreRe){
      $q = $worstScore;
    }
    $fixed.=$q;
    $fixedLength++;
    # ensure that qual isn't longer than seq
    last if($fixedLength >= $length);
  }
  for(my $i=$fixedLength; $i<=$length; $i++){
    $fixed.=$worstScore;
  }

  return $fixed;
}

sub usage{
  "$0: Fix fastq files
  Usage: $0 file.fastq[.gz] > out.fastq
  "
}
