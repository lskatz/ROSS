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

exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help)) or die $!;

  my($fastq)=@ARGV;
  die usage() if(!$fastq || $$settings{help});
  die "ERROR: I could not find fastq at $fastq" if(!-e $fastq);

  my $reads=joinReads($fastq,$settings);

  return 0;
}

sub joinReads{
  my($fastq,$settings)=@_;

  my(@joinedRead,@fwdRead,@revRead);
  my $fh=openFastq($fastq,$settings);
  while(my $fwdRead=<$fh>.<$fh>.<$fh>.<$fh>){
    my $revRead=<$fh>.<$fh>.<$fh>.<$fh>;
    
    calculateOverlap($fwdRead,$revRead,$settings);
  }

  return {joined=>\@joinedRead, fwd=>\@fwdRead, rev=>\@revRead};
}

sub calculateOverlap{
  my($fwdRead,$revRead,$settings)=@_;

}
