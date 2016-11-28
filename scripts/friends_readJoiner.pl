#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Statistics::Descriptive;
use POSIX qw/floor ceil/;
use File::Basename qw/basename fileparse/;
use List::Util qw/min/;

use FindBin qw/$RealBin/;
use lib "$RealBin/../lib/perl5";
use Friends qw/logmsg mktempdir openFastq/;


# Initialize the components to do some binary operations
# https://medium.com/@keithwhor/nbeam-how-i-wrote-an-ultra-fast-dna-sequence-alignment-algorithm-in-javascript-c199e936da#.ur3qwakpq
my %NT=(
  A  => 2**3,
  C  => 2**0,
  G  => 2**1,
  T  => 2**2,
  GAP=> 2**4,
  #N  => 2**5,
);

# Ambiguity codes for binary operations
my %AMB=(
  U  => $NT{T},
  W  => $NT{A}|$NT{T},   # weak
  S  => $NT{C}|$NT{G},   # strong
  K  => $NT{G}|$NT{T},   # keto
  M  => $NT{A}|$NT{C},   # amino
  B  => $NT{C}|$NT{G}|$NT{T}, # not A
  V  => $NT{A}|$NT{C}|$NT{G}, # not T or U
  D  => $NT{A}|$NT{G}|$NT{T}, # not C
  H  => $NT{A}|$NT{C}|$NT{T}, # not G
  N  => $NT{A}|$NT{C}|$NT{G}|$NT{T},
);
for my $amb(keys %AMB){
  $NT{$amb}=$AMB{$amb};
}

# Trick to get the hash inversion:
# http://perldoc.perl.org/functions/reverse.html
my %REV_NT=reverse(%NT);

exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help prefix=s overlap=i numcpus=i)) or die $!;
  $$settings{numcpus}||=1;
  $$settings{overlap}||=15;
  $$settings{prefix} ||="";

  my($fastq)=@ARGV;
  die usage() if(!$fastq || $$settings{help});
  die "ERROR: I could not find fastq at $fastq" if(!-e $fastq);

  my $reads=joinReads($fastq,$settings);

  # Write out the fwd/rev reads
  if($$settings{prefix}){
    open(FWD,">","$$settings{prefix}.fwd.fastq") or die "ERROR: could not write to $$settings{prefix}.fwd.fastq: $!";
    for(@{$$reads{fwd}}){
      print FWD $_;
    }
    close FWD;

    open(REV,">","$$settings{prefix}.rev.fastq") or die "ERROR: could not write to $$settings{prefix}.rev.fastq: $!";
    for(@{$$reads{rev}}){
      print REV $_;
    }
    close REV;
  }

  return 0;
}

sub joinReads{
  my($fastq,$settings)=@_;

  my(@fwdRead,@revRead);
  my $fh=openFastq($fastq,$settings);
  while(my $fwdRead=<$fh>.<$fh>.<$fh>.<$fh>){
    my $revRead=<$fh>.<$fh>.<$fh>.<$fh>;
    
    # Calculate and print the merged read, but if it doesn't 
    # work, then save the forward and reverse
    if(! calculateOverlap($fwdRead,$revRead,$settings) ){
      push(@fwdRead,$fwdRead);
      push(@revRead,$revRead);
    }
  }

  return {fwd=>\@fwdRead, rev=>\@revRead};
}

sub calculateOverlap{
  my($fwdRead,$revRead,$settings)=@_;
  chomp($fwdRead,$revRead);
  
  # Grab sequences
  my($idF,$sequenceF,$plusF,$qualF)=split("\n",$fwdRead);
  my($idR,$sequenceR,$plusR,$qualR)=split("\n",$revRead);
  my $readLengthF=length($sequenceF);
  my $readLengthR=length($sequenceR);

  # Reverse-complement the reverse read so that things can match up
  # in a straightforward way
  $sequenceR=~tr/ACGTacgt/TGCAtgca/;
  my @sequenceR=reverse(split(//,$sequenceR));
  $qualR=reverse($qualR);

  # Convert the forward and reverse read to binary
  my(@binarySequenceF,@binarySequenceR);
  for(my $i=0;$i<$readLengthF;$i++){
    push(@binarySequenceF,
      $NT{ substr($sequenceF,$i,1) }
    )
  }
  for(my $i=0;$i<$readLengthR;$i++){
    # For the reverse read, it is already in an array.
    push(@binarySequenceR,
      $NT{ $sequenceR[$i] }
    )
  }

  # Where does the matching happen?
  my $index=searchForFirstMatch(\@binarySequenceF,\@binarySequenceR,$$settings{overlap});
  
  # Merge the read
  my $mergedRead=mergeReads($index,\@binarySequenceF,\@binarySequenceR,$qualF,$qualR,$$settings{overlap});

  if($$mergedRead{joined}){
    print "$idF\n$$mergedRead{joined}\n+\n$$mergedRead{joinedQual}\n";
    return 1;
  } else {
    return 0;
  }
}

sub searchForFirstMatch{
  my($binarySequenceF,$binarySequenceR,$overlap)=@_;

  # Search for the first match
  my $searchSpan=min(scalar(@$binarySequenceF),scalar(@$binarySequenceR))-$overlap;
  SEARCH: 
  for(my $i=0;$i<$searchSpan;$i++){
    # Search over position 14 down to position 0, etc
    my $endPos=$i+$overlap-1;
    for(my $j=$endPos;$j>=$i;$j--){
      # Use binary operations to make this match go faster
      if(($$binarySequenceR[$j] & $$binarySequenceF[$j]) == 0){
        $i=$endPos;
        next SEARCH;
      }
    }

    return $i;
  }

  # If no match is found, return -1
  return -1;
}
sub mergeReads{
  my($index,$binarySequenceF,$binarySequenceR,$qualF,$qualR,$overlap)=@_;

  my %return=(
    fwd=>$binarySequenceF,rev=>$binarySequenceR,
  );

  if($index == -1){
    return \%return;
  }

  # The merged read starts with whatever is not going to be joined
  # from the forward read, plus the overlap length that was
  # calculated from searchForFirstMatch().
  my @mergedBinary=@$binarySequenceF[0..$index-1+$overlap];
  my $mergedQual=substr($qualF,0,$index+$overlap);
  my @qualF=split(//,$qualF);
  my @qualR=split(//,$qualR);
  
  # Figure out how far the match goes
  my $lengthF=scalar(@$binarySequenceF);
  my $lengthR=scalar(@$binarySequenceR);
  my $mergeSpan=min($lengthF,$lengthR);
  my $i=$index+$overlap;
  while($i<$mergeSpan){
    last if(($$binarySequenceF[$i] & $$binarySequenceR[$i]) == 0);

    $i++;
  }
  # TODO recalculate quality scores
  push(@mergedBinary, @$binarySequenceF[$index+$overlap..$i-1]);
  $mergedQual.=substr($qualF,$index+$overlap, $i-$index-$overlap);

  # Finally, push the rest of the reverse read onto the merged read
  push(@mergedBinary, splice(@$binarySequenceR, $i));
  #push(@mergedBinary, @$binarySequenceR[$i..
  $mergedQual.=substr($qualR,$i);

  # TODO if the merge includes ALL of the reverse read,
  # such as in the diagram below, then the rest of the 
  # forward read should be included.
  #    F -------------------
  #    R        --------
  #if($i >= $lengthR && $mergeSpan < $lengthF){
  #  push(@mergedBinary,splice(@$binarySequenceF,$i));
  #  $mergedQual.=substr($qualF,$i);
  #  my $fwdRead=join("", map{$REV_NT{$_}} @$binarySequenceF);
  #  my $revRead=join("", map{$REV_NT{$_}} @$binarySequenceR);
  #  my $mergedSequence=join("", map{$REV_NT{$_}} @mergedBinary);
  #  die Dumper [scalar(@mergedBinary),length($mergedQual),$fwdRead,$revRead,$mergedSequence];
  #}

  my $mergedSequence="";
  $mergedSequence.=$REV_NT{$_} for(@mergedBinary);

  #my $fwdRead=join("", map{$REV_NT{$_}} @$binarySequenceF);
  #my $revRead=join("", map{$REV_NT{$_}} @$binarySequenceR);
  #die "\n$fwdRead\n$revRead\n$mergedSequence\n$mergedQual";

  @return{qw(joined start stop joinedQual)}=($mergedSequence,$index,$i-1,$mergedQual);
  return \%return;
}

sub usage{
  local $0=basename $0;
  "$0: join fwd/rev reads when then overlap
  Usage: $0 reads.fastq[.gz] > joined.fastq

  --prefix  ''        If supplied, the unjoined reads will be in 
                      prefix.fwd.fastq and prefix.rev.fastq

  --overlap 15        Minimum number of nucleotides for an
                      overlap
  "
}

