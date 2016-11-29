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
);

# Ambiguity codes for binary operations
my %AMB=(
  U  => $NT{T},
  W  => $NT{A}|$NT{T},        # weak
  S  => $NT{C}|$NT{G},        # strong
  K  => $NT{G}|$NT{T},        # keto
  M  => $NT{A}|$NT{C},        # amino
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
  GetOptions($settings,qw(help prefix=s overlap=i numcpus=i tolerance=i)) or die $!;
  $$settings{numcpus}||=1;
  $$settings{overlap}||=15;
  $$settings{prefix} ||="";
  $$settings{tolerance}||=0;

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
  $sequenceR=reverse($sequenceR);
  $qualR=reverse($qualR);

  # Convert the forward and reverse read to binary
  my(@binarySequenceF,@binarySequenceR);
  for(my $i=0;$i<$readLengthF;$i++){
    push(@binarySequenceF,
      $NT{ substr($sequenceF,$i,1) }
    )
  }
  for(my $i=0;$i<$readLengthR;$i++){
    push(@binarySequenceR,
      $NT{ substr($sequenceR,$i,1) }
    )
  }

  # Where does the matching happen?
  my($indexF,$indexR)=searchForFirstMatch(\@binarySequenceF,\@binarySequenceR,$$settings{overlap});
  
  # Merge the read
  my $mergedRead=mergeReads($indexF,$indexR,\@binarySequenceF,\@binarySequenceR,$qualF,$qualR,$$settings{overlap},$$settings{tolerance});

  if($$mergedRead{joined}){
    print "$idF\n$$mergedRead{joined}\n+\n$$mergedRead{joinedQual}\n";
    return 1;
  } else {
    return 0;
  }
}

# Return the coordinate of EACH read instead of
# assuming that they line up exactly.
sub searchForFirstMatch{
  my($binarySequenceF,$binarySequenceR,$overlap)=@_;

  my $lengthF=scalar(@$binarySequenceF)-$overlap;
  my $lengthR=scalar(@$binarySequenceR)-$overlap;

  # Search for the match, where $f and $r are coordnates
  # on the fwd and rev reads.
  SEARCH_F:
  for(my $f=0;$f<$lengthF;$f++){
    SEARCH_R:
    for(my $r=0;$r<$lengthR;$r++){
      # Go backwards from positions $f and $r + $overlap
      # to see if there is an overlap.
      my $i;
      for($i=$overlap;$i>=0;$i--){
        #last if($i-$r > $lengthR || $i-$f > $lengthF);
        # Test if there is a match using binary operator.
        # If there is no match, then jump ahead by $overlap
        # length and try again.
        if(($$binarySequenceR[$r+$i] & $$binarySequenceF[$f+$i]) == 0){
          #$f+=$overlap;
          $r+=$overlap-1; # advance by overlap - 1 because $r will advance by 1 in the next iteration
          next SEARCH_R;
        }
      }

      last if($i >= 0); # make sure we made it through the loop ok
      #my $stopF=$overlap+$f;
      #my $stopR=$overlap+$r;
      return ($f,$r);
    }
  }

  return (-1,-1);
}

sub mergeReads{
  my($indexF,$indexR,$binarySequenceF,$binarySequenceR,$qualF,$qualR,$overlap,$tolerance)=@_;

  my %return=(
    fwd=>$binarySequenceF,rev=>$binarySequenceR,
  );

  if($indexF == -1){
    return \%return;
  }

  # TODO figure out if F or R starts first
  #
  # There are four scenarios if there are $indexF && $indexR:
  #   F -----------
  #          -----------R
  #
  #   F      -----------
  #   R ----------
  #   
  #   F -----------------
  #   R        ---
  #
  #   F        ---
  #   R -----------------
  #
  #   and a special fifth
  #   F -----------------
  #   R -----------------

  # The merged read starts with whatever is not going to be joined
  # from the forward read, plus the overlap length that was
  # calculated from searchForFirstMatch().

  my $lengthF=scalar(@$binarySequenceF);
  my $lengthR=scalar(@$binarySequenceR);
  # Start positions for matching can go ahead by $overlap
  # to reduce redundancy since this was already done in the
  # previous subroutine.
  my $f=$indexF+$overlap-1;
  my $r=$indexR+$overlap-1;
  while($f < $lengthF && $r < $lengthR){
    # Use a binary operator to see if there is a match
    if(($$binarySequenceF[$f] & $$binarySequenceR[$r]) == 0){
      $f--;
      $r--;
      last;
    }

    $f++;
    $r++;
  }

  # TODO
  # Right now if the merge is in the middle of the read, the first
  # part of the unmerged read will come from F and the unmerged
  # at the tail will come from R.  However, the highest quality
  # should be used.

  # The range of matching is $indexF to $f and $indexR to $r.
  my @mergedBinarySequence;
  if($indexF <= $indexR){
    push(@mergedBinarySequence,@$binarySequenceF[0..$indexF-1]);
  } elsif($indexR < $indexF){
    push(@mergedBinarySequence,@$binarySequenceR[0..$indexR-1]);
  }
  
  # Accept the merged sequence
  push(@mergedBinarySequence,@$binarySequenceF[$indexF..$f]);

  # Finish off the merge
  if($r<=$f){
    push(@mergedBinarySequence,@$binarySequenceR[$r+1..$lengthR-1]);
  }elsif($f < $r){
    push(@mergedBinarySequence,@$binarySequenceF[$f+1..$lengthF-1]);
  }

  my $marginF=$lengthF-$f;
  my $marginR=$lengthR-$r;
  print join("\n",
    "$indexF-$f($marginF)  $indexR-$r($marginR)",
    " " x $marginR . binaryToDna(\@mergedBinarySequence),
    " " x $marginF . binaryToDna($binarySequenceF),
    " " x $marginR . binaryToDna($binarySequenceR)
  )."\n";
  die "TODO more validation";

  my $mergedSequence=binaryToDna(\@mergedBinarySequence);

  #@return{qw(joined start stop joinedQual)}=($mergedSequence,$index,$i-1,$mergedQual);

  return \%return;
}

sub binaryToDna{
  my($binarySequence,$settings)=@_;
  my $dna="";
  for(@$binarySequence){
    $dna.=$REV_NT{$_};
  }
  return $dna;
}

sub usage{
  local $0=basename $0;
  "$0: join fwd/rev reads when then overlap
  Usage: $0 reads.fastq[.gz] > joined.fastq

  --prefix    ''      If supplied, the unjoined reads will be in 
                      prefix.fwd.fastq and prefix.rev.fastq

  --overlap   15      Minimum number of nucleotides for an
                      overlap
  "
}

