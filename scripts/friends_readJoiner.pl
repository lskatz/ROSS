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

  # This funtion prints any reads that can be joined
  # and then returns unjoined reads
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
  }else{
    logmsg "Forward and reverse reads were not requested";
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
    if(! printMergedRead($fwdRead,$revRead,$settings) ){
      push(@fwdRead,$fwdRead);
      push(@revRead,$revRead);
    }
  }

  return {fwd=>\@fwdRead, rev=>\@revRead};
}

sub printMergedRead{
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

  if($$mergedRead{joinedSequence}){
    print "$idF-joined\n$$mergedRead{joinedSequence}\n+\n$$mergedRead{joinedQual}\n";
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
      for($i=$overlap-1;$i>=0;$i--){
        #last if($i-$r > $lengthR || $i-$f > $lengthF);
        # Test if there is a match using binary operator.
        # If there is no match, then jump ahead by $overlap
        # length and try again.
        if(($$binarySequenceR[$r+$i] & $$binarySequenceF[$f+$i]) == 0){
          # advance by what we've looked at so far
          if($i-2 > 1){
            $r+=$i-2;
          }
          next SEARCH_R;
        }
      }

      # $i is -1 if we made it through the loop ok
      last if($i != -1);
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
    #fwd=>binaryToDna($binarySequenceF),rev=>binaryToDna($binarySequenceR),joinedSequence=>"",
    joinedSequence=>"",
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
  #
  #   Do not accept the second case, where F and R are opposite
  #   of what you'd expect.

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
  
  # If the rev read is inside of the fwd read,
  # and the entire read is the match to the fwd,
  # return the more complete read.
  if(
       $indexR == 0       # the rev-read match starts with the beginning of the read
    && $lengthR == $r     # the rev-read match goes to the end of the read
  ){
    @return{qw(joinedSequence joinedQual)}=(binaryToDna($binarySequenceF),$qualF);
    return \%return;
  }
  # If the fwd read is inside of the rev read,
  # and the entire read is the match to the rev,
  # return the more complete read.
  elsif(
       $indexF == 0       # the fwd-read match starts with the beginning of the read
    && $lengthF == $f     # the fwd-read match goes to the end of the read
  ){
    @return{qw(joinedSequence joinedQual)}=(binaryToDna($binarySequenceR),$qualR);
    return \%return;
  }
  # If the match doesn't extend to the end of the fwd read,
  # then this is also not quite a match. Return false.
  # Same with not not matching on the start of the reverse read.
  elsif($lengthF != $f || $indexR != 0){
    return \%return;
  }

  # Subtract one more time to make the range inclusive
  $f--;
  $r--;

  # TODO
  # Right now if the merge is in the middle of the read, the first
  # part of the unmerged read will come from F and the unmerged
  # at the tail will come from R.  However, the highest quality
  # should be used.
  my @phredF=split(//,$qualF);
  my @phredR=split(//,$qualR); # QualR has already been reversed to match the order of nts
  for(@phredF,@phredR){
    $_=ord($_)-33;
  }

  # The range of matching is $indexF to $f and $indexR to $r.
  my (@mergedBinarySequence,@mergedPhred);
  if($indexR <= $indexF){ # lower indexR indicates that the fwd read is first
    push(@mergedBinarySequence,@$binarySequenceF[0..$indexF-1]);
    push(@mergedPhred,@phredF[0..$indexF-1]);
  } elsif($indexR < $indexF){
    push(@mergedBinarySequence,@$binarySequenceR[0..$indexR-1]);
    push(@mergedPhred,@phredR[0..$indexR-1]);
  }
  
  # Accept the merged sequence
  push(@mergedBinarySequence,@$binarySequenceF[$indexF..$f]);
  push(@mergedPhred,@phredR[$indexF..$f]); # for now, use fwd quality scores

  # Finish off the merge.
  # If the $r is less than $f, then we have this situation
  # F -------------
  # R     --------------
  if($r<=$f){
    # Add the rest of the rev read
    push(@mergedBinarySequence,@$binarySequenceR[$r+1..$lengthR-1]);
    push(@mergedPhred,@phredR[$r+1..$lengthR-1]);
  }elsif($f < $r){
    # Add the rest of the fwd read
    push(@mergedBinarySequence,@$binarySequenceF[$f+1..$lengthF-1]);
    push(@mergedPhred,@phredF[$f+1..$lengthF-1]);
  }

  # DEBUG
  #my $padding=$indexF-$indexR;
  #my $str="";
  #$str.=substr($_,-1) for(0..40);
  #$str.="\n      ";
  #$str.=substr($_,-1) for(0..40);
  #die "\n$indexF-$f $indexR-$r\n$str\n".binaryToDna($binarySequenceF)."\n". " " x $padding . binaryToDna($binarySequenceR)."\n";
  
  #die binaryToDna(\@mergedBinarySequence)."\n".binaryToDna(\@mergedPhred) if(@mergedBinarySequence != @mergedPhred);
  
  @return{qw(joinedSequence joinedQual)}=(binaryToDna(\@mergedBinarySequence),phredToQual(\@mergedPhred));

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

sub phredToQual{
  my($phred,$settings)=@_;
  my $qual;
  for(@$phred){
    $_//=0;
    $qual.=chr($_ + 33);
  }
  return $qual;
}

sub mergePhred{
  my($phred1,$phred2)=@_;
  my @mergedPhred;
  my $num=@$phred1;
  for(my $i=0;$i<$num;$i++){
    my $p1=10**((-$$phred1[$i])/10);
    my $p2=10**((-$$phred2[$i])/10);
    push(@mergedPhred,
      -10 * log($p1*$p2)/log(10)
    );
  }
  return \@mergedPhred;
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

