#!/usr/bin/env perl

# Friends.pm: library for the Random Operations on Sequences Suite (ROSS)
# Author: Lee Katz <lkatz@cdc.gov>
# Some code taken from CG-Pipeline
# TODO openSam()

package Friends;
require 5.12.0;

use strict;
use warnings;

use File::Basename qw/basename fileparse dirname/;
use File::Temp ('tempdir');
use Data::Dumper;
use IO::Uncompress::Gunzip;

use threads;
use threads::shared;
use Thread::Queue;

use Exporter qw/import/;
our @EXPORT_OK = qw(
             logmsg mktempdir median stdev fullPathToExec isFasta readMfa alnum openFastq
             @fastqExt @fastaExt @bamExt @vcfExt @richseqExt @sffExt @samExt
           );

our @fastqExt=qw(.fastq.gz .fastq .fq .fq.gz);
our @fastaExt=qw(.fasta .fna .faa .mfa .fas .fa);
our @bamExt=qw(.sorted.bam .bam);
our @vcfExt=qw(.vcf.gz .vcf);
our @richseqExt=qw(.gbk .gbf .gb .embl);
our @sffExt=qw(.sff);
our @samExt=qw(.sam .bam);

our $fhStick :shared; # Helps us open only one file at a time


# TODO if 'die' is imported by a script, redefine
# sig die in that script as this function.
$SIG{'__DIE__'} = sub { 
  my $e = $_[0]||""; 
  $e =~ s/(at [^\s]+? line \d+\.$)/\nStopped $1/; 
  my $caller=(caller(1))[3] || (caller(0))[3];
  die("$0: $caller: $e"); 
};

# Centralized logmsg
#sub logmsg {print STDERR "$0: ".(caller(1))[3].": @_\n";}
sub logmsg {
  local $0=basename $0;
  my $parentSub=(caller(1))[3] || (caller(0))[3];
  $parentSub=~s/^main:://;

  # Find the thread ID and stringify it
  my $tid=threads->tid;
  $tid=($tid) ? "(TID$tid)" : "";

  my $msg="$0: $parentSub$tid: @_\n";

  print STDERR $msg;
}



sub mktempdir(;$) {
	my ($settings) = @_;
	my $tempdir_path = File::Spec->join(File::Spec->tmpdir(), (split("::",(caller(1))[3]))[1].".$$.XXXXX");
	my $tempdir = tempdir($tempdir_path, CLEANUP => !($$settings{keep}));
	return $tempdir;
}


# TODO: optimize this
# TODO: setting: unbiased vs. biased estimator
# Can use Statistics::Descriptive.
sub stdev {
	my @list = @_;
	return undef if @list < 1;
	return undef if @list < 1;
	my $mean = sum(@list) / @list;
	return sqrt(sum(map { ($_-$mean)**2 } @list) / @list);
}

# TODO: setting: averaging allowed
# Can use Statistics::Descriptive.
sub median {
	my @list = @_;
	return undef if @list < 1;
#	$list = $$list[0] if @$list == 1 and ref($$list[0]) == 'ARRAY';

	my @slist = sort {$a <=> $b} @list;
	my $count = scalar @list;
	if ($count % 2) {
		return $slist[int($count/2)];
	} else {
		return ($slist[$count/2] + $slist[$count/2 - 1]) / 2;
	}
}

# Return a value in the list such that $percentile percent of values in the list
# are below that value.
# TODO: setting: accept a user-defined sorting function and return index not value
sub percentile($$) {
	my ($list, $percentile) = @_;
	die("argument percentile undefined") if not defined $percentile;
	die("illegal value of argument percentile") if $percentile<0 or $percentile>100;
	return undef if not defined $list or @$list < 1;

	my @slist = sort {$a <=> $b} @$list;
	my $count = scalar @$list;
	return $slist[min(int($count*$percentile/100), $#slist)];
}

# If argument is an executable in the current path, returns the full path to it, otherwise returns undef
# arguments: warn_on_error
sub fullPathToExec($;$) {
	my ($executable,$settings) = @_;
	my $fullpath;
	for ("", split(/:/, $ENV{PATH})) {
    my $path=$_."/".$executable;
		if (-x $path && -f $path) { $fullpath = File::Spec->rel2abs($path); last; }
	}
  if(! -x $fullpath){
	  my $errStr="Error finding full path to executable ($executable)";
    logmsg $errStr if($$settings{warn_on_error});
    die $errStr if(!$$settings{warn_on_error});
  }
	return $fullpath;
}

# TODO: this needs more work
sub isFasta($) {
	my ($file) = @_;
	my $char;
	open(FH, "< $file") or die("Could not open file $file for reading: ".$!);
	read(FH, $char, 1); close FH;
	if ($char eq ">") { return 1; } else { return 0; }
}

# TODO: error checking
# Usage: $seqs = readMfa($seqfile); foreach $seqname (keys %$seqs) { $seq = $$seqs{$seqname}; }
# Usage: (values %{readMfa($seqfile)})[0]
# TODO: return second array with ordering
sub readMfa($;$) {
	my ($mfa_file, $settings) = @_;
	open(FH, '<', $mfa_file) or die("Could not open file $mfa_file for reading: ".$!);
	my %seqs;
	my ($cur_seq_hdr, $seq);
	while (<FH>) {
		if (/^\>\s*(.+)/) {
			$seqs{$cur_seq_hdr} = $seq if $cur_seq_hdr;
			undef $seq;
			$cur_seq_hdr = ($$settings{first_word_only} ? (split /\s/, $1)[0] : $1);
		} else {
			chomp; 
      s/\s//g unless($$settings{keep_whitespace});
      $seq .= uc $_;
		}
	}
	close FH;
	$seqs{$cur_seq_hdr} = $seq if $cur_seq_hdr;
	die("Error: No sequences found in $mfa_file") unless %seqs;
	return \%seqs;
}

sub printSeqsToFile($$;$) {
	my ($seqs, $outfile, $settings) = @_;
	open(OUT, '>', $outfile) or die "Unable to open file $outfile for writing: $!";
	my @seqnames = keys %$seqs;
	@seqnames = sort alnum @seqnames
		if $$settings{order_seqs_by_name};
	@seqnames = sort {length($$seqs{$a}) <=> length($$seqs{$b})} @seqnames
		if $$settings{order_seqs_by_length};
	foreach my $seqname (@seqnames) {
		my $seq = $$seqs{$seqname};
		$seq =~ s/(.{80})/$1\n/g;
		$seq .= "\n" unless $seq =~ /\n$/;
		print OUT ">$seqname\n$seq";
	}
	close OUT;
	return $outfile;
}

# Opens a fastq file in a thread-safe way.
# Returns a file handle.
sub openFastq{
  my($fastq,$settings)=@_;

  my $fh;

  my($name,$dir,$ext)=fileparse($fastq,@fastqExt);

  if(!grep(/$ext/,@fastqExt)){
    die "ERROR: could not read $fastq as a fastq file";
  }

  # Open the file in different ways, depending on if it
  # is gzipped or if the user has gzip installed.
  lock($fhStick);
  if($ext =~/\.gz$/){
    # use binary gzip if we can... why not take advantage
    # of the compiled binary's speedup?
    if(-e "/usr/bin/gzip"){
      open($fh,"gzip -cd $fastq | ") or die "ERROR: could not open $fastq for reading!: $!";
    }else{
      $fh=new IO::Uncompress::Gunzip($fastq) or die "ERROR: could not read $fastq: $!";
    }
  } else {
    open($fh,"<",$fastq) or die "ERROR: could not open $fastq for reading!: $!";
  }

  return $fh;
}

1;
