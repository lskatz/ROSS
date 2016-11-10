#!/usr/bin/env perl

# ROSS.pm: library for the Random Operations on Sequences Suite
# Author: Lee Katz <lkatz@cdc.gov>
# Some code taken from CG-Pipeline

package ROSS;
require 5.12.0;

use strict;
use warnings;

use File::Basename qw/basename fileparse dirname/;
use File::Temp ('tempdir');
use Data::Dumper;

use threads;
use threads::shared;
use Thread::Queue;

use Exporter;
our @EXPORT_OK = qw(
             logmsg mktempdir median stdev fullPathToExec isFasta readMfa alnum
           );


# TODO if 'die' is imported by a script, redefine
# sig die in that script as this function.
local $SIG{'__DIE__'} = sub { my $e = $_[0]; $e =~ s/(at [^\s]+? line \d+\.$)/\nStopped $1/; die("$0: ".(caller(1))[3].": ".$e); };

sub logmsg {my $FH = $AKUtils::LOG || *STDOUT; print $FH "$0: ".(caller(1))[3].": @_\n";}

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

# Function to help sort alphanumerically
# Usage: @sorted=sort alum @unsorted;
	foreach my $seqname (sort alnum keys %$seqs) {
sub alnum {
	my ($i);
	my ($len1, $len2) = (length($a), length($b));
	for ($i = 0; ($i < $len1) && ($i < $len2); ++$i) {
		my $c1 = substr($a, $i, 1);
		my $c2 = substr($b, $i, 1);
		($c1 =~ /^\d/o) || ($c2 =~ /^\d/o) || ($c1 ne $c2) and last;
	}
	my $a_r = ($i < $len1) ? substr($a, $i) : "";
	my $b_r = ($i < $len2) ? substr($b, $i) : "";
	my ($a_n, $a_s) = ($a_r =~ /^(\d+)(.*)$/o);
	my ($b_n, $b_s) = ($b_r =~ /^(\d+)(.*)$/o);
	return (defined($a_n) && defined($b_n)) ?
	(($a_n <=> $b_n) || ($a_s cmp $b_s)) : ($a cmp $b);
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

1;
