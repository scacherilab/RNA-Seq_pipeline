#!/usr/bin/perl -w
# $Id: merge.plx,v 1.6 2003/06/12 17:59:01 pchines Exp $

use strict;
use Getopt::Long;
use Pod::Usage;

=head1 NAME

merge.plx - merge files together based on common key column

=head1 SYNOPSIS

merge.plx [options] file1name:col file2name:col

=head1 DESCRIPTION

This program will merge together delimited text files; both files must use
the same delimiter, in this version.  In order to merge the
files there must be one column in each file where the values correspond to
one another.  The output will contain a line with the concatenated values
from both files where the files all have the same key value in common.  Thus,
merge performs an intersection between the files, and outputs a single line
for each intersection, containing all of the data from all of the files.

Typically each key value will appear only once in each file.  In this version
of the merge program, however, multiple values may be handled in one of three
ways.  See options B<-a>, B<-f> and B<-l> below.  If none of these options
are specified and multiple values for the same key are encountered, the
program will terminate with an error.

=cut

use vars qw(@FILES %KEY_COLUMN $DELIM $CASE_INSENSITIVE
            $USE_FIRST $USE_LAST $USE_ALL $OUTER_JOIN $NEWLINE);

process_commandline();
my $rh_cumulative;
foreach my $file (@FILES) {
    my $rh = read_file($file, $KEY_COLUMN{$file});
    $rh_cumulative = intersect_and_concatenate($rh_cumulative, $rh);
}
foreach my $key (sort keys %$rh_cumulative) {
    if (ref $rh_cumulative->{$key}) {
        print "$_\n" foreach (@{ $rh_cumulative->{$key} });
    }
    else {
        print $rh_cumulative->{$key}, "\n";
    }
}

## End MAIN

sub intersect_and_concatenate {
    my ($rh_cum, $rh_add) = @_;
    my $nfields;
    foreach my $key (keys %$rh_add) {
        my @fields = split($DELIM, $rh_add->{$key});
        $nfields = scalar(@fields);
        last;
    }
    my %new;
    if ($rh_cum) {
        foreach my $key (keys %$rh_cum) {
            if (exists $rh_add->{$key}) {
                $new{$key} = concatenate($rh_cum->{$key}, $rh_add->{$key});
            }
            elsif ($OUTER_JOIN) {
                $new{$key} = concatenate($rh_cum->{$key},
                    $DELIM x ($nfields-1) );
            }
        }
    }
    else {
        %new = %$rh_add;
    }
    return \%new;
}

sub concatenate {
    my ($val1, $val2) = @_;
    my $newval = [];
    if (ref($val1) && ref($val2)) {
        foreach my $v1 (@$val1) {
            foreach my $v2 (@$val2) {
                push @$newval, join($DELIM, $v1, $v2);
            }
        }
    }
    elsif (ref $val1) {
        foreach my $v1 (@$val1) {
            push @$newval, join($DELIM, $v1, $val2);
        }
    }
    elsif (ref $val2) {
        foreach my $v2 (@$val2) {
            push @$newval, join($DELIM, $val1, $v2);
        }
    }
    else {
        $newval = join($DELIM, $val1, $val2);
    }
    return $newval;
}

sub read_file {
    my ($file, $col) = @_;
    my %data;
    open(FILE, "<$file") || die "Can't open $file, $!\n";
    local $/ = guess_newline(*FILE);
    $NEWLINE ||= $/;
    while (<FILE>) {
        chomp;
        my $key = (split($DELIM))[$col-1];
        if ($CASE_INSENSITIVE) {
            $key = uc($key);
        }
        if ($data{$key}) {
            if ($USE_FIRST) {
                # Leave current data alone
            }
            elsif ($USE_LAST) {
                $data{$key} = $_;
            }
            elsif ($USE_ALL) {
                if (ref $data{$key}) {
                    push @{ $data{$key} }, $_;
                }
                else {
                    $data{$key} = [ $data{$key}, $_ ];
                }
            }
            else {
                die "Encountered key value '$key' multiple times at line "
                    .__LINE__." of file '$file'.\n"
                    ."Use option -a, -f or -l to determine how to handle.\n";
            }
        }
        else {
            $data{$key} = $_;
        }
    }
    close FILE;
    return \%data;
}

sub guess_newline {
    my ($fh) = @_;
    my $newline;
    my $text;
    if (read($fh, $text, 4096)) {
        if ($text =~ m/(\015?\012|\015)/) {
            $newline = $1;
        }
    }
    seek $fh, 0, 0;
    return $newline;
}

sub process_commandline {
    my $help;
    $DELIM = "\t";
    $USE_ALL = $USE_FIRST = $USE_LAST = 0;
    $OUTER_JOIN = 0;
    GetOptions(
        'help'  => \$help,
        'delimiter=s' => \$DELIM,
        'i'     => \$CASE_INSENSITIVE,
        'all'   => \$USE_ALL,
        'first' => \$USE_FIRST,
        'last'  => \$USE_LAST,
        'outer' => \$OUTER_JOIN,
        ) || pod2usage(1);
    pod2usage(-verbose => 2) if $help;
    pod2usage(1) if $USE_ALL + $USE_FIRST + $USE_LAST > 1;

    @FILES = @ARGV;
    if (@FILES < 2) {
        pod2usage("Need at least 2 files to merge\n");
    }
    for my $file (@FILES) {
        pod2usage("Must supply column numbers of key columns") unless $file =~ /:\d+$/;
        pod2usage("No other ':' allowed in file names") if $file =~ /:.*:\d+$/;
    }
    %KEY_COLUMN = split(/:/, join(":", @FILES));
    @FILES = map { s/:\d+$//;$_ } @FILES;
}

=head1 OPTIONS

=over 4

=item B<-a|--all>

Where there are multiple lines that share the same key value, keep all of
them in the output.  This results in a Cartesian product, i.e. if there are
two rows in the first file and three in the second file, there will be six
rows in the output.  Mutually exclusive with options B<-f> and B<-l>.

=item B<-f|--first>

Where there are multiple lines that share the same key value, keep only the
first of these encountered in the input.  Mutually exclusive with options
B<-a> and B<-l>.

=item B<-l|--last>

Where there are multiple lines that share the same key value, keep only the
last of these encountered in the input.  Mutually exclusive with options
B<-a> and B<-f>.

=item B<-D|--delimiter>

Set the column delimiter (defaults to tab).

=item B<-o|--outer>

Perform an outer join, i.e. include all records from earlier files, even if
a matching record does not appear in later files.

=item B<-h|--help>

Display this help information.

=item B<-i>

Match key columns without regard to case.

=back

=head1 AUTHOR

Peter Chines - pchines@nhgri.nih.gov

=cut

