#!/usr/bin/perl

our $VERSION = 0.03;

# SCRIPT REMOVES SEQUENCES SHORTER THAN 200 bases
# SCRIPT REMOVES ADAPTERS


# NOTE: I rewrote the following line in adapters.txt:
# Aala.03.053215  3076    2727..2755,2776..2799
# to:
# Aala.03.053215  3076    2727..2799
# 
# Also rewrote
# Aala.03.257735  414     1..25,379..414  adaptor:NGB01030.1
# to
# Aala.03.257735  414     1..25  adaptor:NGB01030.1
# 
# and take care of the end of this in the script
# see occurrences of $MULTI_ADAPTOR_NOT_AFFECTING_GENE

use strict;
use warnings;
use autodie;
use FileHandle;
use JFR::Fasta;
use Data::Dumper;

our $SCF_W_GENE_SPANNING_INTERNAL_ADAPTER = 'Aala.03.045289';
our $MULTI_ADAPTOR_NOT_AFFECTING_GENE     = 'Aala.03.257735';

our $FA = '/bwdata1/jfryan/15-ALATINA/00-DATA/Aala.03';
our $CONTAM = 'adaptors.txt';
our $GFF3 = '/bwdata1/jfryan/15-ALATINA/02-MASURCA.02/04-AUGUSTUS/Alatina_alata.masurca.v3.gff3';
our $OUTFA = 'Aala.06';
our $OUTGFF = 'Aala.06.gff3';

# note the following gens have adaptor overlap:
# Aala.03.005638 Aala.03.005638 Aala.03.005638 Aala.03.005638
our %BEG_ADAPTOR_W_GENE_OVERLAP = ('Aala.06.005638' => 1);
MAIN: {
    my $rh_contam = get_contam_data($CONTAM);  
    my $rh_gff3    = get_gff3($GFF3);
    my $fp = JFR::Fasta->new($FA);
    my $out = FileHandle->new($OUTFA, 'w');
    my $gffout = FileHandle->new($OUTGFF, 'w');
    while (my $rec = $fp->get_record()) {
        my $id = JFR::Fasta->get_def_w_o_gt($rec->{'def'});
        $id =~ m/^Aala.03\.(\d+)/ or die "unexpceted ID: $id";
        my $newid = "Aala.06.$1";
        next unless (length($rec->{'seq'}) > 200);

        die "unexpected" if ($rh_contam->{$id} && $rec->{'seq'} =~ m/^[Nn]/);
        $rec->{'seq'} =~ s/[Nn]+\s*$//;

        if ($rec->{'seq'} =~ m/^([Nn]+)/) {
            my $offset = length($1);
            $rec->{'seq'} =~ s/^[Nn]+//;
            print $out ">$newid\n$rec->{'seq'}\n";
            print_beg_adjusted_gff($rh_gff3->{$id},$gffout,$newid,$offset);
        } elsif ($rh_contam->{$id}) {
            my $pos = beg_mid_or_end($rh_contam->{$id});
            if ($pos eq 'beg') {
                my $offset = 0;
                if ($id eq $MULTI_ADAPTOR_NOT_AFFECTING_GENE) {
                    my $end = 414 - 379 + 1; # deleted from adaptors.txt
                    $rec->{'seq'} =~ m/^.{$rh_contam->{$id}->[2]}(.*).{$end}$/;
                    print $out ">$newid\n$1\n";
                } else {
                    $rec->{'seq'} =~ m/^.{$rh_contam->{$id}->[2]}(.*)/;
                    my $seq = $1;
                    if ($seq =~ m/^([Nn]+)/) {
                        $offset += length($1);
                        $seq =~ s/^[Nn]+//;
                    }
                    $seq =~ s/[Nn]+\s*$//;
                    print $out ">$newid\n$seq\n";
                }
                $offset += $rh_contam->{$id}->[2] - $rh_contam->{$id}->[1] +1;
                print_beg_adjusted_gff($rh_gff3->{$id},$gffout,$newid,$offset);
            } elsif ($pos eq 'mid') {
#ID => [TIG_LEN, START, END]
# Aala.03.045289	1657	1188..1207
                my $num5 = $rh_contam->{$id}->[1] - 1;
                my $num3 = $rh_contam->{$id}->[0] - $rh_contam->{$id}->[2];
                my $numa = $rh_contam->{$id}->[2] - $rh_contam->{$id}->[1] + 1;
                $rec->{'seq'} =~ m/^(\S{$num5})\S{$numa}(\S{$num3})$/ or die "unexpected";
                my $seq5 = $1;
                my $seq3 = $2;
                $seq5 =~ s/([Nn]+)\s*$//;
                my $trim5_back = length($1);
                $seq3 =~ s/^([Nn]+)//;
                my $trim3_front = length($1);
                print $out ">${newid}a\n$seq5\n";
                print $out ">${newid}b\n$seq3\n";
                foreach my $ra_r (@{$rh_gff3->{$id}}) { 
                    if ($ra_r->[3] < $num5) {
                        $ra_r->[0] = "${newid}a";
                        my $gffline = join "\t", @{$ra_r};
                        print $gffout "$gffline\n";
                    } else {
                        my $offset = $rh_contam->{$id}->[2];
                        $ra_r->[0] = "${newid}b";
                        $ra_r->[3] -= $offset;
                        $ra_r->[4] -= $offset;
                        my $gffline = join "\t", @{$ra_r};
                        print $gffout "$gffline\n";
                    }
                }
            } elsif ($pos eq 'end') {
                my $nchar = $rh_contam->{$id}->[2] - $rh_contam->{$id}->[1] + 1;
                $rec->{'seq'} =~ m/^(.*)\S{$nchar}/;
                print $out ">$newid\n$1\n";
                foreach my $ra_r (@{$rh_gff3->{$id}}) { 
                    $ra_r->[0] = $newid;
                    my $gffline = join "\t", @{$ra_r};
                    print $gffout "$gffline\n";
                }
            } else {
                die "unexpected - not beg or end";
            }
        } else {
            print $out ">$newid\n$rec->{'seq'}\n";
            foreach my $ra_r (@{$rh_gff3->{$id}}) {
                $ra_r->[0] = $newid;
                my $gffline = join "\t", @{$ra_r};
                print $gffout "$gffline\n";
            }
        }
    }
}

sub print_beg_adjusted_gff {
    my $ra_gff3 = shift;
    my $gffout  = shift;
    my $newid   = shift;
    my $offset  = shift;
    foreach my $ra_r (@{$ra_gff3}) {
        $ra_r->[0] = $newid;
        $ra_r->[3] -= $offset;
        $ra_r->[4] -= $offset;
        if ($BEG_ADAPTOR_W_GENE_OVERLAP{$newid}) {
            next if ($ra_r->[2] eq 'stop_codon');
            $ra_r->[3] = 1 if ($ra_r->[3] < 1);
        }
        warn "error in: $newid\n" if ($ra_r->[3] < 1 || $ra_r->[4] < 1);
        die "unexpected-3:$ra_r->[3]\n\$offset = $offset" if ($ra_r->[3] < 1);
        die "unexpected-4:$ra_r->[4]\n\$offset = $offset" if ($ra_r->[4] < 1);
        my $gffline = join "\t", @{$ra_r};
        print $gffout "$gffline\n";
    }
}

sub beg_mid_or_end {
    my $ra_c = shift;
    if ($ra_c->[1] == 1) {
        return "beg";
    } elsif ($ra_c->[0] == $ra_c->[2]) {
        return "end";
    } else {
        return "mid";
    }
}

sub get_gff3 {
    my $gff3 = shift;
    my %data = ();
    my $fh = FileHandle->new($gff3,'r');
    while (my $line = <$fh>) {
        if ($line =~ m/^#/) {
            chomp $line;
            push @{$data{'comment'}}, $line;
        } else {
            chomp $line;
            my @f = split /\t/, $line;
            push @{$data{$f[0]}}, \@f;
        }
    }
    return \%data;
}

# INFILE is contams.txt which is formatted as follows:
# ID			TIG_LEN	START..END
# Aala.03.000471	810	1..31
#
# returns ref to %contam = (ID => [TIG_LEN, START, END],...)
sub get_contam_data {
    my $file = shift;
    my %contam = ();
    my $fh = FileHandle->new($CONTAM, 'r');
    while (my $line = <$fh>) {
        chomp $line;
        my @f = split /\t/, $line;
        my @g = split /\./, $f[2];
        die "unexpected: > 1 adapter in a seq" if ($contam{$f[0]});
        $contam{$f[0]} = [$f[1],$g[0],$g[2]];
    }
    return \%contam;
}
